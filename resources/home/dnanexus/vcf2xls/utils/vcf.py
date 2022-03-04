from pathlib import Path
import re
import sys
from typing import Union
import urllib.parse
from weakref import ref

import numpy as np
import pandas as pd

from .columns import splitColumns
from .filters import filter


class vcf():
    """
    Functions to handle reading and manipulating vcf data

    Called in the order:

    read() -> filter() -> drop_columns() -> reorder() -> merge() -> rename()

    Attributes
    ----------
    args : argparse.Namespace
        arguments passed from command line
    vcfs : list
        list of dataframe(s) of vcf data read in and formatted
    refs : list
        list of genome reference files used for given VCFs
    csq_fields : list
        list of VEP CSQ fields parsed from header
    total_vcf_rows : int
        value to keep total number of rows read in to ensure we don't drop
        any unless --filter is used and --keep is not and => intentionally
        dropping filtered rows
    expanded_vcf_rows : int
        value to track total rows expanded out when multiple transcript
        annotations for one variant are present, resulting in one row per
        transcript annotation per variant in resultant dataframe
    dtypes : dict
        common columns present in annotated vcfs and appropriate dtype to apply
    vcfs : list of pd.DataFrame
        list of dataframes read in from self.args.vcfs
    filtered_rows : pd.DataFrame
        dataframe of all rows dropped from all vcfs
    """

    def __init__(self, args) -> None:
        self.args = args
        self.vcfs = []
        self.refs = []
        self.csq_fields = []
        self.total_vcf_rows = 0
        self.expanded_vcf_rows = 0
        self.filtered_rows = pd.DataFrame()
        self.dtypes = {
            "CHROM": str,
            "POS": int,
            "ID": str,
            "REF": str,
            "ALT": str,
            "QUAL": str,
            "FILTER": str,
            "FORMAT": str,
        }


    def process(self) -> None:
        """
        Function to call all methods in vcf() for processing given VCFs and
        formatting ready to write to output file

        Calls methods in following order:

            - self.read()
            - splitColumns.info()
            - splitColumns.csq()
            - splitColumns.format_fields()
            - self.filter()
            - self.drop_columns()
            - self.reorder()
            - self.merge()
            - self.rename()
            - self.remove_nan()
        """
        # read in the each vcf and apply formatting
        for vcf in self.args.vcfs:
            vcf_df, csq_fields = self.read(vcf)

            # split out INFO column, CSQ fields and FORMAT/SAMPLE column values
            # to individual columns in dataframe
            vcf_df, expanded_vcf_rows = splitColumns.csq(vcf_df, csq_fields)
            vcf_df = splitColumns.info(vcf_df)
            vcf_df = splitColumns.format_fields(vcf_df)

            if 'COSMIC' in vcf_df.columns:
                # handle known bug in VEP annotation where it duplicates COSMIC
                vcf_df['COSMIC'] = vcf_df['COSMIC'].apply(
                    lambda x: '&'.join(set(x.split('&')))
                )

            # TO REMOVE, JUST FOR TESTING SINCE WE HAVE MULTIPLE TRANSCRIPT ANNOTATIONS
            vcf_df = vcf_df[~vcf_df.duplicated(subset=['CHROM', 'POS', 'ID', 'REF', 'ALT'], keep='first')]
            vcf_df = vcf_df.reset_index()

            self.expanded_vcf_rows += expanded_vcf_rows

            # set correct dtypes, required for setting numeric & object types
            # to ensure correct filtering
            vcf_df = self.set_types(vcf_df)

            self.vcfs.append(vcf_df)

        print((
            f"\nTotal variants from {len(self.vcfs)} "
            f"vcf(s): {self.total_vcf_rows}\n"
        ))
        if self.expanded_vcf_rows > 0:
            print(f"Total rows expanded from vcfs: {self.expanded_vcf_rows}")

        if self.args.print_columns:
            self.print_columns()

        if self.args.filter:
            # filter dfs of variants against cmd line specified filters
            filters = filter(self.args, self.vcfs)

            filters.validate()
            filters.build()
            filters.filter()

            # update class with any changes from filtering
            self.vcfs, self.args, self.filtered_rows = filters.vcfs, \
                filters.args, filters.filtered_rows

        if self.args.exclude or self.args.include:
            self.drop_columns()

        if self.args.reorder:
            self.order_columns()

        if self.args.merge:
            self.merge()

        self.add_hyperlinks()
        self.rename_columns()
        self.format_strings()

        # run checks to ensure we haven't unintentionally dropped variants
        # self.verify_totals()

        print("\nSUCCESS: Finished munging variants from vcf(s)\n")


    def read(self, vcf) -> pd.DataFrame:
        """
        Reads given vcf into pd.DataFrame and parses header

        Parameters
        ------
        vcf : str
            path to vcf file to use

        Returns
        -------
        vcf_df : pandas.DataFrame
            dataframe of all variants
        csq_fields : list
            VEP CSQ fields read from header
        """
        sample = Path(vcf).stem
        if '_' in vcf:
            sample = sample.split('_')[0]

        print(f"\n\nReading in vcf {vcf}\n")

        header, columns = self.parse_header(vcf)
        csq_fields = self.parse_csq_fields(header)
        self.parse_reference(header)
        self.parse_field_types(header)

        # read vcf into pandas df
        vcf_df = pd.read_csv(
            vcf, sep='\t', comment='#', names=columns,
            dtype=self.dtypes, compression='infer'
        )

        self.total_vcf_rows += len(vcf_df.index)  # update our total count
        print(f"Total rows in current VCF: {len(vcf_df.index)}")
        print(f"Total rows of all vcfs read in: {self.total_vcf_rows}\n")

        if self.args.add_name:
            # add sample name from filename as 1st column
            vcf_df.insert(loc=0, column='sampleName', value=sample)

        return vcf_df, csq_fields


    def parse_header(self, vcf) -> Union[list, list]:
        """
        Read in header lines of given vcf to list, returning the list and the
        vcf column names

        Parameters
        ----------
        vcf : str
            vcf filename to read header from

        Returns
        -------
        header : list
            list of header lines read from vcf
        columns : list
            column names from vcf
        """
        with open(vcf) as fh:
            # read in header of vcf
            header = []
            for line in fh.readlines():
                if line.startswith('#'):
                    header.append(line.rstrip('\n'))
                else:
                    break

        columns = [x.strip('#') for x in header[-1].split()]
        columns[-1] = 'SAMPLE'

        return header, columns


    def parse_reference(self, header) -> None:
        """
        Parse reference file used from VCF header

        Parameters
        ----------
        header : list
            lines of vcf header
        """
        ref = next(
            iter([x for x in header if x.startswith('##reference')]), None
        )
        if ref:
            if ref not in self.refs:
                # add reference file if found and same not already in list
                self.refs.append(Path(ref).name)


    def parse_csq_fields(self, header) -> list:
        """
        Parse out csq field names from vcf header.

        Kind of horrible way to parse the CSQ field names but this is how VEP
        seems to have always stored them (6+ years +) in the header as:
        ['##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence \
        annotations from Ensembl VEP. Format: SYMBOL|VARIANT_CLASS|...

        Parameters
        ----------
        header : list
            list of header lines read from vcf

        Returns
        -------
        csq_fields : list
            list of CSQ field names parsed from vcf header, used to assign as
            column headings when splitting out CSQ data for each variant
        """
        csq_fields = [x for x in header if x.startswith("##INFO=<ID=CSQ")]
        csq_fields = csq_fields[0].split("Format: ")[-1].strip('">').split('|')

        return csq_fields


    def parse_field_types(self, header) -> None:
        """
        Parse the types from the vcf header for INFO and FORMAT fields, used to
        update self.dtypes object used to correctly set types on dataframe(s)

        Parameters
        ----------
        header : list
            list of header lines read from vcf
        """
        dtypes = {
            "Integer": int,
            "Float": float,
            "Flag": bool,
            "String": str,
            "Character": object
        }

        for line in header:
            if line.startswith('##INFO') or line.startswith('##FORMAT'):
                name = re.search(
                    "ID=.+?(?=,)", line).group().replace("ID=", "")
                type = re.search(
                    "Type=.+?(?=,)", line).group().replace("Type=", "")

                self.dtypes[name] = dtypes[type]


    def set_types(self, vcf_df) -> pd.DataFrame:
        """
        Sets appropriate dtypes on given df of variants

        Parameters
        ----------
        vcf_df : pd.DataFrame
            dataframe of all variants from a vcf

        Returns
        -------
        vcf_df : pd.DataFrame
            dataframe of all variants from a vcf with dtypes set
        """
        # # first set any empty strings to np.nan values to not break types
        # vcf_df = vcf_df.replace('', np.nan)

        # # get any AF and AC columns that should be floats
        # int_columns = [
        #     x for x in vcf_df.columns if '_AF' in x or '_AC' in x
        # ]
        # for col in int_columns:
        #     self.dtypes[col] = float

        # filter all dtypes to just those columns in current df
        df_dtypes = {
            k: v for k, v in self.dtypes.items() if k in list(vcf_df.columns)
        }

        vcf_df = vcf_df.astype(df_dtypes, errors='ignore')

        return vcf_df


    def add_hyperlinks(self) -> None:
        """
        Format column value as an Excel hyperlink, currently supports:
            - ClinVar
            - Cosmic
        """
        urls = {
                "clinvar": "https://www.ncbi.nlm.nih.gov/clinvar/variation/",
                "cosmic": "https://cancer.sanger.ac.uk/cosmic/search?q=",
            }

        # some URLs are build specific, infer which to use from build in header
        reference = self.refs[0].lower()
        if 'grch37' in reference or 'hg19' in reference:
            urls.update({
                "gnomad_af": "https://gnomad.broadinstitute.org/region/CHROM-POS?dataset=gnomad_r2_1",
                "gnomadg_af": "https://gnomad.broadinstitute.org/region/CHROM-POS?dataset=gnomad_r2_1"
            })
        elif 'grch38' in reference or 'hg38' in reference:
            urls.update({
                "gnomad_af": "https://gnomad.broadinstitute.org/region/CHROM-POS?dataset=gnomad_r3",
                "gnomadg_af": "https://gnomad.broadinstitute.org/region/CHROM-POS?dataset=gnomad_r3"
            })

        for idx, vcf in enumerate(self.vcfs):
            for col in vcf.columns:
                if urls.get(col.lower(), None):
                    # column has a linked url => add appropriate hyperlink
                    self.vcfs[idx][col] = self.vcfs[idx].apply(
                        lambda x: self.make_hyperlink(
                            col, urls.get(col.lower()), x
                        ), axis=1
                    )


    def make_hyperlink(self, column, url, value):
        """
        Return Excel formatted hyperlink from given url and value

        Parameters
        ----------
        column : string
            current column for adding URL to
        url : string
            URL string to add value(s) to
        value : pd.Series
            current row values to use for formatting of URL

        Returns
        -------
        str
            url string formatted as Excel hyperlink
        """
        if value[column]:
            if 'gnomad' in column.lower():
                # handle gnomad differently as it requires chrom and pos in URL
                # instead of just the value adding to the end
                chrom = str(value.CHROM.replace('chr', ''))
                url = url.replace('CHROM', chrom).replace('POS', str(value.POS))
            else:
                # other URLs with value appended to end
                url = f'{url}{value[column]}'

            return f'=HYPERLINK("{url}", "{value[column]}")'


    def format_strings(self) -> None:
        """
        Fix formatting of string values with different encoding and nans
        """
        for idx, vcf in enumerate(self.vcfs):
            # pass through urllib unqoute and UTF-8 to fix any weird symbols
            vcf = vcf.applymap(
                lambda x: urllib.parse.unquote(x).encode('UTF-8').decode() if type(x) == str else x
            )

            # remove any nans that are strings
            vcf = vcf.applymap(
                lambda x: x.replace('nan', '') if x == 'nan' and type(x) == str else x
            )

            self.vcfs[idx] = vcf


    def print_columns(self) -> None:
        """
        Simple method to just print the columns from each vcf and exit.

        Useful for identify what columns are present in INFO and CSQ fields
        for using --include, --exclude and --reorder arguments
        """
        for name, vcf in zip(self.args.vcfs, self.vcfs):
            print(f"Columns for {Path(name).name}: ")
            print(f"\n\t{list(vcf.columns)}\n\n")

        sys.exit(0)


    def drop_columns(self) -> None:
        """
        If `--exclude` or `--include` passed, drop given columns
        (or inverse of) from vcf data if they exist.

        If `--include` passed will take the given list of columns and drop the
        remaining columns not specified from all dataframes

        If `--exclude` passed will take the given list of columns and drop
        from all dataframes

        Raises
        ------
        AssertionError
            Raised when columns specified with --include / --exclude are not
            present in one or more of the dataframes
        """
        for idx, vcf in enumerate(self.vcfs):
            if self.args.include:
                # include passed => select all columns not specified to drop
                to_drop = list(
                    set(vcf.columns.tolist()) - set(self.args.include)
                )
            else:
                to_drop = self.args.exclude

            # sense check given exclude columns are in the vcfs
            assert all(column in vcf.columns for column in to_drop), (
                "Column(s) specified with --exclude/--include not present in "
                "one or more of the given vcfs. \n\nValid column names: "
                f"{vcf.columns.tolist()}. \n\nColumns specified: {to_drop}"
            )

            self.vcfs[idx].drop(to_drop, axis=1, inplace=True, errors='ignore')


    def order_columns(self) -> None:
        """
        Reorder columns by specified order from `--reorder` argument, any not
        specified will retain original order after reorder columns

        Raises
        ------
        AssertionError
            Raised when columns specified with --reorder are not
            present in one or more of the dataframes
        """
        for idx, vcf in enumerate(self.vcfs):
            vcf_columns = list(vcf.columns)

            # sense check given exclude columns is in the vcfs
            for x in self.args.reorder:
                assert x in vcf_columns, (
                    f"\n\nColumn '{x}' specified with --reorder not "
                    "present in one or more of the given vcfs. Valid column "
                    f"names: \n{vcf.columns}"
                )

            [vcf_columns.remove(x) for x in self.args.reorder]
            column_order = self.args.reorder + vcf_columns

            self.vcfs[idx] = vcf[column_order]


    def rename_columns(self) -> None:
        """
        Rename columnns from key value pairs passed from --rename argument,
        also remove underscores from all names for nicer reading

        Raises
        ------
        AssertionError
            Raised when columns specified with --rename do not exist in more
            or more of the vcfs columns
        """
        for idx, vcf in enumerate(self.vcfs):
            if self.args.rename:
                # sense check given reorder keys are in the vcfs
                assert [
                    x for x in vcf.columns for x in self.args.rename.keys()
                ], (
                    f"Column(s) specified with --rename not present in one or "
                    f"more of the given vcfs. Valid column names: "
                    f"{vcf.columns}. Column names passed to --rename: "
                    f"{self.args.rename.keys()}"
                )
                self.vcfs[idx].rename(
                    columns=dict(self.args.rename.items()), inplace=True
                )

            # remove underscores from all names
            self.vcfs[idx].columns = [
                x.replace('_', ' ') for x in self.vcfs[idx].columns
            ]


    def merge(self) -> None:
        """
        Merge all variants into one big dataframe, should be used with
        --add_name argument if provenance of variants in merged dataframe
        is important
        """
        if not self.args.keep:
            self.vcfs = [pd.concat(self.vcfs).reset_index(drop=True)]
        else:
            # keep filtered df seperate to write to seperate sheet
            self.vcfs[:-1] = [pd.concat(self.vcfs).reset_index(drop=True)]


    def verify_totals(self) -> None:
        """
        Verify total variants in resultant dataframe(s) match what was read in,
        unless --filter has been applied and --keep has not.

        Totals we have tracked to check on:

        - self.total_vcf_rows -> total rows read in from vcf before modifying
        - self.expanded_rows -> total rows expanded out when CSQ contains
            multiple transcripts, and each transcript for each variant is split
            to individual rows
        - self.filtered_rows -> total rows filtered out to separate df by
            filters passed with --filter
        - total_rows_to_write -> total of all rows that will be written to the
            xlsx across all sheets


        Raises
        ------
        AssertionError
            Raised when the total rows we have tracked don't seem to equal what
            is going to be written to file
        """
        print("Verifying total variants being written to file is correct")
        if not self.args.merge:
            # haven't merged to one df => count all
            total_rows_to_write = sum([len(df.index) for df in self.vcfs])
        else:
            total_rows_to_write = len(self.vcfs[0])

        if self.args.filter:
            total_filtered_rows = len(self.filtered_rows)

        # totals we tracked in previous methods
        tracked_total = int(self.total_vcf_rows) + int(self.expanded_vcf_rows)

        print(f"\nTotal variants identified:")
        print(f"\tTotal rows read in from vcf(s): {self.total_vcf_rows}")
        print((
            f"\tTotal rows expanded: "
            f"{self.expanded_vcf_rows}"
        ))
        print(
            f"\tTotal rows filtered with --filter: {len(self.filtered_rows)}"
        )
        if not self.args.keep and self.args.filter:
            print(
                "\t\t--keep not passed, these filtered rows will be dropped"
            )
        elif self.args.keep and self.args.filter:
            print(
                "\t\t--keep passed, filtered variants will be written to file"
            )
        else:
            print("\t\tNo variants excluded as --filter not specified")

        print(f"\tTotal rows going to write to file: {total_rows_to_write}")
        print(f"\tTotal rows we have tracked: {tracked_total}")


        if not self.args.keep:
            # not keeping filtered rows => check that total we would write if
            # they were included is what we expect from reading in & expanding
            total_rows_to_write += int(len(self.filtered_rows))

        assert tracked_total == total_rows_to_write, (
            "Total rows to be written to file don't appear to equal what has "
            "been tracked, this suggests we have dropped some variants..."
        )

