import os
from pathlib import Path
import subprocess
import sys
from typing import Union
import urllib.parse

import pandas as pd

from .columns import splitColumns
from .filters import filter


class vcf():
    """
    Functions to handle reading and manipulating vcf data

    Attributes
    ----------
    args : argparse.Namespace
        arguments passed from command line
    vcfs : list
        list of dataframe(s) of vcf data read in and formatted
    refs : list
        list of genome reference files used for given VCFs
    total_vcf_rows : int
        value to keep total number of rows read in to ensure we don't drop
        any unless --filter is used and --keep is not and => intentionally
        dropping filtered rows
    expanded_vcf_rows : int
        value to track total rows expanded out when multiple transcript
        annotations for one variant are present, resulting in one row per
        transcript annotation per variant in resultant dataframe
    filtered_rows : list
        list of dataframes of rows dropped from vcfs
    """

    def __init__(self, args) -> None:
        self.args = args
        self.vcfs = []
        self.refs = []
        self.total_vcf_rows = 0
        self.expanded_vcf_rows = 0
        self.filtered_vcfs = []


    def process(self) -> None:
        """
        Function to call all methods in vcf() for processing given VCFs and
        formatting ready to write to output file

        Calls methods in following order:

            - self.filter() (optional with --filter)
            - self.read()
            - splitColumns.info()
            - splitColumns.format_fields()
            - self.drop_columns()
            - self.reorder()
            - self.merge()
            - self.rename()
        """
        # read in the each vcf, optionally filter, and then apply formatting
        for vcf in self.args.vcfs:
            # first split multiple transcript annotation to separate VCF
            # records, and separate CSQ fields to separate INFO fields
            self.bcftools_pre_process(vcf)

            if self.args.filter:
                # filter vcf against specified filters using bcftools
                filter(self.args).filter('decomposed.tmp.vcf')

                # filters.filter() writes temp filtered vcf containing the
                # filtered variants to read into df
                keep_df = self.read('filtered.tmp.vcf', Path(vcf).stem)
                os.remove('filtered.tmp.vcf')
                _, columns = self.parse_header(vcf)

                # get filtered out rows and read back to new df
                filtered_df = filter(self.args).get_filtered_rows(
                    'decomposed.tmp.vcf', keep_df, columns
                )

                # split out INFO and FORMAT column values to individual
                # columns in dataframe
                keep_df = splitColumns().split(keep_df)
                filtered_df = splitColumns().split(filtered_df)

                self.vcfs.append(keep_df)
                self.filtered_vcfs.append(filtered_df)
            else:
                # not filtering vcf, read in full vcf and split out INFO and
                # FORMAT/SAMPLE column values to individual columns in df
                vcf_df = self.read('decomposed.tmp.vcf', Path(vcf).stem)
                vcf_df = splitColumns().split(vcf_df)
                self.vcfs.append(vcf_df)

            # delete tmp vcf from splitting CSQ str in bcftools_pre_process()
            os.remove('decomposed.tmp.vcf')

        if self.args.print_columns:
            self.print_columns()

        if self.args.merge:
            self.vcfs = self.merge(self.vcfs)

        if self.args.filter and self.args.keep:
            # merge all filtered dataframes to one and add to list of vcfs for
            # doing colum operations and writing to Excel file
            self.filtered_vcfs = self.merge(self.filtered_vcfs)
            self.vcfs.append(self.filtered_vcfs[0])
            self.args.sheets.append('filtered')

        if self.args.exclude or self.args.include:
            self.drop_columns()

        if self.args.reorder:
            self.order_columns()

        self.add_hyperlinks()
        self.rename_columns()
        self.format_strings()
        self.strip_csq_prefix()

        # run checks to ensure we haven't unintentionally dropped variants
        # self.verify_totals()

        print("\nSUCCESS: Finished munging variants from vcf(s)\n")


    def bcftools_pre_process(self, vcf):
        """
        Decompose multiple transcript annotation to individual records, and
        split VEP CSQ string to individual INFO keys. Adds a 'CSQ_' prefix
        to each field extracted from the CSQ string to stop potential conflicts
        with existing INFO fields, which is then stripped before writing
        to the Excel file

        Parameters
        ------
        vcf : str
            path to vcf file to use

        Outputs
        -------
        decomposed.tmp.vcf : file
            vcf file output from bcftools

        Raises
        ------
        AssertionError
            Raised when non-zero exit code returned by bcftools
        """
        cmd = (
            f"bcftools +split-vep --columns - -a CSQ -p 'CSQ_' -d {vcf} | "
            f"bcftools annotate -x INFO/CSQ -o decomposed.tmp.vcf"
        )

        decomposed_vcf = subprocess.run(
            cmd, shell=True,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )

        assert decomposed_vcf.returncode == 0, (
            f"Error in filtering VCF with bcftools. VCF: {vcf}. "
            f"Exitcode:{decomposed_vcf.returncode}."
            f"bcftools filter command used: {self.args.filter}"
        )


    def read(self, vcf, sample) -> pd.DataFrame:
        """
        Reads given vcf into pd.DataFrame and parses header

        Parameters
        ------
        vcf : str
            path to vcf file to use
        sample : str
            name of vcf, used for adding name to df if --add_name passed

        Returns
        -------
        vcf_df : pandas.DataFrame
            dataframe of all variants
        """
        print(f"\n\nReading in vcf {vcf} for sample {sample}\n")

        if '_' in sample:
            sample = sample.split('_')[0]

        header, columns = self.parse_header(vcf)
        self.parse_reference(header)

        if self.args.print_header:
            self.print_header(header)

        # read vcf into pandas df
        vcf_df = pd.read_csv(
            vcf, sep='\t', comment='#', names=columns, compression='infer'
        )

        self.total_vcf_rows += len(vcf_df.index)  # update our total count
        print(f"Total rows in current VCF: {len(vcf_df.index)}")
        print(f"Total rows of all vcfs read in: {self.total_vcf_rows}\n")

        if self.args.add_name:
            # add sample name from filename as 1st column
            vcf_df.insert(loc=0, column='sampleName', value=sample)

        return vcf_df


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


    def add_hyperlinks(self) -> None:
        """
        Format column value as an Excel hyperlink, currently supports:
            - ClinVar
            - Cosmic
        """
        urls = {
                "clinvar": "https://www.ncbi.nlm.nih.gov/clinvar/variation/",
                "cosmic": "https://cancer.sanger.ac.uk/cosmic/search?q=",
                "hgmd": "https://my.qiagendigitalinsights.com/bbp/view/hgmd/pro/mut.php?acc=",
                "mastermind_mmid3": "https://mastermind.genomenon.com/detail?mutation="
            }

        # some URLs are build specific, infer which to use from build in header
        reference = self.refs[0].lower()
        if '37' in reference or 'hg19' in reference:
            urls.update({
                "gnomad_af": "https://gnomad.broadinstitute.org/variant/CHROM-POS-REF-ALT?dataset=gnomad_r2_1",
                "gnomadg_af": "https://gnomad.broadinstitute.org/variant/CHROM-POS-REF-ALT?dataset=gnomad_r2_1"
            })
        elif '38' in reference:
            urls.update({
                "gnomad_af": "https://gnomad.broadinstitute.org/variant/CHROM-POS-REF-ALT?dataset=gnomad_r3",
                "gnomadg_af": "https://gnomad.broadinstitute.org/variant/CHROM-POS-REF-ALT?dataset=gnomad_r3"
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
        if not value[column] or pd.isna(value[column]) or value[column] == 'nan':
            # no value to build hyperlink
            return

        if 'gnomad' in column.lower():
            # handle gnomad differently as it requires chrom, pos, ref and
            # alt in URL instead of just the value adding to the end
            chrom = str(value.CHROM.replace('chr', ''))
            url = url.replace('CHROM', chrom)
            url = url.replace('POS', str(value.POS))
            url = url.replace('REF', value.REF)
            url = url.replace('ALT', value.ALT)
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
                lambda x: urllib.parse.unquote(x).encode('UTF-8').decode() \
                    if type(x) == str else x
            )

            # remove any nans that are strings
            vcf = vcf.applymap(
                lambda x: x.replace('nan', '') \
                    if x == 'nan' and type(x) == str else x
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


    def print_header(self, header) -> None:
        """
        Simple method to print vcf header(s) after splitting CSQ string with
        bcftools to show all available fields and their types

        Parameters
        ----------
        header : list
            vcf header as list read from file
        """
        [print(x) for x in header]
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


    def strip_csq_prefix(self) -> None:
        """
        Strip CSQ prefix added by bcftools -split-vep from column names

        Any conflicts in names with already present columns will retain prefix
        """
        for idx, vcf in enumerate(self.vcfs):
            # strip prefix from column name if present and not already a column
            self.vcfs[idx].columns = [
                x.replace('CSQ ', '', 1) if (
                    x.startswith('CSQ ') and x.replace('CSQ ', '') not in vcf.columns
                ) else x for x in vcf.columns
            ]


    def merge(self, vcfs) -> None:
        """
        Merge all variants into one big dataframe, should be used with
        --add_name argument if provenance of variants in merged dataframe
        is important
        """
        return [pd.concat(vcfs).reset_index(drop=True)]


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
