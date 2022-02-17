from pathlib import Path
import sys
from typing import Union
import urllib.parse

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
    dtypes : dict
        common columns present in annotated vcfs and appropriate dtype to apply
    vcfs : list of pd.DataFrame
        list of dataframes read in from self.args.vcfs
    filtered_rows : pd.DataFrame
        dataframe of all rows dropped from all vcfs
    """

    def __init__(self, args) -> None:
        self.args = args
        self.refs = []
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
            "DP": int,
            "SYMBOL": str,
            "VARIANT_CLASS": str,
            "Consequence": str,
            "EXON": str,
            "HGVSc": str,
            "HGVSp": str,
            "gnomAD_AF": float,
            "gnomADg_AF": float,
            "CADD_PHRED": float,
            "Existing_variation": str,
            "ClinVar": str,
            "ClinVar_CLNDN": str,
            "ClinVar_CLNSIG": str,
            "COSMIC": str,
            "Feature": str,
            "STR": bool,
            "RU": bool,
            "Prev_AC": pd.Int16Dtype(),
            "Prev_NS": pd.Int16Dtype()
        }


    def process(self) -> None:
        """
        Function to call all methods in vcf() for processing given VCFs and
        formatting ready to write to output file
        """
        # read in the vcfs
        self.vcfs = [self.read(x) for x in self.args.vcfs]

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

        self.rename_columns()

        self.remove_nan()

        # run checks to ensure we haven't unintentionally dropped variants
        self.verify_totals()

        print("\nSUCCESS: Finished munging variants from vcf(s)\n")


    def read(self, vcf) -> pd.DataFrame:
        """
        Reads given vcf into pd.DataFrame, calls following methods:

        - self.parse_header()
        - splitColumns.parse_csq_fields()
        - splitColumns.info()
        - splitColumns.csq()
        - splitColumns.format_fields()

        Parameters
        ------
        vcf : str
            path to vcf file to use

        Returns
        -------
        vcf_df : pandas.DataFrame
            dataframe of all variants
        header : list
            vcf header lines
        """
        sample = Path(vcf).stem
        if '_' in vcf:
            sample = sample.split('_')[0]

        print(f"\n\nReading in vcf {vcf}\n")

        header, columns = self.parse_header(vcf)
        self.parse_reference(header)

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

        # split out INFO column values and CSQ
        vcf_df['CSQ'] = vcf_df['INFO'].apply(lambda x: x.split('CSQ=')[-1])
        csq_fields = splitColumns.parse_csq_fields(header)
        vcf_df = splitColumns.info(vcf_df)
        vcf_df, self.expanded_vcf_rows = splitColumns.csq(
            vcf_df, csq_fields
        )
        vcf_df = splitColumns.format_fields(vcf_df)

        vcf_df = self.set_types(vcf_df)

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
        # first set any empty strings to pd.NA values to not break types
        vcf_df = vcf_df.replace('', np.nan)

        # get any AF and AC columns that should be floats
        int_columns = [
            x for x in vcf_df.columns if '_AF' in x or '_AC' in x
        ]
        for col in int_columns:
            self.dtypes[col] = float

        # filter all dtypes to just those columns in current df
        df_dtypes = {
            k: v for k, v in self.dtypes.items() if k in list(vcf_df.columns)
        }

        return vcf_df.astype(df_dtypes, errors='ignore')


    def remove_nan(self) -> None:
        """
        Remove NaN values from all dataframes, cast everything to a string
        first as this method is called before writing and it reliably
        allows us to identify and replace all NaN values.
        """
        for idx, vcf in enumerate(self.vcfs):
            vcf = vcf.astype(str)
            # pass through urllib unqoute and UTF-8 to fix any weird symbols
            vcf = vcf.applymap(
                lambda x: urllib.parse.unquote(x).encode('UTF-8').decode()
            )

            for col in vcf.columns:
                vcf[col] = vcf[col].apply(
                    lambda x: x.replace('nan', '') if x == 'nan' else x
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

            # sense check given exclude columns is in the vcfs
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
            assert [x for x in vcf.columns for x in self.args.reorder], (
                "Column '{x}' specified with --reorder not "
                "present in one or more of the given vcfs. Valid column "
                f"names: {vcf.columns}"
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
                self.vcfs[idx].rename(columns=dict(self.args.rename.items()))

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

