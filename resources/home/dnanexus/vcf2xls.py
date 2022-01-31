import argparse
from pathlib import Path
import re
import sys

import numpy as np
from openpyxl.styles import Alignment, Border, colors, Font, Side
import pandas as pd
from pygments import highlight



class vcf():
    """
    Functions to handle reading and manipulating vcf data
    """

    def __init__(self, args) -> None:
        self.args = args
        # possible columns and appropriate dtypes
        self.dtypes = {
            "CHROM": str,
            "POS": int,
            "REF": str,
            "ALT": str,
            "DP": int
        }
        # read in the vcfs
        self.vcfs = [self.read(x) for x in args.vcfs]

        if args.filter:
            # filters vcfs against passed parameters
            self.validate_filters()
            self.build_filters()
            self.filter()

        if args.exclude or args.include:
            self.drop_columns()

        if args.merge:
            self.merge()


    def read(self, vcf):
        """
        Reads given vcf into pd.DataFrame

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
        sample = Path(vcf).stem.split('_')[0]

        with open(vcf) as fh:
            # read in header of vcf
            header = []
            for line in fh.readlines():
                if line.startswith('#'):
                    header.append(line.rstrip('\n'))
                else:
                    break

        # split last line of header out header to use as column names
        columns = [x.strip('#') for x in header[-1].split()]

        # kind of horrible way to parse the CSQ field names but this is how VEP
        # seems to have always stored them (6+ years +) in the header as:
        # ['##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence \
        # annotations from Ensembl VEP. Format: SYMBOL|VARIANT_CLASS|...
        csq_fields = [x for x in header if x.startswith("##INFO=<ID=CSQ")]
        csq_fields = csq_fields[0].split("Format: ")[-1].strip('">').split('|')

        # read vcf into pandas df
        vcf_df = pd.read_csv(
            vcf, sep='\t', comment='#', names=columns
        )

        if self.args.add_name:
            # add sample name from filename as 1st column
            vcf_df.insert(loc=0, column='sampleName', value=sample)

        print(vcf_df)

        info_df, info_keys = self.split_info(vcf_df, csq_fields)

        for col in info_keys:
            # add all info values to main vcf df
            vcf_df[col] = info_df[col]

        vcf_df.drop('INFO', axis=1, inplace=True)  # drop INFO as we fully split it out

        vcf_df = self.set_types(vcf_df)

        return vcf_df


    def split_info(self, vcf_df, csq_fields):
        """
        Given a vcf of data read in to a df, split out the INFO column to
        all separate values
        """
        # split CSQ values to own column
        vcf_df['CSQ'] = vcf_df['INFO'].apply(lambda x: x.split('CSQ=')[-1])

        # variants with multiple transcript annotation will have duplicate CSQ
        # data that is comma sepparated => expand this to multiple rows
        columns = list(vcf_df.columns)
        columns.remove('CSQ')
        vcf_df = vcf_df.set_index(columns).apply(
            lambda x: x.str.split(',').explode()
        ).reset_index()

        # split each CSQ value to own columns
        vcf_df[csq_fields] = vcf_df.CSQ.str.split('|', expand=True)

        # split out rest of INFO into respective columns, everything will be
        # in key=value pairs or single values that are flags

        # get the key value pairs of INFO data
        info_pairs = [
            x.split('CSQ=')[0].split(';') for x in vcf_df['INFO'].tolist()
        ]
        info_pairs = [[x for x in lst if x] for lst in info_pairs]

        # get unique list of keys from key=value pairs in INFO
        info_keys = sorted(list(set([
            x.split('=')[0] if '=' in x else x for lst in info_pairs for x in lst
        ])))
        info_keys = [x for x in info_keys if x]  # can end up with empty strings

        info_values = []

        for variant in info_pairs:
            # for every variants values, split them out to dict to add to df
            pair_values = {}

            for value in variant:
                if '=' in value:
                    # key value pair
                    key, value = value.split('=')
                else:
                    # Flag value present (e.g RU, STR)
                    key, value = key, True

                pair_values[key] = value

            info_values.append(pair_values)

        # build df of values to add to main df
        info_df = pd.DataFrame(
            info_values, columns=info_keys
        )

        if 'COSMIC' in vcf_df.columns:
            # handle known bug in VEP annotation where it duplicates COSMIC
            vcf_df['COSMIC'] = vcf_df['COSMIC'].apply(
                lambda x: '&'.join(set(x.split('&')))
            )

        return info_df, info_keys


    def set_types(self, vcf_df):
        """
        Sets appropriate dtypes on given df of variants
        """
        # first set any empty strings to pd.NA values to not break types
        vcf_df = vcf_df.replace('', np.nan)

        # get any AF and AC columns that should be floats
        int_columns = [
            x for x in vcf_df.columns if '_AF' in x or '_AC' in x
        ]
        for col in int_columns:
            self.dtypes[col] = float

        return vcf_df.astype(self.dtypes)


    def validate_filters(self):
        """
        Validate filters passed for filtering variants down
        """
        for filter in self.args.filter:
            # check a valid operand passed
            assert len(re.findall('>|<|>=|<=|==|!=', filter)) == 1, (
                f"invalid operand passed in filter: {filter}"
            )

        for vcf in self.vcfs:
            # for each filter, check the specified column is in all the vcfs
            for filter in self.args.filter:
                assert re.split(r'>|<|>=|<=|==|!=', filter)[0] in vcf.columns, (
                    f"Column specified in filter '{filter}' not in vcf "
                    f"columns: {vcf.columns}"
                )


    def build_filters(self):
        """
        Returns cmd line passed filters as list of lists of each filter
        expression for filtering df of variants.


        """
        field_value = [
            re.split(r'>|<|>=|<=|==|!=', x) for x in self.args.filter
        ]
        operator = [
            re.findall('>|<|>=|<=|==|!=', x) for x in self.args.filter
        ]

        self.filters = [
            [x[0], y[0], x[1]] for x, y in zip(field_value, operator)
        ]


    def filter(self):
        """
        Apply filters passed to dfs of variants
        """
        # build list of indices of variants to filter out against specified
        # filters, then apply filter to df, retaining filtered rows
        self.filtered_rows = pd.DataFrame()

        for idx, vcf in enumerate(self.vcfs):
            all_filter_idxs = []
            for filter in self.filters:
                col, operator, value = filter[0], filter[1], filter[2]
                if pd.api.types.is_numeric_dtype(vcf[f'{filter[0]}']):
                    # check column we're filtering is numeric and set types
                    value = int(value)
                else:
                    # string values have to be wrapped in quotes from np.where()
                    value = f"'{value}'"

                # get row indices to filter
                filter_idxs = eval((
                    f"np.where(vcf['{col}'].apply("
                    f"lambda x: x {operator} {value}))[0]"
                ))
                all_filter_idxs.extend(filter_idxs)

            all_filter_idxs = sorted(all_filter_idxs)

            # apply the filter, assign back the filtered df
            self.filtered_rows = self.filtered_rows.append(
                vcf.loc[all_filter_idxs]
            )
            self.vcfs[idx] = vcf.drop(all_filter_idxs)


    def drop_columns(self):
        """
        If --exclude or --include passed, drop given columns (or inverse of)
        from vcf data if they exist.
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
            assert [x for x in vcf.columns for x in to_drop], (
                "Column '{x}' specified with --exclude not present in one or more "
                "of the given vcfs."
            )
            self.vcfs[idx].drop(to_drop, axis=1, inplace=True)


    def merge(self):
        """
        Merge all variants into one big sorted dataframe
        """
        self.vcfs = pd.concat(self.vcfs).reset_index(drop=True)


class excel():
    """
    Functions for wrangling variant data into spreadsheet formatting and
    writing output file
    """
    def __init__(self, args, vcfs) -> None:
        self.args = args
        self.vcfs = vcfs
        self.writer = pd.ExcelWriter(f"test.xlsx", engine='openpyxl')
        self.workbook = self.writer.book

        self.write_variants()
        self.write_summary()
        self.set_font()


        self.workbook.save(f"test.xlsx")


    def write_summary(self):
        """
        Write summary sheet to excel file
        """
        self.summary = self.workbook.create_sheet('summary')
        self.summary.cell(1, 1).value = "Summary"


    def write_variants(self):
        """

        """
        with self.writer:
            # add variants
            self.vcfs.to_excel(self.writer, sheet_name="variants", index=False)


    def set_font(self):
        """
        Set font to all cells in sheet to Calibri
        """
        for ws in self.workbook:
            for cells in ws.rows:
                for cell in cells:
                    cell.font = Font(name="Calibri")


def parse_args() -> argparse.Namespace:
    """
    Parse command line arguments

    Returns
    -------
    args : Namespace
        Namespace of passed command line argument inputs
    """
    parser = argparse.ArgumentParser(
        'Turns an annotated vcf into an xlsx for human viewing'
    )

    parser.add_argument(
        '-v', '--vcfs', nargs='+',
        help='annotated vcf file'
    )
    parser.add_argument(
        '-e', '--exclude', nargs='+',
        help=(
            'columns in vcf to EXCLUDE from output, by default all INFO and '
            'CSQ fields are expanded to their own columns'
        )
    )
    parser.add_argument(
        '-i', '--include', nargs='+',
        help=(
            'columns in vcf to INCLUDE from output, by default all INFO and '
            'CSQ fields are expanded to their own columns'
        )
    )
    parser.add_argument(
        '-f', '--filter', nargs='+',
        help=(
            'columns on which to filter out variants. Format should be '
            'as <column><operator><value> (gnomAD_AF<0.02). Supported '
            'operands are >|<|>=|<=|==|!='
        )
    )
    parser.add_argument(
        '-k', '--keep', action='store_true',
        help=(
            'pass when using --filter to keep filtered variants in a '
            'separated "filtered" tab'
        )
    )
    parser.add_argument(
        '-s', '--add_name', action='store_true',
        help='Add sample name from filename as first column'
    )
    parser.add_argument(
        '-m', '--merge', action='store_true',
        help='Merge multiple vcfs into one dataframe of variants to display'
    )
    parser.add_argument(
        '-a', '--analysis', required=False,
        help='name of analysis to display in summary'
    )
    parser.add_argument(
        '-w', '--workflow', required=False,
        help='id of workflow to display in summary'
    )
    parser.add_argument(
        '-p', '--panel', required=False,
        help='panel name to display in summary'
    )

    return parser.parse_args()


def main():
    args = parse_args()
    vcf_handler = vcf(args)
    excel_handler = excel(args, vcf_handler.vcfs)


if __name__ == "__main__":
    main()
