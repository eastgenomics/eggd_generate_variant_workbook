import argparse
from pathlib import Path
import re
import sys
from typing import Union

import numpy as np
from openpyxl.styles import Alignment, Border, colors, Font, Side
import pandas as pd


class vcf():
    """
    Functions to handle reading and manipulating vcf data

    Attributes
    ----------
    args : argparse.Namespace
        arguments passed from command line
    dtypes : dict
        common columns present in annotated vcfs and appropriate dtype to apply
    vcfs : list of pd.DataFrame
        list of dataframes read in from self.args.vcfs
    """

    def __init__(self, args) -> None:
        self.args = args
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

        # read in the vcfs
        self.vcfs = [self.read(x) for x in args.vcfs]

        if args.print_columns:
            self.print_columns()

        if args.filter:
            # filters vcfs against passed parameters
            self.validate_filters()
            self.build_filters()
            self.filter()

        if args.exclude or args.include:
            self.drop_columns()

        if args.reorder:
            self.order_columns()

        if args.rename:
            self.rename_columns()

        if args.merge:
            self.merge()

        self.remove_nan()

        print("Finished munging variants from vcf(s)")


    def read(self, vcf) -> pd.DataFrame:
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
        print(f"Reading in vcf {vcf}")

        header, columns = self.parse_header(vcf)
        csq_fields = self.parse_csq_fields(header)

        # read vcf into pandas df
        vcf_df = pd.read_csv(
            vcf, sep='\t', comment='#', names=columns,
            dtype=self.dtypes, compression='infer'
        )

        if self.args.add_name:
            # add sample name from filename as 1st column
            vcf_df.insert(loc=0, column='sampleName', value=sample)

        # split out INFO column values and CSQ
        vcf_df = self.split_info(vcf_df)
        vcf_df = self.split_csq(vcf_df, csq_fields)

        # drop INFO and CSQ as we fully split them out
        vcf_df.drop(['INFO', 'CSQ'], axis=1, inplace=True)

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


    def split_info(self, vcf_df) -> pd.DataFrame:
        """
        Splits out the INFO column of vcf to all separate values

        Parameters
        ----------
        vcf_df : pd.DataFrame
            dataframe of all variants from a vcf

        Returns
        -------
        vcf_df : pd.DataFrame
            dataframe of all variants from a vcf with separated INFO fields
        """
        # get the key value pairs of INFO data
        info_pairs = [
            x.split('CSQ=')[0].split(';') for x in vcf_df['INFO'].tolist()
        ]
        info_pairs = [[x for x in lst if x] for lst in info_pairs]

        # get unique list of keys from key=value pairs in INFO
        info_keys = sorted(list(set([
            x.split('=')[0] if '=' in x else x for pair in info_pairs for x in pair
        ])))
        info_keys = [x for x in info_keys if x]  # can end up with empty string

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

        for col in info_keys:
            # add all info values to main vcf df
            vcf_df[col] = info_df[col]

        return vcf_df


    def split_csq(self, vcf_df, csq_fields) -> pd.DataFrame:
        """
        Split out CSQ string to separate fields to get annotation

        Parameters
        ----------
        vcf_df : pd.DataFrame
            dataframe of all variants from a vcf
        csq_fields : list
            list of CSQ field names parsed from vcf header, used to assign as
            column headings when splitting out CSQ data for each variant

        Returns
        -------
        vcf_df : pd.DataFrame
            dataframe of all variants from a vcf with separated CSQ fields
        """
        vcf_df['CSQ'] = vcf_df['INFO'].apply(lambda x: x.split('CSQ=')[-1])

        # variants with multiple transcript annotation will have duplicate CSQ
        # data that is comma sepparated => expand this to multiple rows, if no
        # ',' present rows will remain unaffacted (i.e. one transcript)
        columns = list(vcf_df.columns)
        columns.remove('CSQ')
        vcf_df = vcf_df.set_index(columns).apply(
            lambda x: x.str.split(',').explode()
        ).reset_index()

        # split each CSQ value to own columns
        vcf_df[csq_fields] = vcf_df.CSQ.str.split('|', expand=True)

        if 'COSMIC' in vcf_df.columns:
            # handle known bug in VEP annotation where it duplicates COSMIC
            vcf_df['COSMIC'] = vcf_df['COSMIC'].apply(
                lambda x: '&'.join(set(x.split('&')))
            )

        return vcf_df


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
        first as this method is called before writing and it reliably allows
        us to identify and replace all NaN values
        """
        for idx, vcf in enumerate(self.vcfs):
            vcf = vcf.astype(str)

            for col in vcf.columns:
                vcf[col] = vcf[col].apply(
                    lambda x: x.replace('nan', '') if x == 'nan' else x
                )

            self.vcfs[idx] = vcf


    def validate_filters(self) -> None:
        """
        Validate filters passed for filtering variants down.

        These come from args.filter and wil be in the format
        <column><operator><value> (i.e. "Consequence!=synonymous"), currently
        supports the operators >, <, >=, <=, == and !=

        Method will check that a valid operator has been passed, and that the
        given column is present in all columns of all vcfs
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


    def build_filters(self) -> None:
        """
        Formats cmd line passed filters as list of lists of each filter
        expression for filtering df of variants.

        This will separate out the passed filters from the format
        "Consequence!=synonymous" -> ["Consequence", "!=", "synonymous"]
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


    def filter(self) -> None:
        """
        Apply filters passed to ech dataframe of variants.

        Filters first checked in self.validate_filters then formatted in
        self.build filters ready to be passed to np.where()

        Currently wrapped in an eval() call which is not ideal but works for
        interpreting the operator passed as a string from the cmd line
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

            print((
                f"Applied the filters {self.filters} to vcf, filtered out "
                f"{len(all_filter_idxs)} rows"
            ))


    def print_columns(self) -> None:
        """
        Simple method to just print the columns from each vcf and exit.

        Useful for identify what columns are present in INFO and CSQ fields
        for using --include, --exclude and --reorder arguments
        """
        for name, vcf in zip(self.args.vcfs, self.vcfs):
            print(f"Columns for {Path(name).name}: ")
            print(f"\n\t{list(vcf.columns)}\n\n")

        sys.exit()


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
            assert [x for x in vcf.columns for x in to_drop], (
                "Column '{x}' specified with --exclude/--include not "
                "present in one or more of the given vcfs. Valid column "
                f"names: {vcf.columns}"
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
        Rename columnns from key value pairs passed from --rename argument

        Raises
        ------
        AssertionError
            Raised when columns specified with --rename do not exist in more
            or more of the vcfs columns
        """
        for idx, vcf in enumerate(self.vcfs):
            # sense check given reorder keys are in the vcfs
            assert [x for x in vcf.columns for x in self.args.rename.keys()], (
                f"Column(s) specified with --rename not present in one or "
                f"more of the given vcfs. Valid column names: {vcf.columns}."
                f"Column names passed to --rename: {self.args.rename.keys()}"
            )
            self.vcfs[idx].rename(columns=dict(self.args.rename.items()))


    def merge(self) -> None:
        """
        Merge all variants into one big dataframe, should be used with
        --add_name argument if provenance of variants in merged dataframe
        is important
        """
        self.vcfs = [pd.concat(self.vcfs).reset_index(drop=True)]


class excel():
    """
    Functions for wrangling variant data into spreadsheet formatting and
    writing output file

    Attributes
    ----------
    args : argparse.Namespace
        arguments passed from command line
    vcfs : list of pd.DataFrame
        list of dataframes formatted to write to file from vcf() methods
    writer : pandas.io.excel._openpyxl.OpenpyxlWriter
        writer object for writing Excel data to file
    workbook : openpyxl.workbook.workbook.Workbook
        openpyxl workbook object for interacting with per-sheet writing and
        formatting of output Excel file

    Outputs
    -------
    {args.output}.xlsx : file
        Excel file with variants written to, name passed from command line or
        inferred from input vcf name if not specified
    """
    def __init__(self, args, vcfs) -> None:
        print(f"Writing to output file: {args.output}")
        self.args = args
        self.vcfs = vcfs
        self.writer = pd.ExcelWriter(args.output, engine='openpyxl')
        self.workbook = self.writer.book

        self.write_variants()
        self.write_summary()
        self.set_font()

        self.workbook.save(args.output)


    def write_summary(self) -> None:
        """
        Write summary sheet to excel file
        """
        self.summary = self.workbook.create_sheet('summary')
        self.summary.cell(1, 1).value = "Summary"


    def write_variants(self) -> None:
        """
        Writes all variants from dataframe(s) to sheet(s) specified in
        self.args.sheets.

        If sheet names not specified, these will be set as "variant" where one
        dataframe is being written, or the vcf filename prefix if there are
        more than one dataframes to write
        """
        total_rows = sum([len(x) for x in self.vcfs])
        print(f"Writing {total_rows} rows to output xlsx file")
        with self.writer:
            # add variants
            for sheet, vcf in zip(self.args.sheets, self.vcfs):
                vcf.to_excel(
                    self.writer, sheet_name=sheet, index=False
                )


    def set_font(self) -> None:
        """
        Set font to all cells in sheet to Calibri

        Default is Times New Roman and it is ugly
        """
        for ws in self.workbook:
            for cells in ws.rows:
                for cell in cells:
                    cell.font = Font(name="Calibri")



class parsePairs(argparse.Action):
    """
    Simple class method for enabling passing of key value pairs to argparse
    as this is not natively supported by argparse

    Used for passing in names to rename columns in output Excel
    """
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, dict())
        for value in values:
            key, value = value.split('=')
            getattr(namespace, self.dest)[key] = value


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
        help='Annotated vcf file'
    )
    parser.add_argument(
        '-e', '--exclude', nargs='+',
        help=(
            'Columns in vcf to EXCLUDE from output, by default all INFO and '
            'CSQ fields are expanded to their own columns'
        )
    )
    parser.add_argument(
        '-i', '--include', nargs='+',
        help=(
            'Columns in vcf to INCLUDE from output, by default all INFO and '
            'CSQ fields are expanded to their own columns'
        )
    )
    parser.add_argument(
        '-r', '--reorder', required=False, nargs='+',
        help=(
            'Set order for columns in output vcf, any not specified will be '
            'appended to the end'
        )
    )
    parser.add_argument(
        '-z', '--rename', nargs='*', action=parsePairs,
        help=(
            'Pass pairs of {column_name}={new_column_name} for renaming '
            'columns in output excel, should be passed as '
            '--rename CHROM=chr POS=pos REF=ref...'
        )
    )
    parser.add_argument(
        '-f', '--filter', nargs='+',
        help=(
            'Columns on which to filter out variants. Format should be '
            'as <column><operator><value> (gnomAD_AF<0.02). Supported '
            'operands are >|<|>=|<=|==|!='
        )
    )
    parser.add_argument(
        '-k', '--keep', action='store_true',
        help=(
            'Pass when using --filter to keep filtered variants in a '
            'separated "filtered" tab'
        )
    )
    parser.add_argument(
        '-n', '--add_name', action='store_true',
        help='Add sample name from filename as first column'
    )
    parser.add_argument(
        '-s', '--sheets', nargs='+',
        help=(
            'Names to use for multiple sheets, these MUST be the same number '
            'as the number of vcfs passed and in the same order. If not '
            'given, if there is 1 vcf passed the sheet will be named '
            '"variants", if multiple the name prefix of the vcf will be used'
        )
    )
    parser.add_argument(
        '-o', '--output', required=False,
        help=(
            'output name prefix for file, if more than 1 vcf passed a name '
            'must be specified. If only 1 vcf passed and no output name, the '
            'vcf filename prefix will be used'
        )
    )
    parser.add_argument(
        '-m', '--merge', action='store_true',
        help='Merge multiple vcfs into one dataframe of variants to display'
    )
    parser.add_argument(
        '-a', '--analysis', required=False,
        help='Name of analysis to display in summary'
    )
    parser.add_argument(
        '-w', '--workflow', required=False,
        help='ID of workflow to display in summary'
    )
    parser.add_argument(
        '-p', '--panel', required=False,
        help='panel name to display in summary'
    )
    parser.add_argument(
        '--print-columns', required=False, action='store_true',
        help=(
            'Print total columns of all vcfs that will be output to the xlsx. '
            'Useful to identify what will be in the output to include / exclude'
        )
    )

    args = parser.parse_args()

    if not args.output:
        if len(args.vcfs) > 1:
            raise RuntimeError((
                "More than one vcf passed but no output name specified with "
                "--output"
            ))
        else:
            args.output = Path(args.vcfs[0]).name.replace(
                '.vcf', '').replace('.gz', '')

    args.output += ".xlsx"

    if not args.sheets:
        if len(args.vcfs) > 1 and not args.merge:
            # sheet names not specified for > 1 vcf passed => use vcf names
            args.sheets = [Path(x).name.split('_')[0] for x in args.vcfs]
        else:
            # one vcf => name it variants
            args.sheets = ["variants"]
    else:
        assert len(args.vcfs) == len(args.sheets), (
            "Different number of sheets specified to total vcfs passed. Number of "
            f"vcf passed: {len(args.vcfs)}. Number of sheet names passed: "
            f"{len(args.sheets)}"
        )

    return args


def main():
    args = parse_args()

    vcf_handler = vcf(args)
    excel(args, vcf_handler.vcfs)


if __name__ == "__main__":
    main()
