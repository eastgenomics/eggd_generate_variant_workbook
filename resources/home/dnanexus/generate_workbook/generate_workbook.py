import argparse
import os
from pathlib import Path
import sys

import pandas as pd

from utils.excel import excel
from utils.vcf import vcf


class arguments():
    """
    Functions for handling and parsing of command line arguments
    """
    def __init__(self):
        """
        Parse command line arguments with argparse, do some sense checking
        and formatting of arguments passed
        """
        self.args = self.parse_args()

        # sense checking and setting / formatting of arguments
        self.check_output()
        self.parse_output()
        self.set_sheet_names()
        self.verify_sheets()
        self.check_include_exclude()


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


    def parse_args(self) -> argparse.Namespace:
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
            help='Annotated vcf file(s)'
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
            '-z', '--rename', nargs='*', action=self.parsePairs,
            help=(
                'Pass pairs of {column_name}={new_column_name} for renaming '
                'columns in output excel, should be passed as '
                '--rename CHROM=chr POS=pos REF=ref...'
            )
        )
        parser.add_argument(
            '-f', '--filter', type=str,
            help=(
                'bcftools filter command to filter variants against'
            )
        )
        parser.add_argument(
            '-t', '--types', nargs='*', action=self.parsePairs, help=(
                'Pass pairs of column=type for modifying types in VCF header,'
                'this should be used when the given type in the VCF header is '
                'incorrect and needs correcting to allow filtering (i.e. '
                'AF fields wrongly set as strings. Types of all fields may be '
                'inspected with --print_header.'
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
            '--keep_tmp', action='store_true',
            help=(
                'Pass to keep intermediary split vcf, and filtered vcf if '
                '--filter specified'
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
            '--out_dir', required=False, default=os.getcwd(),
            help="path to where to output report"
        )
        parser.add_argument(
            '-m', '--merge', action='store_true',
            help='Merge multiple vcfs into one dataframe of variants to display'
        )
        parser.add_argument(
            '--summary', required=False,
            help='summary sheet to include, must be one of: dias'
        )
        parser.add_argument(
            '--analysis', default=('', ''), nargs=2,
            help='Name and ID of analysis to display in summary'
        )
        parser.add_argument(
            '--workflow', default=('', ''), nargs=2,
            help='Name and ID of workflow to display in summary'
        )
        parser.add_argument(
            '--panel', default='',
            help='panel name to display in summary'
        )
        parser.add_argument(
            '--clinical_indication', default='',
            help="clinical indication to write into summary sheet"
        )
        parser.add_argument(
            '--sample', default='',
            help='name of sample to display in summary report'
        )
        parser.add_argument(
            '--print_columns', required=False, action='store_true',
            help=(
                'Print total columns of all vcfs that will be output to the xlsx. '
                'Useful to identify what will be in the output to include/exclude'
            )
        )
        parser.add_argument(
            '--print_header', required=False, action='store_true',
            help=(
                'Print header from first given VCF. Useful for identifying '
                'available INFO/FORMAT/CSQ fields and their types.'
                'n.b. CSQ fields are prefixed with "CSQ_", this is stripped '
                'before writing to the Excel file.'
            )
        )

        return parser.parse_args()


    def check_output(self) -> None:
        """
        Check if args.output specified where multiple VCFs passed, if not
        raise RuntimeError

        Raises
        ------
        RuntimeError
            Raised when multiple VCFs passed with no output name specified
        """
        if len(self.args.vcfs) > 1 and not self.args.output:
            raise RuntimeError((
                "More than one vcf passed but no output name "
                "specified with --output"
            ))


    def parse_output(self) -> None:
        """
        Strip potential messy extensions from input VCF name, then set output
        to include outdir for writing output file
        """
        if not self.args.output:
            self.args.output = Path(
                self.args.vcfs[0]).name.replace('.vcf', '').replace('.gz', '')

        self.args.output = (
            f"{Path(self.args.out_dir)}/{self.args.output}.xlsx"
        )


    def check_include_exclude(self) -> None:
        """
        Use of --inlcude / --exclude is mutually exclusive, therefore check
        and raise exception if both are passed

        Raises
        ------
        AssertionError
            Raised when include and exclude arguments are both specified
        """
        assert not (self.args.include and self.args.exclude), (
            "Both --include and --exclude passed, these arguments are "
            "mutually exclusive."
        )


    def verify_sheets(self) -> None:
        """
        Check if total number of sheets matches number of VCFs where
        args.sheets is passed

        Raises
        ------
        AssertionError
            Raised when number of VCFs does not match the number of sheets
            given with args.sheets
        """
        assert len(self.args.vcfs) == len(self.args.sheets), (
                "Different number of sheets specified to total vcfs passed. "
                f"Number of vcf passed: {len(self.args.vcfs)}. Number of "
                f"sheet names passed: {len(self.args.sheets)}"
            )


    def set_sheet_names(self) -> None:
        """
        Sets list of names for naming output sheets in Excel file
        """
        if not self.args.sheets:
            if len(self.args.vcfs) > 1 and not self.args.merge:
                # sheet names not specified for > 1 vcf passed => use vcf names
                self.args.sheets = [
                    Path(x).name.split('_')[0] for x in self.args.vcfs
                ]
            else:
                # one vcf (or merged) => name it variants
                self.args.sheets = ["variants"]


def main():
    parser = arguments()

    # read in and process vcf(s)
    vcf_handler = vcf(parser.args)
    vcf_handler.process()

    # generate output Excel file
    excel_handler = excel(parser.args, vcf_handler.vcfs, vcf_handler.refs)
    excel_handler.generate()


if __name__ == "__main__":
    main()