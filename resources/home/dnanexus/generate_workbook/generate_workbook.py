import argparse
import os
from pathlib import Path
import sys

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

        print(f"Arguments passed: ", ''.join([
            f"\n\t\t{' : '.join((str(x), str(self.args.__dict__[x])))}"
            for x in self.args.__dict__
        ]))


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


    class joinList(argparse.Action):
        """
        Merge a list input to a single string

        Used where panels and clinical indication may be passed in as a list
        and need to format as a string to display in summary of workbook
        """
        def __call__(self, parser, namespace, values, option_string=None):
            setattr(namespace, self.dest, ' '.join(values))


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
        group = parser.add_mutually_exclusive_group()
        group.add_argument(
            '-e', '--exclude', nargs='+',
            help=(
                'Columns from vcf to EXCLUDE from output, by default all INFO '
                'and CSQ fields are expanded to their own columns'
            )
        )
        group.add_argument(
            '-i', '--include', nargs='+',
            help=(
                'Columns in vcf to INCLUDE from output, by default all INFO '
                'and CSQ fields are expanded to their own columns'
            )
        )
        parser.add_argument(
            '-r', '--reorder', required=False, nargs='+',
            help=(
                'Set order for columns in output workbook, any not specified '
                'will be appended to the end'
            )
        )
        parser.add_argument(
            '-z', '--rename', nargs='*', action=self.parsePairs,
            help=(
                'Pass pairs of {column_name}={new_column_name} for renaming '
                'columns in output excel, should be passed as '
                '--rename CHROM=chr POS=pos REF=ref... (use --print_columns '
                'to see valid column names)'
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
            '-c', '--add_comment_column', action='store_true',
            help='Add empty comment column to end of each sheet'
        )
        parser.add_argument(
            '-s', '--sheets', nargs='+',
            help=(
                'Names to use for multiple sheets, these MUST be the same '
                'number as the number of vcfs passed and in the same order. '
                'If not given, if there is 1 vcf passed the sheet will be '
                'named "variants", if multiple the name prefix of the vcf '
                'will be used'
            )
        )
        parser.add_argument(
            '-o', '--output', required=False,
            help=(
                'output name prefix for file, if more than 1 vcf passed a '
                'name must be specified. If only 1 vcf passed and no output '
                'name, the vcf filename prefix will be used'
            )
        )
        parser.add_argument(
            '--out_dir', required=False, default=os.getcwd(),
            help="path to output directory"
        )
        parser.add_argument(
            '-m', '--merge', action='store_true',
            help=(
                'Merge multiple vcfs into one combined sheet of '
                'variants to display'
            )
        )
        parser.add_argument(
            '--summary', required=False,
            help='summary sheet to include, must be one of: dias'
        )
        parser.add_argument(
            '--human_filter',  nargs='+', action=self.joinList,
            help=(
                'String to add to summary sheet with humanly readable form of '
                'the given filter string. No checking is done of this matching'
                ' the actual filter(s) used.'
            )
        )
        parser.add_argument(
            '--acmg', action='store_true',
            help='add extra ACMG reporting template sheet'
        )
        parser.add_argument(
            '--job_id', required=False,
            help='Job ID of eggd_generate_workbook to add to Dias summary'
        )
        parser.add_argument(
            '--workflow', default=('', ''), nargs=2,
            help='Name and ID of workflow to display in summary'
        )
        parser.add_argument(
            '--panel', nargs='+', action=self.joinList,
            help='panel name to display in summary'
        )
        parser.add_argument(
            '--clinical_indication', nargs='+', action=self.joinList,
            help="clinical indication to write into summary sheet"
        )
        parser.add_argument(
            '--sample', default='',
            help='name of sample to display in summary report'
        )
        parser.add_argument(
            '--print_columns', required=False, action='store_true',
            help=(
                'Print total columns of all vcfs that will be output to the '
                'workbook. Useful to identify what will be in the output to '
                'include/exclude'
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


    def verify_sheets(self) -> None:
        """
        Check if total number of sheets matches number of VCFs

        Raises
        ------
        AssertionError
            Raised when number of VCFs does not match the number of sheets
            given with args.sheets

        AssertionError
            Raised when --merge and --sheets given, but more than one sheet
            name specified
        """
        if not self.args.merge:
            assert len(self.args.vcfs) == len(self.args.sheets), (
                "Different number of sheets specified to total vcfs passed. "
                f"Number of vcf passed: {len(self.args.vcfs)}. Number of "
                f"sheet names passed: {len(self.args.sheets)}"
            )
        else:
            # merging => should be at most one name passed (or will default
            # to 'variants' in set_sheet_names())
            if self.args.sheets:
                assert len(self.args.sheets) == 1, (
                    "--merge specified but more than one user defined sheet "
                    "name passed. Either pass one name, or do not pass the "
                    "--sheets arg and the default will be used ('variants')"
                )


    def set_sheet_names(self) -> None:
        """
        Sets list of names for naming output sheets in Excel file
        """
        if not self.args.sheets:
            if len(self.args.vcfs) > 1 and not self.args.merge:
                # sheet names not specified for > 1 vcf passed => use vcf names
                self.args.sheets = [
                    Path(x).name.split('_')[0] if '_' in x else
                    Path(x).stem.replace('.vcf', '') for x in self.args.vcfs
                ]
            elif self.args.filter:
                # one vcf (or merged) and filtering
                self.args.sheets = ["included"]
            else:
                # one vcf (or merged) and NOT filtering => name it variants
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
