import fileinput
import os
import re
import subprocess
import sys
from typing import Union

import pandas as pd

from utils.columns import splitColumns


class filter():
    """"
    Functions related to filtering dataframes of variants from args.filter,
    called during processing of variants in vcf.process()

    Attributes
    ----------
    args : argparse.Namespace
        arguments passed from command line
    vcfs : list of pd.DataFrame
        list of dataframes read in from self.args.vcfs
    filtered_rows : pd.DataFrame
        dataframe of all rows dropped from all vcfs
    """
    def __init__(self, args) -> None:
        self.args = args


    def filter(self, split_vcf, filter_vcf, columns) -> None:
        """
        Filter given vcf using bcftools

        Parameters
        ----------
        split_vcf : pathlib.PosixPath
            path to vcf file to filter
        filter_vcf : str
            name for output filtered vcf

        Outputs
        -------
        filter_vcf : file
            vcf file with PASS/EXCLUDE added to FILTER columns

        Raises
        ------
        AssertionError
            Raised when non-zero exitcode returned from bcftools annotate
        """
        # first check if any types have been specified to modify before filter
        if self.args.types:
            new_header = self.modify_header_types(split_vcf)
            self.write_header(split_vcf, new_header)

        # add string EXLCUDE to sites that don't meet given filter
        self.args.filter = self.args.filter.replace(
            'filter', 'filter --soft-filter \"EXCLUDE\" -m +'
        )

        # write to temporary vcf files to read from with vcf.read()
        command = f"{self.args.filter} {split_vcf} -o {filter_vcf}"

        print(
            f"\nFiltering {split_vcf} with the command: \n\t{command}\n"
        )

        output = subprocess.run(command, shell=True, capture_output=True)

        assert output.returncode == 0, (
            f"\n\tError in filtering VCF with bcftools\n"
            f"\n\tVCF: {split_vcf}\n"
            f"\n\tExitcode:{output.returncode}\n"
            f"\n\tbcftools filter command used: {self.args.filter}\n"
            f"\n\t{output.stderr.decode()}"
        )


    def split_include_exclude(self, variant_df) -> Union[pd.DataFrame, pd.DataFrame]:
        """
        Split out given dataframe of variants by FILTER column when filtered
        through self.filter().

        self.filter() applies a soft filter with bcftools and appends the term
        EXCLUDE to the end of the FILTER column fields, this is used to
        identify which variants met the given filter and should be excluded.

        Parameters
        ----------
        variant_df : pd.DataFrame
            dataframe of variants read from vcf output by self.filter()

        Returns
        -------
        include_df : pd.DataFrame
            dataframe of variants passing specified filter(s)
        exclude_df : pd.DataFrame
            dataframe of variants not passing specified filter(s)
        """
        include_df = variant_df[~variant_df['FILTER'].str.endswith(
            'EXCLUDE')].reset_index(drop=True)
        exclude_df = variant_df[variant_df['FILTER'].str.endswith(
            'EXCLUDE')].reset_index(drop=True)

        # remove our added soft filter tag from FILTER column
        exclude_df['FILTER'] = exclude_df['FILTER'].apply(
            lambda x: x.replace('EXCLUDE', '')
        )

        return include_df, exclude_df



    def verify_total_variants(self, split_vcf, include_df, exclude_df) -> None:
        """
        Verify no variants are dropped from filtering by checking total
        included and excluded dataframe rows against input VCF

        Parameters
        ----------
        split_vcf : string
            filename of vcf used for filtering
        include_df : pd.DataFrame
            dataframe of variants retained from bcftools filter
        exclude_df : pd.DataFrame
            dataframe of variants excluded by bcftools filter

        Raises
        ------
        AssertionError
            Raised when total variants in the include and exclude dataframe
            do not equal the total variants in the input vcf
        """
        print(f"\nVerifying total variants after filtering\n")
        output = subprocess.run(
            f"zgrep -v '^#' {split_vcf} | wc -l", shell=True,
            capture_output=True
        )

        assert output.returncode == 0, (
            f"\n\tError in reading total rows from vcf used for filtering\n"
            f"\n\tExit code: {output.returncode}\n"
            f"\n\tCommand: zgrep -v '^#' {split_vcf} | wc -l\n"
            f"\n\t{output.stderr.decode()}"
        )

        vcf_total = int(output.stdout.decode())

        print(
            f"\tTotal variants in input vcf: {vcf_total}\n"
            f"\tTotal variants included: {len(include_df)}\n"
            f"\tTotal variants excluded: {len(exclude_df)}\n\n"
        )

        assert vcf_total == len(include_df) + len(exclude_df), (
            "Total variants in input VCF does not match what is included + "
            "excluded"
        )


    def modify_header_types(self, vcf) -> list:
        """
        Reads header from vcf and returns header as a list with specified
        column types modified

        Types for INFO fields can be wrongly set in VEP, and we allow for
        column=type pairs to be passed via --types to change these to allow
        for the correct type to be set before filtering. Most commonly these
        are fields wrongly set as strings

        Parameters
        ----------
        vcf : string
            vcf file to read header from

        Returns
        -------
        new_header : list
            lines of header with modified types
        """
        new_header = []

        with open(vcf) as fh:
            for line in fh:
                if not line.startswith('#'):
                    # end of header
                    break
                for column, type in self.args.types.items():
                    if line.startswith(f'##INFO=<ID={column},'):
                        # correct line for given column, find type and replace
                        curr_type = re.search(
                            'Type=(Integer|Float|Flag|Character|String)', line
                        )
                        if curr_type:
                            # found valid type in INFO string to replace
                            line = line.replace(
                                curr_type.group(), f"Type={type.capitalize()}"
                            )
                new_header.append(line.rstrip('\n'))

        return new_header


    def write_header(self, vcf, new_header) -> None:
        """
        Write the new header with modified types to the temporary vcf

        Parameters
        ----------
        vcf : string
            vcf file to modify
        new_header : list
            lines of header with modified types
        """
        with fileinput.FileInput(vcf, inplace=True) as fh:
            for idx, line in enumerate(fh):
                if idx < len(new_header):
                    sys.stdout.write(f"{new_header[idx]}\n")
                else:
                    # end of header
                    sys.stdout.write(line)
