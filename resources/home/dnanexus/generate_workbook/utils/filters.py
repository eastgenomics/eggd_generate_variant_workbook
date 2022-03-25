import fileinput
import re
import subprocess
import sys

import numpy as np
import pandas as pd

pd.options.mode.chained_assignment = None


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


    def filter(self, split_vcf, filter_vcf) -> None:
        """
        Filter given vcf using bcftoolsd

        Parameters
        ----------
        split_vcf : pathlib.PosixPath
            path to vcf file to filter
        filter_vcf : str
            name for output filtered vcf

        Outputs
        -------
        filter_vcf : file
            vcf file filtered with bedtools and specified filters

        Raises
        ------
        AssertionError
            Raised when non-zero exitcode returned from bcftools annotate
        """
        # first check if any types have been specified to modify before filter
        if self.args.types:
            new_header = self.modify_header_types(split_vcf)
            self.write_header(split_vcf, new_header)

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


    def get_filtered_rows(self, vcf, keep_df, columns) -> pd.DataFrame():
        """
        Given dataframe of variants passing filter from vcf, return a
        dataframe of the filtered out variants

        Parameters
        ----------
        vc : str
            name of full vcf passed to .filter() to read all variants from
        keep_df : pd.DataFrame
            dataframe of variants passing filters
        columns : list
            column names read from header of vcf

        Returns
        -------
        filtered_out_df : pd.DataFrame
            dataframe of filtered out variants
        """
        # read full unfiltered vcf into pandas df
        full_df = pd.read_csv(
            vcf, sep='\t', comment='#', names=columns, compression='infer'
        )

        # need to find filtered rows without FILTER column since this is
        # modified with bcftools
        columns.remove('FILTER')

        # store the current dtypes to set later, set both to str for merging
        dtypes = full_df.dtypes.to_dict()
        keep_df = keep_df.astype(str)
        full_df = full_df.astype(str)

        # get the rows only present in original vcf => filtered out rows
        filtered_out_df = full_df.merge(
            keep_df, how='left', on=columns, indicator=True
        )
        filtered_out_df = filtered_out_df.query('_merge == "left_only"')

        # drop unneeded column and rename filter
        filtered_out_df.drop(['_merge', 'FILTER_y'], axis=1, inplace=True)
        filtered_out_df.rename(columns={'FILTER_x': 'FILTER'}, inplace=True)
        filtered_out_df.reset_index(inplace=True, drop=True)

        filtered_out_df = filtered_out_df.astype(dtypes)

        return filtered_out_df


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
