import fileinput
import os
import re
import subprocess
import sys

import pandas as pd


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


    def get_filtered_rows(self, split_vcf, filter_vcf, columns) -> pd.DataFrame():
        """
        Given dataframe of variants passing filter from vcf, return a
        dataframe of the filtered out variants

        Parameters
        ----------
        split_vcf : str
            filename of vcf of all variants
        filter_vcf : str
            filename of vcf of filtered out variants
        columns : list
            column names read from header of vcf

        Returns
        -------
        filtered_out_df : pd.DataFrame
            dataframe of filtered out variants
        """
        os.makedirs('tmp', exist_ok=True)

        # index both vcfs with tabix
        split_output = subprocess.run(
            f"tabix -f {split_vcf}", shell=True, capture_output=True)
        filter_output = subprocess.run(
            f"tabix -f {filter_vcf}", shell=True, capture_output=True)

        assert split_output.returncode == 0 and filter_output.returncode == 0, (
            f"\nError in indexing VCF(s)\nExit code for {split_vcf}: "
            f"{split_output.returncode}\nExit code for {filter_vcf}: "
            f"{filter_output.returncode}.\nstderr:\n{split_output.stderr.decode()}"
            f"\n{filter_output.stderr.decode()}"
        )

        # use bcftools isec to find excluded variants from bcftools filter
        isec_command = f"bcftools isec -p tmp {split_vcf} {filter_vcf}"
        isec_output = subprocess.run(isec_command, shell=True, capture_output=True)

        assert isec_output.returncode == 0, (
            f"\nError in bcftools isec\nReturncode: {isec_output.returncode}"
            f"\n{isec_output.stderr.decode()}"
        )

        # variants excluded will be in the 0000.vcf
        filtered_out_df = pd.read_csv(
            'tmp/0000.vcf', sep='\t', comment='#', names=columns,
            compression='infer'
        )

        if self.args.add_name:
            # add sample name as first column
            sample = split_vcf.replace('.vcf', '').replace('.gz', '')
            if '_' in sample:
                sample = sample.split('_')[0]

            filtered_out_df.insert(loc=0, column='sampleName', value=sample)

        # tidy up bcftools isec output
        os.remove('tmp/0000.vcf')
        os.remove('tmp/0001.vcf')
        os.remove('tmp/0002.vcf')
        os.remove('tmp/0003.vcf')
        os.remove('tmp/README.txt')
        os.remove('tmp/sites.txt')

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
