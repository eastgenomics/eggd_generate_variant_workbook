import fileinput
import os
import re
import subprocess
import sys

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

        # add string EXLCUDE to sites that don't meet given filter
        self.args.filter = self.args.filter.replace(
            'filter', 'filter --soft-filter \"EXCLUDE\"'
        )

        print(self.args.filter)
        # sys.exit()

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

        df = pd.read_csv(
            filter_vcf, sep='\t', comment='#', names=columns, compression='infer'
        )

        print(f"\nAll variants:\n")
        print(df)

        include_df = df[df['FILTER'] != 'EXCLUDE']
        exclude_df = df[df['FILTER'] == 'EXCLUDE']

        include_df = splitColumns().split(include_df)
        exclude_df = splitColumns().split(exclude_df)

  

        sys.exit()


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

        Raises
        ------
        AssertionError
            Raised when non zero exit code returned from subprocess.run(),
            either when indexing with tabix or intersecting with bcftools isec

        Returns
        -------
        filtered_out_df : pd.DataFrame
            dataframe of filtered out variants
        """
        # os.makedirs('tmp', exist_ok=True)

        # # index both vcfs with tabix
        # split_output = subprocess.run(
        #     f"tabix -f {split_vcf}", shell=True, capture_output=True)
        # filter_output = subprocess.run(
        #     f"tabix -f {filter_vcf}", shell=True, capture_output=True)

        # assert split_output.returncode == 0 and filter_output.returncode == 0, (
        #     f"\nError in indexing VCF(s)\nExit code for {split_vcf}: "
        #     f"{split_output.returncode}\nExit code for {filter_vcf}: "
        #     f"{filter_output.returncode}.\nstderr:\n{split_output.stderr.decode()}"
        #     f"\n{filter_output.stderr.decode()}"
        # )

        # # use bcftools isec to find excluded variants from bcftools filter
        # isec_command = f"bcftools isec -p tmp {split_vcf} {filter_vcf}"
        # isec_output = subprocess.run(isec_command, shell=True, capture_output=True)

        # assert isec_output.returncode == 0, (
        #     f"\nError in bcftools isec\nReturncode: {isec_output.returncode}"
        #     f"\n{isec_output.stderr.decode()}"
        # )

        # # variants excluded will be in the 0000.vcf
        # filtered_out_df = pd.read_csv(
        #     'tmp/0000.vcf', sep='\t', comment='#', names=columns,
        #     compression='infer'
        # )

        # if self.args.add_name:
        #     # add sample name as first column
        #     sample = split_vcf.replace('.vcf', '').replace('.gz', '')
        #     if '_' in sample:
        #         sample = sample.split('_')[0]

        #     filtered_out_df.insert(loc=0, column='sampleName', value=sample)

        # # tidy up bcftools isec output
        # output_files = [
        #     'tmp/0000.vcf', 'tmp/0001.vcf', 'tmp/0002.vcf', 'tmp/0003.vcf',
        #     'tmp/README.txt', 'tmp/sites.txt'
        # ]
        # # for file in output_files:
        # #     os.remove(file)

        all_variants = pd.read_csv(
            split_vcf, sep='\t', comment='#', names=columns, compression='infer'
        )

        include_variants = pd.read_csv(
            filter_vcf, sep='\t', comment='#', names=columns, compression='infer'
        )
        
        include_variants = splitColumns().split(include_variants)
        all_variants = splitColumns().split(all_variants)
        print(include_variants.columns)
        # sys.exit()

        # print(include_variants['BaseQRankSum'].unique().tolist())

        # for col in include_variants.columns.tolist():
        #     print(col)
        #     print(include_variants[col].unique().tolist()[:3])
        #     print(all_variants[col].unique().tolist()[:3])
        #     print(' ')
        # sys.exit()

        # print(split_vcf)
        # print(filter_vcf)
        # with open(split_vcf, 'r') as fh:
        #     all_variants = fh.readlines()
        #     all_variants = [x for x in all_variants if not x.startswith('#')]
        
        # with open(filter_vcf) as fh:
        #     include_variants = fh.readlines()
        #     include_variants = [x for x in include_variants if not x.startswith('#')]
        
        print(len(all_variants))
        print(len(include_variants))

        # print(type(all_variants[0]))
        # print(include_variants[0])

        # v = all_variants.copy()

        # for x in include_variants:
        #     if x in v:
        #         v.remove(x)

        # print(len(v))
        # sys.exit()


        # print(all_variants)
        # print(include_variants)

        # columns = all_variants.columns.tolist()
        # columns.remove('FILTER')
        # print(columns)

        # all_variants.drop('FILTER',1,inplace=True)
        # include_variants.drop('FILTER',1,inplace=True)
        df = pd.merge(
            all_variants, include_variants, how='outer', suffixes=(None, '_RIGHT'),
            on=['CHROM', 'POS', 'REF', 'ALT', 'CSQ_Consequence', 'CSQ_Feature'], indicator='exist'
        )
        print(df)
        df = df.loc[df['exist'] == 'left_only']

        print(df)

        for column in df.columns:
            if column.endswith('_RIGHT'):
                df.drop(column, 1, inplace=True)
        df.drop('exist', 1, inplace=True)
        
        print(df)
        sys.exit()
        # all_variants = pd.concat([all_variants, include_variants])

        # all_variants.drop('FILTER', 1, inplace=True)

        # # print(all_variants)

        # all_variants = all_variants.drop_duplicates(keep=False)
        # print(all_variants)

        # for x in df[df['POS'] == 23408923]['INFO'].tolist():
        #     print(x)
        #     print(' ')
        
        # print('here')
        # print(' ')
        # for x in include_variants[include_variants['POS'] == 23408923]['INFO']:
        #     print(x)
        sys.exit()

        """
        filter - 383
        all - 23236
        out - 22853
        """


        return filtered_out_df


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
