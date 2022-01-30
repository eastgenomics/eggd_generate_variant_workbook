import argparse
from collections import OrderedDict
from distutils.log import info
import gzip
import os
import re
import sys

import numpy as np
import pandas as pd



class vcf():
    """
    Functions to handle reading and manipulating vcf data
    """

    def __init__(self, args) -> None:
        # returns list of variants from each vcf
        self.args = args
        self.vcfs = [self.read(x) for x in args.vcfs]

        if args.filter:
            # filters vcfs against passed parameters
            self.validate_filters()
            self.build_filters()
            self.filter()


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

        # kind horrible way to parse the CSQ fields but this is how VEP
        # seems to have always stored them (6+ years +) in the header as:
        # ['##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence \
        # annotations from Ensembl VEP. Format: SYMBOL|VARIANT_CLASS|...
        csq_fields = [x for x in header if x.startswith("##INFO=<ID=CSQ")]
        csq_fields = csq_fields[0].split("Format: ")[-1].strip('">').split('|')

        # read vcf into pandas df
        vcf_df = pd.read_csv(
            vcf, sep='\t', comment='#', names=columns
        )

        # split out INFO into respective columns, everything will be in key=value
        # pairs except for the CSQ fields from above. Therefore, we can split those
        # out to separate columns, then split the rest of pairs to own columns
        csq_values = [
            x.split('CSQ=')[-1].split('|') for x in vcf_df['INFO'].tolist()
        ]

        for idx, name in enumerate(csq_fields):
            # add each column of data to df by index
            vcf_df[name] = [value[idx] for value in csq_values]

        # get the key value pairs of INFO data
        info_pairs = [
            x.split('CSQ=')[0].split(';') for x in vcf_df['INFO'].tolist()
        ]
        info_pairs = [[x for x in lst if x] for lst in info_pairs]

        # get unique list of keys from key=value pairs in INFO
        info_keys = sorted(list(set([
            x.split('=')[0] if '=' in x else x for lst in info_pairs for x in lst
        ])))
        info_keys = [x for x in info_keys if x]

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

        vcf_df.drop('INFO', axis=1, inplace=True)  # drop INFO as we fully split it out
        pd.set_option('display.max_rows', None)  # or 1000

        vcf_df['Prev_AC'] = vcf_df['Prev_AC'].apply(lambda x: 0 if x == '' else x)

        print(vcf_df['Prev_AC'])

        vcf_df = vcf_df.astype({'Prev_AC': int})


        return vcf_df


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

        for vcf in self.vcfs:
            for filter in self.filters:
                
                filter_idxs = eval((
                    f"vcf['{filter[0]}'].apply("
                    f"lambda x: x {filter[1]} {int(filter[2])})"
                ))


        



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
        '-e', '--exclude', nargs='?',
        help=(
            'columns in vcf to exclude from output, by default all INFO and '
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
    vcf(args)

    print(args)

    sys.exit()

    header, vcf_df = read_vcf(args.vcf)

    print(vcf_df.columns)


if __name__ == "__main__":
    main()
