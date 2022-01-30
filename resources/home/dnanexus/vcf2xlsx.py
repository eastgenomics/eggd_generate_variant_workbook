import argparse
from collections import OrderedDict
import gzip
import os
import sys

import pandas as pd

VCF_HEADER = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']

def read_vcf(vcf):
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

    print(csq_fields)

    # read vcf into pandas df
    vcf_df = pd.read_csv(
        vcf, sep='\t', comment='#', names=columns
    )

    print(vcf_df)

    # split out INFO into respective columns, everything will be in key=value
    # pairs except for the CSQ fields from above. Therefore, we can split those
    # out to separate columns, then split the rest of pairs to own columns

    # separate CSQ value to list of lists, then add to df with field names
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



    # get unique list of keys from key=value pairs in INFO, omit if they have
    # no value present
    info_keys = sorted(list(set([
        x.split('=')[0] for l in info_pairs for x in l if '=' in x
    ])))

    print(info_keys)

    print(info_pairs[1])

    # sys.exit()

    info_values = [
        [x.split('=')[1] for x in sub_list if '=' in x and x] for sub_list in info_pairs
    ]

    for x in info_pairs:
        print(x)


    sys.exit()

    for x in info_values:
        print(x[11])

    # add key=value pair INFO fields to df as separate columns
    for idx, name in enumerate(info_keys):
        print(idx, name)
        vcf_df[name] = [value[idx] for value in info_values]

    print(len(info_keys))

    # print(info_values)

    print(vcf_df)

    print(vcf_df.columns)

    # for x in csq_values:
    #     print(x)

    sys.exit()


    return header, vcf_df







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
        '-v', '--vcf', help='annotated vcf file'
    )
    parser.add_argument(
        '-a', '--analysis', required=False,
        help='name of analysis to display in summary'
    )
    parser.add_argument(
        '-w', '--workflow', required=False,
        help='id of workflow to display in summary'
    )


    return parser.parse_args()


def main():
    args = parse_args()

    header, vcf_df = read_vcf(args.vcf)

    print(vcf_df)


if __name__ == "__main__":
    main()
