import re
import subprocess
import sys
from typing import Union

import numpy as np
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

    @staticmethod
    def switch_include_exclude(filter_str, x, y) -> str:
        """
        Given a bcftools filter string, switches -e and -i to do inverse
        of specified filter

        Parameters
        ----------
        filter_str : str
            string of filters to switch

        Returns
        ----------
        filter_str : str
            string of filters
        """
        return y.join(part.replace(y, x) for part in filter_str.split(x))


    def filter(self, vcf) -> Union[pd.DataFrame, pd.DataFrame]:
        """
        Filter given vcf using bcftools

        Parameters
        ----------
        vcf : pathlib.PosixPath
            path to vcf file to filter

        Returns
        -------
        vcf_df : pd.DataFrame
            dataframe of variants to retain
        filtered_vcf_rows : pd.DataFrame
            dataframe of variants filtered out
        """
        filter_keep = self.args.filter
        filter_out = self.switch_include_exclude(self.args.filter, '-i', '-e')

        # write to temporary vcf files to read from with vcf.read()
        filter_keep += f" {vcf} > tmp1.vcf"
        filter_out += f" {vcf} > tmp2.vcf"

        keep_variants = subprocess.run(
            filter_keep, shell=True,
            capture_output=True
        )

        out_variants = subprocess.run(
            filter_out, shell=True,
            capture_output=True, encoding='UTF-8'
        )

        assert keep_variants.returncode == 0, (
            f"\n\tError in filtering VCF with bcftools\n"
            f"\n\tVCF: {vcf}\n"
            f"\n\tExitcode:{keep_variants.returncode}\n"
            f"\n\tbcftools filter command used: {self.args.filter}\n"
            f"\n\t{keep_variants.stderr.decode()}"
        )
        assert out_variants.returncode == 0, (
            f"\n\tError in filtering VCF with bcftools\n"
            f"\n\tVCF: {vcf}.n"
            f"\n\tExitcode:{out_variants.returncode}\n"
            f"\n\tbcftools filter command used: {out_variants}\n"
            f"\n\t{out_variants.stderr.decode()}"
        )


    def retain_variants(self, vcf_df, filter_idxs) -> list:
        """
        Given a vcf and list of indices, check if any to filter are in the
        given list of variants to never filter and remove if so.

        Regions to not apply filters to are taken from tsv file passed with
        self.args.always_keep, and coordinates in file are inclusive.

        Parameters
        ----------
        vcf_df : pd.DataFrame
            dataframe of variants from VCF
        filter_idxs : list
            list of df indices selected to filter out

        Returns:
        filter_idxs : list
            list of df indices selected to filter out
        """
        retain_idxs = []

        for _, row in self.args.always_keep.iterrows():
            # get indices of df of any variants to always retain, then
            # drop these from the list of filter indices
            retain_idxs.extend(
                np.where((
                        vcf_df['CHROM'] == row['chrom']
                    ) & (
                        vcf_df['POS'] >= row['start']
                    ) & (
                        vcf_df['POS'] <= row['end']
                ))[0]
            )

        filter_idxs = list(set(filter_idxs) - set(retain_idxs))

        return filter_idxs

