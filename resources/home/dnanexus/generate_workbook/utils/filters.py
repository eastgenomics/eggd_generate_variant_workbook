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
    def __init__(self, args=None) -> None:
        self.args = args


    def filter(self, vcf) -> None:
        """
        Filter given vcf using bcftools

        Parameters
        ----------
        vcf : pathlib.PosixPath
            path to vcf file to filter

        Outputs
        -------
        filtered.tmp.vcf : file
            vcf file filtered with bedtools and specified filters

        Raises
        ------
        AssertionError
            Raised when non-zero exitcode returned from bcftools annotate
        """
        # write to temporary vcf files to read from with vcf.read()
        filter = f"{self.args.filter} {vcf} > filtered.tmp.vcf"

        output = subprocess.run(filter, shell=True, capture_output=True)

        assert output.returncode == 0, (
            f"\n\tError in filtering VCF with bcftools\n"
            f"\n\tVCF: {vcf}\n"
            f"\n\tExitcode:{output.returncode}\n"
            f"\n\tbcftools filter command used: {self.args.filter}\n"
            f"\n\t{output.stderr.decode()}"
        )


    def get_filtered_rows(self, vcf, keep_df, columns) -> pd.DataFrame():
        """
        Given filtered vcf dataframe, return a dataframe of the filtered out
        variants

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

        # get the rows only present in original vcf => filtered out rows
        filtered_out_df = full_df.merge(keep_df, how='left', indicator=True)
        filtered_out_df = filtered_out_df.query('_merge == "left_only"')
        filtered_out_df.drop(['_merge'], axis=1, inplace=True)

        return filtered_out_df


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
