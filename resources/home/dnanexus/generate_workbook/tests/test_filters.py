import argparse
import os
from pathlib import Path
from random import shuffle
import shutil
import sys
from typing import Union

import numpy as np
import pandas as pd
import pytest

sys.path.append(os.path.abspath(
    os.path.join(os.path.realpath(__file__), '../../')
))

from utils.columns import splitColumns
from utils.filters import filter
from utils.vcf import vcf
from tests import TEST_DATA_DIR

# initialise vcf class that contains functions for parsing header
vcf_handler = vcf(argparse.Namespace)

# vcf we are using for testing, ~5000 variants with multiple transcript
# annotation for each
TEST_VCF = "NA12878_unittest.vcf"


class TestModifyingFieldTypes():
    """
    Tests for modifying the INFO field types in the header using --types arg
    """
    # test data vcf
    columns_vcf = os.path.join(TEST_DATA_DIR, TEST_VCF)

    # initialise vcf class with a valid argparse input to allow
    # calling .read()
    vcf_handler = vcf(argparse.Namespace(
        add_name=True, analysis='',
        filter=None, keep=False, merge=False,
        reorder=[], exclude=None, include=None,
        out_dir='', output='', always_keep=pd.DataFrame(),
        panel='', print_columns=False, reads='', rename=None,
        sample='', sheets=['variants'], summary=None, usable_reads='',
        vcfs=[columns_vcf], workflow=('', ''),
        types={'gnomADg_AF': 'Float'}
    ))

    # get header from vcf
    file_header, _ = vcf_handler.parse_header(columns_vcf)
    current_file_header = file_header

    # get new header
    new_header = filter(vcf_handler.args).modify_header_types(columns_vcf)


    def test_type_correctly_modified(self):
        """
        Test that the header is as expected
        """
        header_diff = list(set(self.new_header) - set(self.file_header))

        correct_diff = [(
            '##INFO=<ID=gnomADg_AF,Number=.,Type=Float,Description=\"AF field '
            'from /opt/vep/.vep/gnomad.genomes.r2.0.1.sites.noVEP.vcf.gz\">'
        )]

        assert header_diff == correct_diff, (
            "Header wrongly modified with modify_header_types()"
        )


    def test_header_overwritten_correctly(self):
        """
        Test that writing the new header with modified types back to file
        """
        # make a copy of the vcf to test modifying
        test_vcf = f"{self.columns_vcf.replace('.vcf', '.test_filters.tmp.vcf')}"
        shutil.copy(self.columns_vcf, test_vcf)

        # call method to overwrite with new header
        filter(vcf_handler.args).write_header(test_vcf, self.new_header)

        with open(self.columns_vcf) as fh:
            # read variants from unmodified testing vcf
            current_vcf_variants = [
                x for x in fh.read().splitlines() if not x.startswith('#')
            ]

        with open(test_vcf) as fh:
            # read back header and variants from modified vcf to test
            written_header = [
                x for x in fh.read().splitlines() if x.startswith('#')
            ]
            fh.seek(0)
            new_vcf_variants = [
                x for x in fh.read().splitlines() if not x.startswith('#')
            ]

        assert written_header == self.new_header, (
            "Modified header not correctly written to file"
        )

        assert current_vcf_variants == new_vcf_variants, (
            "Variants modified when modifying types in header"
        )

        os.remove(test_vcf)




class TestFilters():
    """
    Tests for building and applying filters from filters.py to dataframe(s)
    of variants
    """
    # test data vcf
    columns_vcf = os.path.join(TEST_DATA_DIR, TEST_VCF)

    # initialise vcf class with a valid argparse input to allow
    # calling .read()
    vcf_handler = vcf(argparse.Namespace(
        add_name=False, analysis='',
        filter=None, keep=False, merge=False,
        reorder=[], exclude=None, include=None,
        out_dir='', output='', always_keep=pd.DataFrame(),
        panel='', print_columns=False, print_header=False, reads='',
        rename=None, sample='', sheets=['variants'], summary=None,
        usable_reads='', vcfs=[columns_vcf], workflow=('', ''), types=None
    ))

    # names for output vcfs
    test_processed_vcf = columns_vcf.replace('.vcf', '.split.vcf')
    test_filter_vcf = columns_vcf.replace('.vcf', '.filter.vcf')

    # process with bcftools +split-vep ready for filtering
    vcf_handler.bcftools_pre_process(columns_vcf, test_processed_vcf)


    def filter_vcf_and_read(self, filter_str) -> Union[pd.DataFrame, pd.DataFrame]:
        """
        Given a bcftools filter string, run the filter and get the filtered and
        filtered out rows of the vcf as 2 dataframes

        Parameters
        ----------
        filter_str : str
            string of bcftools filter command

        Returns
        -------
        keep_df : pd.DataFrame
            df of filtered variants
        filtered_df : pd.DataFrame
            df of filtered out variants
        """
        self.vcf_handler.args.filter = filter_str
        self.vcf_handler.args.keep = True

        filter_handle = filter(self.vcf_handler.args)

        # apply filter, read in filtered vcf, then get the filtered out rows
        filter_handle.filter(self.test_processed_vcf, self.test_filter_vcf)

        keep_df = self.vcf_handler.read(
            self.test_filter_vcf, Path(self.columns_vcf).stem
        )
        _, columns = vcf_handler.parse_header(self.test_processed_vcf)

        filtered_df = filter_handle.get_filtered_rows(
            self.test_processed_vcf, keep_df, columns
        )

        # split out INFO and FORMAT column values to individual
        # columns in dataframe
        keep_df = splitColumns().split(keep_df)
        filtered_df = splitColumns().split(filtered_df)

        # delete the filtered vcf file
        os.remove(self.test_filter_vcf)

        return keep_df, filtered_df


    def test_correct_rows_filtered_with_include_eq(self):
        """
        Test when using include with equal operator filter is correctly applied
        """
        keep_df, filtered_df = self.filter_vcf_and_read(
            "bcftools filter -i 'CHROM==\"4\"'"
        )

        assert keep_df['CHROM'].unique().astype(str).tolist() == ["4"], (
            "Filtering with '==' operator returned incorrect rows"
        )

        assert all([
            x != "4" for x in filtered_df['CHROM'].astype(str).tolist()
        ]), (
            "Filtering with bcftools filter -i 'CHROM==\"4\"' operator "
            "filtered incorrect rows"
        )


    def test_correct_rows_filtered_with_exclude_eq(self):
        """
        Test when using exclude with equal operator filter is correctly applied
        """
        keep_df, filtered_df = self.filter_vcf_and_read(
            "bcftools filter -e 'CHROM==\"4\"'"
        )

        assert filtered_df['CHROM'].astype(str).unique().tolist() == ["4"], (
            "Filtering bcftools filter -e 'CHROM==\"4\"' operator returned "
            "incorrect rows"
        )

        assert all([
            x != "4" for x in keep_df['CHROM'].astype(str).tolist()
        ]), (
            "Filtering with bcftools filter -e 'CHROM==\"4\"' operator "
            "filtered out incorrect rows"
        )


    def test_correct_rows_filtered_with_exclude_gt(self):
        """
        Test when using exclude with gt operator filter is correctly applied
        """
        keep_df, filtered_df = self.filter_vcf_and_read(
            "bcftools filter -e 'CSQ_gnomAD_AF>0.02'"
        )

        # set '.' to 0 to allow column to be a float for comparing
        keep_df['CSQ_gnomAD_AF'] = keep_df['CSQ_gnomAD_AF'].apply(
            lambda x: '0' if x == '.' else x
        ).astype(float)

        filtered_df['CSQ_gnomAD_AF'] = filtered_df['CSQ_gnomAD_AF'].astype(float)


        assert all(keep_df['CSQ_gnomAD_AF'] <= 0.02), (
            "Filtering bcftools filter -e 'gnomAD_AF>0.02' operator returned "
            "incorrect rows"
        )

        assert all(filtered_df['CSQ_gnomAD_AF'] > 0.02), (
            "Filtering with bcftools filter -e 'gnomAD_AF>0.02' operator "
            "filtered out incorrect rows"
        )




if __name__ == "__main__":

    modify_header = TestModifyingFieldTypes()
    modify_header.test_type_correctly_modified()
    modify_header.test_header_overwritten_correctly()


    t = TestFilters()
    t.test_correct_rows_filtered_with_include_eq()
    t.test_correct_rows_filtered_with_exclude_eq()
    t.test_correct_rows_filtered_with_exclude_gt()
