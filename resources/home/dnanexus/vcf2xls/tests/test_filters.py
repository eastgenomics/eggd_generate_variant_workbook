import argparse
import os
from pathlib import Path
import sys

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



class TestFilters():
    """
    Tests for building and applying filters from filters.py to dataframe(s)
    of variants
    """
    def read_vcf(self):
        """
        Read in vcf with valid argparse NameSpace for testing, call in every
        test function to reset up unmodified vcf to test each method on

        Returns
        -------
        vcf_handler : utils.vcf.vcf
            class instance of vcf from utils
        """
        # test data vcf
        columns_vcf = os.path.join(TEST_DATA_DIR, "column_methods_test.vcf")

        # initialise vcf class with a valid argparse input to allow
        # calling .read()
        vcf_handler = vcf(argparse.Namespace(
            add_name=True, analysis='',
            filter=None, keep=False, merge=False,
            reorder=[], exclude=None, include=None,
            out_dir='', output='', always_keep=pd.DataFrame(),
            panel='', print_columns=False, reads='', rename=None,
            sample='', sheets=['variants'], summary=None, usable_reads='',
            vcfs=[columns_vcf], workflow=('', '')
        ))
        vcf_df, csq_fields = vcf_handler.read(columns_vcf)

        # call methods like in vcf.read() to process df and add to self
        vcf_df, expanded_vcf_rows = splitColumns.csq(vcf_df, csq_fields)
        vcf_df = splitColumns.info(vcf_df)
        vcf_df = splitColumns.format_fields(vcf_df)

        # set correct dtypes, required for setting numeric & object types
        # to ensure correct filtering filtering
        vcf_df = vcf_handler.set_types(vcf_df)
        vcf_handler.vcfs.append(vcf_df)

        return vcf_handler


    def test_building_filters(self):
        """
        Test filters passed to args are correctly built
        """
        vcf_handler = self.read_vcf()

        # set some example filters as passed to args
        vcf_handler.args.filter = [
            "CHROM==5", "gnomAD_AF>0.02", "Consequence!=missense_variant",
            "AN<2", "ReadPosRankSum<=0.1", "FS>=4"
        ]

        # the above should be turned into a list of lists, with column,
        # operand and value as separate elements
        correct_filter_format = [
            ["CHROM", "==", "5"], ["gnomAD_AF", ">", "0.02"],
            ["Consequence", "!=", "missense_variant"], ["AN", "<", "2"],
            ["ReadPosRankSum", "<=", "0.1"], ["FS", ">=", "4"]
        ]

        filter_handler = filter(vcf_handler.args, vcf_handler.vcfs)
        filter_handler.build()

        assert filter_handler.filters == correct_filter_format, (
            "Filters incorrectly parsed with filters.build()"
        )


    def test_keep_argument_adds_extra_df_filtered_rows(self):
        """
        Tests that when --keep is passed, the filtered rows are saved to the
        last dataframe in self.vcfs and sheets list has "filtered" appended
        """
        vcf_handler = self.read_vcf()
        vcf_handler.args.filter = ["CHROM==chr17"]
        vcf_handler.args.keep = True

        filter_handler = filter(vcf_handler.args, vcf_handler.vcfs)
        filter_handler.build()
        filter_handler.filter()


        # check sheets name added to list of names
        assert filter_handler.args.sheets[-1] == "filtered", (
            "sheets list does not have 'filtered' added with --keep passed"
        )

        assert len(filter_handler.vcfs) == 2, (
            'filtered variants df not retained when args.keep specified'
        )


    def test_correct_rows_filtered_with_eq(self):
        """
        Test when using equal operand filter is correctly applied
        """
        vcf_handler = self.read_vcf()
        vcf_handler.args.filter = ["CHROM==chr17"]
        vcf_handler.args.keep = True

        filter_handler = filter(vcf_handler.args, vcf_handler.vcfs)
        filter_handler.build()
        filter_handler.filter()

        assert filter_handler.vcfs[0]['CHROM'].tolist() == ['chr17'], (
            "Filtering with '==' operand returned incorrect rows"
        )

        assert all([x != ['chr17'] for x in filter_handler.vcfs[-1]['CHROM'].tolist()]), (
            "Filtering with '==' operand returned incorrect rows"
        )


    def test_correct_rows_filtered_with_not_eq(self):
        """
        Test when using not equal operand filter is correctly applied
        """
        vcf_handler = self.read_vcf()
        vcf_handler.args.filter = ["CHROM!=chr17"]
        vcf_handler.args.keep = True

        filter_handler = filter(vcf_handler.args, vcf_handler.vcfs)
        filter_handler.build()
        filter_handler.filter()

        assert filter_handler.vcfs[0]['CHROM'].tolist() != ['chr17'], (
            "Filtering with '!=' operand returned incorrect rows"
        )

        assert all([x != ['chr17'] for x in filter_handler.vcfs[-1]['CHROM'].tolist()]), (
            "Filtering with '!=' operand returned incorrect rows"
        )


    def test_correct_rows_filtered_with_gt(self):
        """
        Tests the correct rows are filtered with one filter applied, and no
        rows are dropped
        """
        vcf_handler = self.read_vcf()
        vcf_handler.args.filter = ["DP>1500"]
        vcf_handler.args.keep = True

        total_rows_before_filter = len(vcf_handler.vcfs[0].index)

        filter_handler = filter(vcf_handler.args, vcf_handler.vcfs)
        filter_handler.build()
        filter_handler.filter()

        total_rows_after_filter = sum([len(x.index) for x in filter_handler.vcfs])

        # check no rows dropped
        assert total_rows_after_filter == total_rows_before_filter, (
            "Mismatch between total rows before and after filtering"
        )

        assert all([x > 1500 for x in filter_handler.vcfs[0]['DP'].tolist()]), (
            'Rows not filtered out correctly when filtering with DP > 1500'
        )

        assert all([x <= 1500 for x in filter_handler.vcfs[-1]['DP'].tolist()]), (
            "Rows incorrectly filtered out when filtering with DP > 1500"
        )


    def test_correct_rows_filtered_w_two_filters(self):
        """
        Tests the correct rows are filtered with two filter applied, rows are
        not double counted and filtered twice and no rows are dropped
        """
        print('testing 2 filters')
        vcf_handler = self.read_vcf()
        vcf_handler.args.filter = ["DP>1500", "MPOS<35"]
        vcf_handler.args.keep = True

        print(vcf_handler.vcfs[0])
        total_rows_before_filter = len(vcf_handler.vcfs[0].index)

        filter_handler = filter(vcf_handler.args, vcf_handler.vcfs)
        filter_handler.build()
        filter_handler.filter()

        print(filter_handler.vcfs[0])
        print(filter_handler.vcfs[-1])

        print(filter_handler.vcfs[0].dtypes)

        total_rows_after_filter = sum([len(x.index) for x in filter_handler.vcfs])






if __name__ == "__main__":
    t = TestFilters()
    t.read_vcf()
    t.test_building_filters()
    t.test_correct_rows_filtered_with_gt()
    t.test_keep_argument_adds_extra_df_filtered_rows()
    t.test_correct_rows_filtered_with_eq()
    t.test_correct_rows_filtered_w_two_filters()

