import argparse
from multiprocessing.dummy import Namespace
import os
from pathlib import Path
import sys
from types import new_class
from bleach import clean

import pytest

sys.path.append(os.path.abspath(
    os.path.join(os.path.realpath(__file__), '../../')
))

from utils.vcf import vcf
from utils.columns import splitColumns
from tests import TEST_DATA_DIR

# initialise vcf class that contains functions for parsing header
vcf_handler = vcf(argparse.Namespace)


class TestHeader():
    """
    Tests for checking header attributes are correct from parse_header()
    """
    header_test_vcf = os.path.join(TEST_DATA_DIR, "header_test.vcf")

    # read in header from our test vcf, call functions to parse reference and
    # from read in header
    header, columns = vcf_handler.parse_header(header_test_vcf)
    vcf_handler.parse_reference(header)


    def test_column_names(self):
        """
        Tests the column names read in from VCF are correctly parsed to a list
        """
        print("Checking column names")
        correct_columns = [
            'CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL',
            'FILTER', 'INFO', 'FORMAT', 'SAMPLE'
        ]

        assert self.columns == correct_columns, (
            f"column names wrongly parsed from {self.header_test_vcf}"
        )


    def test_parse_reference(self):
        """
        Tests the reference is correctly parsed from header lines and stored in
        refs attribute.

        Found as line in file header as: ##reference=file://genome/hs37d5.fa
        """
        assert vcf_handler.refs == ['hs37d5.fa']


    def test_only_header_parsed(self):
        """
        Tests that only header lines parsed from file
        """
        assert all([x.startswith('#') for x in self.header])


class TestDataFrameActions():
    """
    Tests for functions that modify dataframes of variants (i.e. reorder(),
    rename(), merge()...)
    """
    # test data vcf
    columns_vcf = os.path.join(TEST_DATA_DIR, "column_methods_test.vcf.gz")

    # names for intermediary vcfs
    split_vcf = f"{Path(columns_vcf).stem}.split.vcf"
    split_vcf_gz = f"{Path(columns_vcf).stem}.split.vcf.gz"

    def read_vcf(self):
        """
        Read in vcf with valid argparse NameSpace for testing, call in every
        test function to reset up unmodified vcf to test each method on

        Returns
        -------
        vcf_handler : utils.vcf.vcf
            class instance of vcf from utils
        """
        # initialise vcf class with a valid argparse input to
        # allow calling .read()
        vcf_handler = vcf(argparse.Namespace(
            add_name=True, analysis='',
            filter=None, keep=False, merge=False,
            reorder=[], exclude=None, include=None,
            out_dir='', output='',
            panel='', print_columns=False, print_header=False, reads='',
            rename=None, sample='', sheets=['variants'], summary=None,
            vcfs=[self.columns_vcf], workflow=('', '')
        ))

        # first split multiple transcript annotation to separate VCF
        # records, and separate CSQ fields to separate INFO fields
        vcf_handler.bcftools_pre_process(self.columns_vcf, self.split_vcf)
        vcf_handler.bgzip(self.split_vcf)

        vcf_df = vcf_handler.read(self.columns_vcf, Path(self.columns_vcf).stem)

        # call methods like in vcf.read() to process df and add to self
        vcf_df = splitColumns().split(vcf_df)
        vcf_handler.vcfs.append(vcf_df)

        vcf_handler.add_hyperlinks()

        return vcf_handler


    def clean_up(self):
        """
        Reading and processing vcf with bcftools for every test creates
        intermediate files, delete these to not end up with a mess
        """
        pass
        os.remove(self.split_vcf)
        os.remove(self.split_vcf_gz)



    def test_drop_columns_exclude(self):
        """
        Test calling drop_columns() with --exclude removes the given columns
        """
        vcf_handler = self.read_vcf()

        # update namespace with columns to test dropping
        vcf_handler.args.exclude = ["SYMBOL", "ClinVar_CLNSIG", "ID"]

        vcf_handler.drop_columns()

        # ensure columns in exclude are no longer in df columns
        columns = sorted(list(
            set(["Allele", "Gene", "HGNC"]) -
            set(vcf_handler.vcfs[0].columns.tolist())
        ))

        self.clean_up()

        assert columns == ["Allele", "Gene", "HGNC"], (
            "Columns specified with --exlcude not dropped from dataframe"
        )


    def test_drop_columns_include(self):
        """
        Test calling drop_columns() with --include retains only the
        specified columns
        """
        vcf_handler = self.read_vcf()

        # update namespace with columns to test dropping
        include_cols = ['CHROM', 'POS', 'REF', 'ALT']
        vcf_handler.args.include = include_cols

        vcf_handler.drop_columns()

        self.clean_up()

        assert vcf_handler.vcfs[0].columns.tolist() == include_cols, (
            "Columns after calling vcf.drop_columns() with --include not "
            "as expected"
        )


    def test_reorder_columns_correct_order(self):
        """
        Tests that when vcf.order_columns() called the columns are in the
        expected order
        """
        vcf_handler = self.read_vcf()

        # capture original order to compare against
        original_cols = vcf_handler.vcfs[0].columns.tolist()

        order_cols = ["CSQ_SYMBOL", "CSQ_ClinVar_CLNSIG", "ID", "CHROM", "POS"]
        vcf_handler.args.reorder = order_cols

        vcf_handler.order_columns()

        self.clean_up()

        assert vcf_handler.vcfs[0].columns.tolist()[:5] == order_cols, (
            "Reordered columns are not in correct order"
        )


    def test_reorder_columns_no_dropped_columns(self):
        """
        Tests when reorder columns, no columns are dropped
        """
        vcf_handler = self.read_vcf()

        # capture original order to compare against
        original_cols = vcf_handler.vcfs[0].columns.tolist()

        order_cols = ["CSQ_SYMBOL", "CSQ_ClinVar_CLNSIG", "ID", "CHROM", "POS"]
        vcf_handler.args.reorder = order_cols

        vcf_handler.order_columns()

        # sort both original and new columns to compare
        original_cols = sorted(original_cols)
        new_cols = sorted(vcf_handler.vcfs[0].columns.tolist())

        self.clean_up()

        assert original_cols == new_cols, (
            "Columns dropped wrongly when reordering"
        )

    @pytest.fixture
    def rename(self):
        """
        Call rename_columns() and store the previous columns, new columns and
        renaming dict in object Namespace to for other rename() tests
        """
        vcf_handler = self.read_vcf()

        prev_columns = vcf_handler.vcfs[0].columns.tolist()

        rename_dict = {
            "CHROM": "chrom",
            "POS": "pos",
        }

        vcf_handler.args.rename = rename_dict
        vcf_handler.rename_columns()

        # get column list after renaming, replacing ' ' in each name with '_'
        # to undo what is done in .rename_columns() for comparing against
        # original names
        new_columns = vcf_handler.vcfs[0].columns.tolist()
        new_columns = [x.replace(' ', '_') for x in new_columns]

        # add rename variable to Namespace to access in other tests of renaming
        rename_vals = Namespace()
        rename_vals.prev_columns = prev_columns
        rename_vals.new_columns = new_columns
        rename_vals.rename_dict = rename_dict

        self.clean_up()

        return rename_vals


    def test_non_rename_columns_unaffacted(self, rename):
        """
        Tests columns not specified with --rename are unaffacted
        """

        prev_non_rename_cols = sorted(
            list(set(rename.prev_columns) - set(rename.rename_dict.keys())))
        new_non_rename_cols = sorted(
            list(set(rename.new_columns) - set(rename.rename_dict.values())))

        assert new_non_rename_cols == prev_non_rename_cols, (
            "Renaming columns has affected more than specified columns"
        )


    def test_renamed_correctly(self, rename):
        """
        Tests that the renamed columns have the specified names
        from args.rename
        """
        renamed_columns = sorted(
            list(set(rename.new_columns) & set(rename.rename_dict.values())))

        assert renamed_columns == sorted(list(rename.rename_dict.values())), (
            "Column names specified in args.rename not present in renamed "
            "column list"
        )


if __name__ == "__main__":
    header = TestHeader()
    header.test_column_names()

    df_actions = TestDataFrameActions()
    df_actions.read_vcf()
    # df_actions.test_exclude()
