import argparse
from multiprocessing.dummy import Namespace
import os
from pathlib import Path
import sys
import pandas as pd
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
    header_test_vcf = os.path.join(TEST_DATA_DIR, "NA12878_unittest.vcf")

    # read in header from our test vcf, call functions to parse reference and
    # from read in header
    header, columns = vcf_handler.parse_header(header_test_vcf)


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


    def test_parse_reference_vep(self):
        """
        Tests the reference is correctly parsed from header lines and stored in
        refs attribute when ##VEP is present in header, reference is stored
        as 'assembly' string in VEP header line
        """
        vcf_handler.refs = []

        vcf_handler.parse_reference(self.header)
        assert vcf_handler.refs == ['GRCh37.p13'], (
            'Reference build not correctly parsed from VEP string in header'
        )


    def test_parse_reference_no_vep(self):
        """
        Tests the reference is correctly parsed from header lines and stored in
        refs attribute when ##VEP is not present in header

        Found as line in file header as: ##reference=file://genome/hs37d5.fa
        """
        vcf_handler.refs = []
        vep_header = [x for x in self.header if not x.startswith('##VEP')]

        vcf_handler.parse_reference(vep_header)
        assert vcf_handler.refs == ['hs37d5.fa'], (
            'Reference build not correctly parsed from header'
        )


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
        vcf_handler.args.exclude = ["POS", "AF", "ID"]

        vcf_handler.drop_columns()

        # ensure columns in exclude are no longer in df columns
        columns = sorted(list(
            set(["AF", "ID", "POS"]) -
            set(vcf_handler.vcfs[0].columns.tolist())
        ))

        self.clean_up()

        assert columns == ["AF", "ID", "POS"], (
            "Columns specified with --exclude not dropped from dataframe"
        )


    def test_drop_columns_include(self):
        """
        Test calling drop_columns() with --include retains only the
        specified columns
        """
        vcf_handler = self.read_vcf()

        # update namespace with columns to test dropping everything else
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

        order_cols = ["ALT", "REF", "ID", "CHROM", "POS"]
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

        order_cols = ["ALT", "REF", "ID", "CHROM", "POS"]
        vcf_handler.args.reorder = order_cols

        # capture original order to compare against with ordered cols removed
        original_cols = vcf_handler.vcfs[0].columns.tolist()
        original_cols = [x for x in original_cols if x not in order_cols]

        vcf_handler.order_columns()

        self.clean_up()

        assert vcf_handler.vcfs[0].columns.tolist()[5:] == original_cols, (
            "Original columns after reordering not in correct order"
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


class TestHyperlinks():
    '''
    Tests to check hyperlinks are generated correctly
    '''
    @staticmethod
    def test_column_creation():
        '''
        Test that the function to create a "DECIPHER" column works
        '''
        # Create test dataframe with no "DECIPHER" columnn
        df = pd.DataFrame([
            {'CHROM': 1, 'POS': 64883298, 'REF': 'T', 'ALT': 'C'},
            {'CHROM': 10, 'POS': 27035066, 'REF': 'C', 'ALT': 'T'},
        ])

        to_add_column = vcf(argparse.Namespace())
        to_add_column.vcfs = [df]
        to_add_column.refs = ['38']  # Set reference = build 38

        # Call add_decipher_column which should add an empty column titled
        # 'DECIPHER' to this dataframe
        vcf.add_decipher_column(to_add_column)

        # Assert statement to check that DECIPHER column exists
        assert "DECIPHER" in to_add_column.vcfs[0].columns, (
            'DECIPHER column not created by make_decipher_columns() function'
        )

    @staticmethod
    def test_decipher_column_added():
        """
        Test to ensure that the DECIPHER column is created when --decipher
        input is specified.
        """
        # Create test dataframe with no "DECIPHER" columnn
        df = pd.DataFrame([
            {'CHROM': 1, 'POS': 64883298, 'REF': 'T', 'ALT': 'C'},
            {'CHROM': 10, 'POS': 27035066, 'REF': 'C', 'ALT': 'T'},
        ])

        # Include --decipher input, all other inputs are false to avoid giving
        # this dataframe to other functions except make_decipher_columns
        should_have_decipher_column = vcf(argparse.Namespace(
            additional_files=False, filter=False, print_columns=False,
            rename=False, vcfs=[], merge=False, include=False, exclude=False,
            reorder=False, decipher=True,  # Set DECIPHER = True
            ))
        should_have_decipher_column.vcfs = [df]
        should_have_decipher_column.refs = ['38']  # Set build = 38

        # Call process() which should call make_decipher_columns()
        vcf.process(should_have_decipher_column)

        # Assert statement to check that DECIPHER column exists
        assert "DECIPHER" in should_have_decipher_column.vcfs[0].columns, (
            'DECIPHER column not created despite --decipher input given'
        )

    @staticmethod
    def test_decipher_column_not_added():
        """
        Test to ensure that the DECIPHER column is not made if --decipher is
        not specified.
        """
        df = pd.DataFrame([
            {'CHROM': 1, 'POS': 64883298, 'REF': 'T', 'ALT': 'C'},
            {'CHROM': 10, 'POS': 27035066, 'REF': 'C', 'ALT': 'T'},
        ])
        # Using the same test dataframe but without --decipher input, all other
        # inputs are false to avoid giving this dataframe to other functions
        # except make_decipher_columns
        should_not_have_decipher_column = vcf(argparse.Namespace(
            additional_files=False, filter=False, print_columns=False,
            rename=False, vcfs=[], merge=False, include=False, exclude=False,
            reorder=False, decipher=False,  # Set DECIPHER = False
            ))
        should_not_have_decipher_column.vcfs = [df]
        should_not_have_decipher_column.refs = ['38']  # Set build = 38

        # Call process() which should not call make_decipher_columns()
        vcf.process(should_not_have_decipher_column)

        # Assert statement to check that DECIPHER column was not added
        assert "DECIPHER" not in should_not_have_decipher_column.vcfs[0].columns, (
            'DECIPHER column created despite --decipher input not given'
        )

    @staticmethod
    def test_decipher_build_37():
        '''
        Test that no DECIPHER column is created if the build is 37
        '''
        # Intialise test dataframe with build 37 genome positions
        df = pd.DataFrame([
            {'CHROM': 1, 'POS': 1273116, 'REF': 'A',
             'ALT': 'G'},
            {'CHROM': 12, 'POS': 103241070, 'REF': 'T',
             'ALT': 'C'},
        ])

        build_37_vcf = vcf(argparse.Namespace(
            additional_files=False, filter=False, print_columns=False,
            rename=False, vcfs=[], merge=False, include=False, exclude=False,
            reorder=False, decipher=True  # Set DECIPHER = True
        ))

        build_37_vcf.vcfs = [df]
        build_37_vcf.refs = ['37']  # Set reference = build 37

        # Call process() which should call make_decipher_column(), however no
        # column should be created as make_decipher_column() has if statement
        # to skip build 37
        vcf.process(build_37_vcf)

        # Assert statement to check that DECIPHER column was not added
        assert "DECIPHER" not in build_37_vcf.vcfs[0].columns, (
            'DECIPHER column created despite variants being in build 37. '
            'DECIPHER does not store build 37 data.'
        )

    @staticmethod
    def test_decipher_links_build_38():
        '''
        Test that the DECIPHER links are generated correctly
        '''
        # Intialise test dataframe with build 38 genome positions
        # DECIPHER column is empty as that is the input to generate hyperlinks
        df = pd.DataFrame([
            {'CHROM': 1, 'POS': 64883298, 'REF': 'T',
             'ALT': 'C', 'DECIPHER': ''},
            {'CHROM': 10, 'POS': 27035066, 'REF': 'C',
             'ALT': 'T', 'DECIPHER': ''},
        ])

        test_vcf = vcf(argparse.Namespace(decipher=True))
        test_vcf.vcf = [df]
        test_vcf.refs = ['38']  # Set reference = build 38

        # Call function to add hyperlinks
        vcf.add_hyperlinks(test_vcf)

        # Define expected string output
        valid_string = (
            '=HYPERLINK("https://www.deciphergenomics.org/sequence-variant/1-6'
            '4883298-T-C", "1-64883298-T-C")'
        )

        # Assert that output is as expected
        assert test_vcf.vcfs[0]["DECIPHER"][0] == valid_string, (
            "DECIPHER link output incorrect"
        )

    @staticmethod
    def test_gnomad_build_37():
        '''
        Test that the gnomAD links are generated correctly for build 37
        '''
        # Intialise test dataframe with build 37 genome positions
        df = pd.DataFrame([
            {'CHROM': 1, 'POS': 1271940, 'REF': 'C',
             'ALT': 'T', 'gnomAD': 0.0108147},
        ])

        test_vcf = vcf(argparse.Namespace())
        test_vcf.vcfs = [df]
        test_vcf.refs = ['37']  # Set reference = build 37

        # Call function to add hyperlinks
        vcf.add_hyperlinks(test_vcf)

        valid_string = (
            '=HYPERLINK("https://gnomad.broadinstitute.org/variant/1-1271940-C'
            '-T?dataset=gnomad_r2_1", 0.0108147)'
            )

        # Assert the output is as expected
        assert test_vcf.vcfs[0]["gnomAD"][0] == valid_string, (
            "gnomAD AF link output incorrect for build 37 input"
        )

    @staticmethod
    def test_gnomad_build_38():
        '''
        Test that the gnomAD links are generated correctly for build 38
        '''
        # Intialise test dataframe with build 38 genome positions
        df = pd.DataFrame([
            {'CHROM': 1, 'POS': 64883298, 'REF': 'T',
             'ALT': 'C', 'gnomADg AF': 0.0004271},
        ])

        test_vcf = vcf(argparse.Namespace())
        test_vcf.vcfs = [df]
        test_vcf.refs = ['38']  # Set reference = build 38

        # Call function to add hyperlinks
        vcf.add_hyperlinks(test_vcf)

        valid_string = (
            '=HYPERLINK("https://gnomad.broadinstitute.org/variant/1-64883298'
            '-T-C?dataset=gnomad_r3", 0.0004271)'
        )

        # Assert the output is as expected
        assert test_vcf.vcfs[0]["gnomADg AF"][0] == valid_string, (
            "gnomAD AF link output incorrect for build 38 input"
        )

if __name__ == "__main__":
    header = TestHeader()
    header.test_column_names()

    df_actions = TestDataFrameActions()

    hyperlinks = TestHyperlinks()
