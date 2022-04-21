import argparse
import os
from pathlib import Path
import shutil
import sys
from typing import Union

import pandas as pd
import pytest

sys.path.append(os.path.abspath(
    os.path.join(os.path.realpath(__file__), '../../')
))

from utils.columns import splitColumns
from utils.filters import filter
from utils.vcf import vcf
from tests import TEST_DATA_DIR

# vcf we are using for testing, ~5000 variants with multiple transcript
# annotation for each
TEST_VCF = "NA12878_unittest.vcf"


class TestModifyingFieldTypes():
    """
    Tests for modifying the INFO field types in the header specified in the
    --types arg.

    This tests the functions from filters.py:
        - modify_header_types()
        - write_header()
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
        filter(self.vcf_handler.args).write_header(test_vcf, self.new_header)

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
    test_vcf = os.path.join(TEST_DATA_DIR, TEST_VCF)

    # names for intermediary vcfs
    split_vcf = f"{Path(test_vcf).stem}.split.vcf"
    split_vcf_gz = f"{Path(test_vcf).stem}.split.vcf.gz"
    filter_vcf = f"{Path(test_vcf).stem}.filter.vcf"
    filter_vcf_gz = f"{Path(test_vcf).stem}.filter.vcf.gz"

    # initialise vcf class with a valid argparse input to allow
    # calling .read()
    vcf_handler = vcf(argparse.Namespace(
        add_name=False, analysis='',
        filter=None, keep=False, merge=False,
        reorder=[], exclude=None, include=None,
        out_dir='', output='', always_keep=pd.DataFrame(),
        panel='', print_columns=False, print_header=False, reads='',
        rename=None, sample='', sheets=['variants'], summary=None,
        usable_reads='', vcfs=[test_vcf], workflow=('', ''), types=None
    ))


    def filter(self, filter_str) -> Union[pd.DataFrame, pd.DataFrame]:
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
        # process with bcftools +split-vep ready for filtering
        _, columns = self.vcf_handler.parse_header(self.split_vcf)
        self.vcf_handler.bcftools_pre_process(self.test_vcf, self.split_vcf)
        self.vcf_handler.bgzip(self.split_vcf)

        self.vcf_handler.args.filter = filter_str
        self.vcf_handler.args.keep = True

        filter_handle = filter(self.vcf_handler.args)

        # filter vcf against specified filters using bcftools
        filter_handle.filter(self.split_vcf_gz, self.filter_vcf, columns)
        self.vcf_handler.bgzip(self.filter_vcf)

        # filters.filter() writes temp filtered vcf containing the
        # filtered variants to read into df
        variant_df = self.vcf_handler.read(self.filter_vcf, Path(self.test_vcf).stem)

        # get filtered out rows and read back to new dfs
        keep_df, filtered_df = filter_handle.split_include_exclude(variant_df)

        # split out INFO and FORMAT column values to individual
        # columns in dataframe
        keep_df = splitColumns().split(keep_df)
        filtered_df = splitColumns().split(filtered_df)

        # delete the filtered vcf file
        os.remove(self.split_vcf)
        os.remove(self.split_vcf_gz)
        os.remove(self.filter_vcf)
        os.remove(self.filter_vcf_gz)

        return keep_df, filtered_df


    def test_filter_with_include_eq(self):
        """
        Test when using include with equal operator filter is correctly applied
        """
        keep_df, filtered_df = self.filter(
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


    def test_filter_with_exclude_eq(self):
        """
        Test when using exclude with equal operator filter is correctly applied
        """
        keep_df, filtered_df = self.filter(
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


    def test_filter_with_exclude_gt(self):
        """
        Test when using exclude with gt operator filter is correctly applied
        """
        keep_df, filtered_df = self.filter(
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


    def test_combined_exclude_float_and_string(self):
        """
        Test filtering on gnomAD at 2%, and filtering out synonymous and
        intronic variants
        """
        keep_df, filtered_df = self.filter(
            "bcftools filter -e 'CSQ_gnomAD_AF>0.01 "
            "| CSQ_gnomADg_AF>0.01 "
            "| CSQ_Consequence=\"synonymous_variant\" "
            "| CSQ_Consequence=\"intron_variant\"'"
        )

        # set '.' to 0 to allow column to be a float for comparing
        keep_df['CSQ_gnomAD_AF'] = keep_df['CSQ_gnomAD_AF'].apply(
            lambda x: '0' if x == '.' else x
        ).astype(float)
        keep_df['CSQ_gnomADg_AF'] = keep_df['CSQ_gnomADg_AF'].apply(
            lambda x: '0' if x == '.' else x
        ).astype(float)
        filtered_df['CSQ_gnomAD_AF'] = filtered_df['CSQ_gnomAD_AF'].apply(
            lambda x: '0' if x == '.' else x
        ).astype(float)
        filtered_df['CSQ_gnomADg_AF'] = filtered_df['CSQ_gnomADg_AF'].apply(
            lambda x: '0' if x == '.' else x
        ).astype(float)

        # check we have correctly filtered variants
        assert all(keep_df['CSQ_gnomAD_AF'] <= 0.01) & \
                all(keep_df['CSQ_Consequence'] != 'synonymous_variant') & \
                    all(keep_df['CSQ_Consequence'] != 'intron_variant'), (
            "Filtering to exclude gnomAD_AF>0.01 and synonymous/intronic "
            "variants did not filter out the correct variants"
        )


        filtered_df['check'] = filtered_df.apply(
                lambda x: x['CSQ_gnomAD_AF'] > 0.01 or \
                x['CSQ_gnomADg_AF'] > 0.01 or \
                x['CSQ_Consequence'] == 'synonymous_variant' or \
                x['CSQ_Consequence'] == 'intron_variant', axis=1
            )

        df = filtered_df[[
            'CHROM', 'POS', 'REF', 'ALT', 'CSQ_SYMBOL', 'CSQ_Feature',
            'CSQ_gnomAD_AF', 'CSQ_gnomADg_AF', 'CSQ_Consequence', 'check']]

        print(df[df['check'] == False])

        print(len(keep_df))
        print(len(filtered_df))
        print(len(keep_df) + len(filtered_df))


        assert all(
            filtered_df.apply(
                lambda x: x['CSQ_gnomAD_AF'] > 0.01 or \
                x['CSQ_gnomADg_AF'] > 0.01 or \
                x['CSQ_Consequence'] == 'synonymous_variant' or \
                x['CSQ_Consequence'] == 'intron_variant', axis=1
            )
        ), (
            "Filtering to exclude gnomAD_AF>0.01 and synonymous/intronic "
            "variants filtered out the wrong variants"
        )



if __name__ == "__main__":

    # modify_header = TestModifyingFieldTypes()
    # modify_header.test_type_correctly_modified()
    # modify_header.test_header_overwritten_correctly()


    t = TestFilters()

    # t.test_filter_with_include_eq()
    # t.test_filter_with_exclude_eq()
    # t.test_filter_with_exclude_gt()
    t.test_combined_exclude_float_and_string()
