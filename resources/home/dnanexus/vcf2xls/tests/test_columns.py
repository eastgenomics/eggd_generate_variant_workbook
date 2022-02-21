import argparse
import os
from pathlib import Path
import subprocess
import sys

import numpy as np
import pandas as pd

import pytest

sys.path.append(os.path.abspath(
    os.path.join(os.path.realpath(__file__), '../../')
))

from utils.vcf import vcf
from utils.columns import splitColumns
from tests import TEST_DATA_DIR


def read_test_vcf(vcf_file):
    """
    Read in test vcf to dataframe using methods from vcf()
    """
    # test data vcf
    columns_vcf = os.path.join(TEST_DATA_DIR, vcf_file)

    # initialise vcf class with a valid argparse input to allow calling .read()
    vcf_handler = vcf(argparse.Namespace(
        add_name=False, analysis='', clinical_indication='', exclude=None,
        filter=None, include=None, keep=False, merge=False,
        out_dir='/home/jethro/Projects/eggd_vcf2xls_nirvana',
        output='/home/jethro/Projects/eggd_vcf2xls_nirvana/NA12878_R29.1_ID.xlsx',
        panel='', print_columns=False, reads='', rename=None, reorder=None,
        sample='', sheets=['variants'], summary=None, usable_reads='',
        vcfs=[columns_vcf], workflow=('', '')
    ))
    vcf_df, csq_fields = vcf_handler.read(columns_vcf)

    return vcf_df, csq_fields


def read_csq_from_vcf() -> list:
    """
    Reads CSQ strings direct from VCF to compare against VCF using
    subprocess with grep, cut and sed

    Returns
    -------
    csq_strings : list
        list of each variants/transcript csq string
    """
    # read in INFO column from VCF to list of records, split the multiple
    # transcripts to separate lines by splitting with sed on commas
    info_strings = subprocess.run((
        f"grep -v '^#' {TEST_DATA_DIR}/multi_transcript_annotation.vcf | "
        "cut -f8 | sed 's/,/\\n/g'")
        , shell=True, capture_output=True).stdout.decode().rstrip('\n').split('\n')

    # retain just CSQ part, results in one item with all CSQ joined by '|'
    csq_strings = [x.split('CSQ=')[-1] for x in info_strings]

    return csq_strings



class TestInfoColumn():
    """
    Tests for splitColumns.info() that splits out key value pairs
    from INFO column to separate columns in dataframe

    INFO values from test vcf:

    DP=1400;ECNT=5;POP_AF=0.001;P_GERMLINE=-0;RPA=11,10;RU=A;STR;TLOD=36.38;
    DP=1400;ECNT=5;POP_AF=0.001;P_GERMLINE=-0;RPA=11,12;RU=A;STR;TLOD=418.64;
    DP=1400;ECNT=5;POP_AF=0.001;P_GERMLINE=-0;RPA=11,13;RU=A;STR;TLOD=21.15;
    DP=1563;ECNT=2;POP_AF=0.001;P_GERMLINE=-0.029;RPA=3,2;RU=CGGCGC;STR;TLOD=3.87;
    DP=2714;ECNT=1;POP_AF=0.001;P_GERMLINE=-0;TLOD=3878.54;
    DP=2965;ECNT=2;POP_AF=0.001;P_GERMLINE=-0;RPA=5,4;RU=GTG;STR;TLOD=87.52;
    DP=3077;ECNT=2;POP_AF=0.001;P_GERMLINE=-0.001;TLOD=5.47;
    DP=2943;ECNT=2;POP_AF=0.001;P_GERMLINE=-0.002;TLOD=5.02;
    DP=2917;ECNT=1;POP_AF=0.001;P_GERMLINE=-0;RPA=4,3;RU=CAC;STR;TLOD=36.68;
    DP=1907;ECNT=1;POP_AF=0.001;P_GERMLINE=-0;TLOD=2904.48;
    DP=1892;ECNT=1;POP_AF=0.001;P_GERMLINE=-0;RPA=8,9;RU=T;STR;TLOD=19.45;
    DP=2358;ECNT=1;POP_AF=0.001;P_GERMLINE=-0;TLOD=74.01;

    Keys from each variant:

    DP ECNT POP_AF P_GERMLINE RPA RU STR TLOD
    DP ECNT POP_AF P_GERMLINE RPA RU STR TLOD
    DP ECNT POP_AF P_GERMLINE RPA RU STR TLOD
    DP ECNT POP_AF P_GERMLINE RPA RU STR TLOD
    DP ECNT POP_AF P_GERMLINE TLOD
    DP ECNT POP_AF P_GERMLINE RPA RU STR TLOD
    DP ECNT POP_AF P_GERMLINE TLOD
    DP ECNT POP_AF P_GERMLINE TLOD
    DP ECNT POP_AF P_GERMLINE RPA RU STR TLOD
    DP ECNT POP_AF P_GERMLINE TLOD
    DP ECNT POP_AF P_GERMLINE RPA RU STR TLOD
    DP ECNT POP_AF P_GERMLINE TLOD

    Unique keys:

    ['STR', 'POP_AF', 'RPA', 'TLOD', 'P_GERMLINE', 'ECNT', 'DP', 'RU']
    """
    # run dataframe through splitColumns.info() to split out INFO column
    vcf_df, _ = read_test_vcf(vcf_file="column_methods_test.vcf")
    vcf_df = splitColumns.info(vcf_df)


    def test_parsed_correct_columns_from_info_records(self) -> None:
        """
        Test that all unique INFO values from all rows correctly parsed and
        set to columns in returned dataframe
        """
        correct_columns = [
            'CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'FORMAT',
            'SAMPLE', 'STR', 'POP_AF', 'RPA', 'TLOD', 'P_GERMLINE',
            'ECNT', 'DP', 'RU'
        ]

        assert sorted(correct_columns) == sorted(self.vcf_df.columns.tolist()), (
            "Number of columns after spliting INFO is incorrect"
        )


    def test_parsed_correct_DP_values(self):
        """
        Test values read into dataframe for DP match the values
        above from the VCF
        """
        correct_dp_values = [
            '1400', '1400', '1400', '1563', '2714', '2965',
            '3077','2943', '2917', '1907', '1892', '2358'
        ]

        assert correct_dp_values == self.vcf_df['DP'].tolist()


    def test_parsed_correct_ECNT_values(self):
        """
        Test values read into dataframe for ECNT match the values
        above from the VCF
        """
        correct_ecnt_values = [
            '5', '5', '5', '2', '1', '2', '2', '2', '1', '1', '1', '1'
        ]

        assert correct_ecnt_values == self.vcf_df['ECNT'].tolist()


    def test_parsed_correct_pop_af_values(self):
        """
        Test values read into dataframe for POP_AF match the values
        above from the VCF
        """
        correct_pop_af_values = [
            '0.001', '0.001', '0.001', '0.001', '0.001', '0.001',
            '0.001', '0.001', '0.001', '0.001', '0.001', '0.001',
        ]

        assert correct_pop_af_values == self.vcf_df['POP_AF'].tolist()


    def test_parsed_correct_p_germline_values(self):
        """
        Test values read into dataframe for P_GERMLINE match the values
        above from the VCF
        """
        correct_p_germline_values = [
            '-0', '-0', '-0', '-0.029', '-0', '-0',
            '-0.001', '-0.002', '-0', '-0', '-0', '-0'
        ]

        assert correct_p_germline_values == self.vcf_df['P_GERMLINE'].tolist()


    def test_parsed_correct_rpa_values(self):
        """
        Test values read into dataframe for RPA match the values
        above from the VCF
        """
        correct_rpa_values = [
            '11,10', '11,12', '11,13', '3,2', np.nan, '5,4',
            np.nan, np.nan, '4,3', np.nan, '8,9', np.nan
        ]


        assert correct_rpa_values == self.vcf_df['RPA'].tolist()


    def test_parsed_correct_ru_values(self):
        """
        Test values read into dataframe for RU match the values
        above from the VCF
        """
        correct_ru_values = [
            'A', 'A', 'A', 'CGGCGC', np.nan, 'GTG', np.nan,
            np.nan, 'CAC', np.nan, 'T', np.nan
        ]

        assert correct_ru_values == self.vcf_df['RU'].tolist()


    def test_parsed_correct_str_values(self):
        """
        Test values read into dataframe for STR match the values
        above from the VCF, since these are just flags and not key=value
        pairs, these get added to the DF as True if present
        """
        correct_str_values = [
            True, True, True, True, np.nan, True,
            np.nan, np.nan, True, np.nan, True, np.nan
        ]

        assert correct_str_values == self.vcf_df['STR'].tolist()


    def test_parsed_correct_tlod_values(self):
        """
        Test values read into dataframe for TLOD match the values
        above from the VCF
        """
        correct_tlod_values = [
            '36.38', '418.64', '21.15', '3.87', '3878.54', '87.52',
            '5.47', '5.02', '36.68', '2904.48', '19.45', '74.01'
        ]

        assert correct_tlod_values == self.vcf_df['TLOD'].tolist()


    def test_info_column_dropped(self):
        """
        Tests if INFO column dropped from dataframe after it is expanded to
        separate columns
        """
        assert "INFO" not in self.vcf_df.columns.tolist(), 'INFO column not dropped'


class TestCSQ():
    """
    Tests for splitColumns.csq that splits out CSQ string to separate
    columns as defined from header. CSQ fields to split out from test vcf:

    ['SYMBOL', 'VARIANT_CLASS', 'Consequence', 'EXON', 'HGVSc', 'HGVSp',
    'gnomAD_AF', 'gnomADg_AF', 'CADD_PHRED', 'Existing_variation', 'ClinVar',
    'ClinVar_CLNDN', 'ClinVar_CLNSIG', 'COSMIC', 'Feature']

    Data transforms as:

    5 variants with multiple comma separated CSQ annotation as:

    CSQ=LEPR|SNV|intron_variant||NM_001003679.3:c.370+16G>T||0.9999|...
    CSQ=LEPR|SNV|missense_variant|6/20|NM_001003679.3:c.668A>G|...
    CSQ=LEPR|SNV|intron_variant||NM_001003679.3:c.850-50A>C||0.4683|...
    CSQ=LEPR|SNV|synonymous_variant|9/20|NM_001003679.3:c.1029T>C|...
    CSQ=CEP170|SNV|upstream_gene_variant||||0.5003|4.16591e-01|16.36|...


    Multiple transcripts split to separate rows if present:

    LEPR|SNV|intron_variant||NM_001003679.3:c.370+16G>T||...
    LEPR|SNV|intron_variant||NM_001003680.3:c.370+16G>T||...
    LEPR|SNV|intron_variant||NM_001198687.2:c.370+16G>T||...
    LEPR|SNV|intron_variant||NM_001198688.1:c.370+16G>T||...
    LEPR|SNV|intron_variant||NM_001198689.2:c.370+16G>T||...2
    LEPR|SNV|intron_variant||NM_002303.6:c.370+16G>T||...


    LEPR|SNV|missense_variant|6/20|NM_001003679.3:c.668A>G|...
    LEPR|SNV|missense_variant|6/20|NM_001003680.3:c.668A>G|...
    LEPR|SNV|missense_variant|5/19|NM_001198687.2:c.668A>G|...
    LEPR|SNV|missense_variant|5/19|NM_001198688.1:c.668A>G|...
    LEPR|SNV|missense_variant|5/19|NM_001198689.2:c.668A>G|.
    LEPR|SNV|missense_variant|6/20|NM_002303.6:c.668A>G|...


    LEPR|SNV|intron_variant||NM_001003679.3:c.850-50A>C||0.4683|...
    LEPR|SNV|intron_variant||NM_001003680.3:c.850-50A>C||0.4683|...
    LEPR|SNV|intron_variant||NM_001198687.2:c.850-50A>C||0.4683|...
    LEPR|SNV|intron_variant||NM_001198688.1:c.850-50A>C||0.4683|...
    LEPR|SNV|intron_variant||NM_001198689.2:c.850-50A>C||0.4683|...
    LEPR|SNV|intron_variant||NM_002303.6:c.850-50A>C||0.4683|...


    LEPR|SNV|synonymous_variant|9/20|NM_001003679.3:c.1029T>C|...
    LEPR|SNV|synonymous_variant|9/20|NM_001003680.3:c.1029T>C|...
    LEPR|SNV|synonymous_variant|8/19|NM_001198687.2:c.1029T>C|...
    LEPR|SNV|synonymous_variant|8/19|NM_001198688.1:c.1029T>C|...
    LEPR|SNV|synonymous_variant|8/19|NM_001198689.2:c.1029T>C|...
    LEPR|SNV|synonymous_variant|9/20|NM_002303.6:c.1029T>C|...


    CEP170|SNV|upstream_gene_variant||||0.5003|4.16591e-01|16.36|...
    CEP170|SNV|upstream_gene_variant||||0.5003|4.16591e-01|16.36|...
    SDCCAG8|SNV|5_prime_UTR_variant|1/20|NM_001350246.2:c.-1159T>A||...
    SDCCAG8|SNV|5_prime_UTR_variant|1/19|NM_001350247.2:c.-1047T>A||...
    SDCCAG8|SNV|5_prime_UTR_variant|1/19|NM_001350248.2:c.-47T>A||...
    SDCCAG8|SNV|5_prime_UTR_variant|1/18|NM_001350249.2:c.-354T>A||...
    SDCCAG8|SNV|5_prime_UTR_variant|1/21|NM_001350251.2:c.-1420T>A||...
    SDCCAG8|SNV|5_prime_UTR_variant|1/18|NM_006642.5:c.-47T>A||0.5003|...
    CEP170|SNV|upstream_gene_variant||||0.5003|4.16591e-01|16.36|...


    Example split out annotation for one variant / transcript:

    SYMBOL | VARIANT_CLASS | Consequence    | EXON | HGVSc \
    ------------------------------------------------------
    LEPR   | SNV           | intron_variant |      | NM_001003679.3:c.850-50A>C \

    HGVSp | gnomAD_AF | gnomADg_AF  | CADD_PHRED | Existing_variation \
    -----------------------------------------------------------------
          | 0.4683    | 4.99514e-01 | 7.575      | rs3762273 \

    ClinVar | ClinVar_CLNDN | ClinVar_CLNSIG | COSMIC | Feature
    -----------------------------------------------------------
    1263741 | not_provided  | Benign         |        | NM_001003679.3

    """
    # read in test vcf with multiple transcript of annotation for each variant
    vcf_df, csq_fields = read_test_vcf(vcf_file="multi_transcript_annotation.vcf")

    # count total number of transcript annotation across all variants
    # number of transcripts will equal number of ',' + 1 for each string
    transcript_count = 0
    tmp_df = pd.DataFrame()
    tmp_df['CSQ'] = vcf_df['INFO'].apply(lambda x: x.split('CSQ=')[-1])

    for idx, row in tmp_df.iterrows():
        transcript_count += row['CSQ'].count(',') + 1

    # call splitColumns.csq() method to split out CSQ column to test against
    split_csq_vcf_df, expanded_vcf_rows = splitColumns.csq(vcf_df, csq_fields)

    # read csq strings direct from vcf to compare against dataframe columns
    csq_strings = read_csq_from_vcf()


    def test_correct_csq(self):
        """
        Test the columns in the dataframe of expanded CSQ match what is read in
        from file
        """
        # dump out all csq columns from dataframe as a '|' joined list
        df_csq_values = self.split_csq_vcf_df[self.csq_fields].agg('|'.join, axis=1).tolist()

        assert df_csq_values == self.csq_strings, (
            "CSQ values from dataframe do not match what is in file"
        )


    def test_correct_total_number_rows_after_csq_explode(self):
        """
        Test that when we call .explode() on CSQ column to split out multiple
        transcripts to separate rows that this results in the correct number
        of rows
        """
        assert len(self.split_csq_vcf_df.index) == self.transcript_count, (
            "Differing total of rows after .explode() on CSQ column"
        )


    def test_all_csq_fields_present_in_dataframe(self):
        """
        Tests that the csq fields parsed from header are all present in the
        resultant dataframe
        """
        assert all(map(
            lambda v: v in self.split_csq_vcf_df.columns.tolist(), self.csq_fields
        )), "CSQ columns parsed from header not present in df after splitColumns.csq()"




if __name__ == "__main__":
    info = TestInfoColumn()
    csq = TestCSQ()
