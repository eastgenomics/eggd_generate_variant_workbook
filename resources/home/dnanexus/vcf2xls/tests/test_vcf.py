import argparse
import os
from pathlib import Path
import sys

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
    # csq feilds from read in header
    header, columns = vcf_handler.parse_header(header_test_vcf)
    vcf_handler.parse_reference(header)
    vcf_handler.parse_csq_fields(header)


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


    def test_parse_csq_fields(self):
        """
        Tests the read in csq fields from header are correct
        """
        correct_csq_fields = [
            "Allele", "Gene", "HGNC", "RefSeq", "Feature", "Consequence",
            "cDNA_position", "Protein_position", "Amino_acids",
            "Existing_variation", "SIFT", "PolyPhen", "HGVSc"
        ]

        assert vcf_handler.csq_fields == correct_csq_fields



if __name__=="__main__":
    header = TestHeader()
    header.test_column_names()

    print(vcf_handler.refs)