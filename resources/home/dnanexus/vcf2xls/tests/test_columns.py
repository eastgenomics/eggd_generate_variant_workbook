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

# test data vcf
columns_vcf = os.path.join(TEST_DATA_DIR, "column_methods_test.vcf")

# initialise vcf class
vcf_handler = vcf(argparse.Namespace)
vcf_df, csq_fields = vcf_handler.read(columns_vcf)

print(vcf_df)



    # # read in header from our test vcf, call functions to parse reference and
    # # csq feilds from read in header
    # header, columns = vcf_handler.parse_header(header_test_vcf)
    # vcf_handler.parse_reference(header)
    # vcf_handler.parse_csq_fields(header)


class TestColumnInfo():
    """
    Tests for columns.info() that splits out key value pairs from INFO column to
    separate columns in dataframe
    """

    def test_() -> None:
        """
        columns.info()
        """
        pass

# if __name__ == "__main__":
    