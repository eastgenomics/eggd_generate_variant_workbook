import argparse
import os
import sys

import pytest

sys.path.append(os.path.abspath(os.path.join(os.path.realpath(__file__), '../../')))

from src.vcf import vcf


header_test_vcf = "test_data/header.vcf"

def test_reading_header(vcf_handler):
    header, columns = vcf.parse_header(header_test_vcf)

    print(header)



if __name__=="__main__":
    vcf_handler = vcf(argparse.Namespace)
    test_reading_header(vcf_handler)