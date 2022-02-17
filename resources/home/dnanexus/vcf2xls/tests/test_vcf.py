import argparse
import os
import sys

import pytest

sys.path.append(os.path.abspath(os.path.join(os.path.realpath(__file__), '../../')))

from utils.vcf import vcf
from utils.columns import splitColumns

sys.path.append('..')
print(sys.path)
header_test_vcf = "./tests/test_data/header.vcf"

def test_reading_header(vcf_handler):
    header, columns = vcf_handler.parse_header(header_test_vcf)

    print(header)



if __name__=="__main__":
    vcf_handler = vcf(argparse.Namespace)
    test_reading_header(vcf_handler)