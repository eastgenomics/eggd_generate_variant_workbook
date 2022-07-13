import os
import sys

sys.path.append(os.path.abspath(
    os.path.join(os.path.realpath(__file__), '../../')
))

from utils.utils import is_numeric


def test_is_numeric():
    """
    Tests that all forms of numerical value are correctly identified as
    being numeric
    """
    numeric_values = ['5', '5.5', '-0.1', '1e-5']
    non_numeric_values = ['5,6', '5, 6', '.', '3/4', ]

    assert all([is_numeric(x) for x in numeric_values]), (
        'some numeric values do not properly evaluate to being numeric'
    )

    assert all([not is_numeric(x) for x in non_numeric_values]), (
        'some noin-numeric values wrongly evaluate to being numeric'
    )
