import os
import sys

from sklearn.tree import DecisionTreeRegressor

sys.path.append(os.path.abspath(
    os.path.join(os.path.realpath(__file__), '../../')
))

from utils.utils import determine_delimeter, is_numeric


def test_is_numeric():
    """
    Tests that all forms of numerical value are correctly identified as
    being numeric, and things we don't want to are treat as non-numeric
    """
    numeric_values = ['5', '5.5', '-0.1', '1e-5', '.1']

    non_numeric_values = [
        '5,6', '5, 6', '.', '3/4', '1..1', '1.1.1', '5ee5', '123e', '-e1',
        '.e1', 'e1', '1.1.1', '01', '0001'
    ]

    assert all([is_numeric(x) for x in numeric_values]), (
        'some numeric values do not properly evaluate to being numeric'
    )

    assert all([not is_numeric(x) for x in non_numeric_values]), (
        'some non-numeric values wrongly evaluate to being numeric'
    )


class TestDetermineDelimeter():
    """
    Tests for utils.determine_delimeter for detecting correct delimeter
    """
    @staticmethod
    def test_comma():
        delimeter = determine_delimeter('this,is,a,comma,separated,string')

        assert delimeter == ',', 'failed to correctly identify comma delimeter'

    @staticmethod
    def test_semicolon():
        delimeter = determine_delimeter('this;is;a;semi colon;separated;string')

        assert delimeter == ';', (
            'failed to correctly identify semi colon delimeter'
        )

    @staticmethod
    def test_tab():
        delimeter = determine_delimeter('this\tis\ta\tab\tseparated\tstring')

        assert delimeter == '\t', 'failed to correctly identify tab delimeter'

    @staticmethod
    def test_space():
        delimeter = determine_delimeter('this is a space separated string')

        assert delimeter == ' ', 'failed to correctly identify space delimeter'

    @staticmethod
    def test_mixed():
        delimeter = determine_delimeter(
            '#this;is;a string with a.mix, of characters\nthat\tshould\tbe'
            '\tidentified\nas\ttab\tdelimited\tbecause\tit\thas\na\tweird'
            '\theader\tline'
        )

        assert delimeter == '\t', 'failed to correctly identify space delimeter'


