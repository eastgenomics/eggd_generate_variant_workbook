import os
import sys


sys.path.append(os.path.abspath(
    os.path.join(os.path.realpath(__file__), '../../')
))

from utils.utils import determine_delimiter, is_numeric


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


class TestDeterminedelimiter():
    """
    Tests for utils.determine_delimiter for detecting correct delimiter
    """
    @staticmethod
    def test_comma():
        delimiter = determine_delimiter('this,is,a,comma,separated,string', [])

        assert delimiter == ',', 'failed to correctly identify comma delimiter'

    @staticmethod
    def test_semicolon():
        delimiter = determine_delimiter('this;is;a;semi colon;separated;string', [])

        assert delimiter == ';', (
            'failed to correctly identify semi colon delimiter'
        )

    @staticmethod
    def test_tab():
        delimiter = determine_delimiter('this\tis\ta\tab\tseparated\tstring', [])

        assert delimiter == '\t', 'failed to correctly identify tab delimiter'

    @staticmethod
    def test_space():
        delimiter = determine_delimiter('this is a space separated string', [])

        assert delimiter == ' ', 'failed to correctly identify space delimiter'

    @staticmethod
    def test_mixed():
        delimiter = determine_delimiter(
            '#this;is;a string with a.mix, of characters\nthat\tshould\tbe'
            '\tidentified\nas\ttab\tdelimited\tbecause\tit\thas\na\tweird'
            '\theader\tline', []
        )

        assert delimiter == '\t', 'failed to correctly identify space delimiter'

    @staticmethod
    def test_tsv_suffix():
        delimiter = determine_delimiter('tsvFileStringWithNodelimiters', ['.tsv'])

        assert delimiter == '\t', 'failed to correctly identify space delimiter'

    @staticmethod
    def test_csv_suffix():
        delimiter = determine_delimiter('csvFileStringWithNodelimiters', ['.csv'])

        assert delimiter == ',', 'failed to correctly identify space delimiter'