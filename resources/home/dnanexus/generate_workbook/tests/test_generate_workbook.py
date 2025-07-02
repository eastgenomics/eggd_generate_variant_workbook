from pathlib import Path
import pytest
import argparse
from unittest.mock import patch
from generate_workbook import arguments

class TestVerifyImages():
    """
    Methods to test assertions in generate_workbook.arguments.verify_images
    """
    # initialise the arguments object, using __new__ method to circumvent
    # the __init__ method in arguments() that sets up argument parsing and
    # runs checks to step through and test each
    # patch sys.argv to an empty list to be able to call argparse.parse_args()
    with patch("sys.argv", []):
        args_obj = object.__new__(arguments)
        args_obj.args = args_obj.parse_args()

    def test_invalid_image_assertion(self):
        """
        Tests that AssertionError correctly raised for non image files
        being passed to --images
        """
        self.args_obj.args.images = [Path(__file__).absolute()]

        with pytest.raises(AssertionError):
            self.args_obj.verify_images()

        self.args_obj.args.images = None


    def test_differing_images_and_image_sheet_names(self):
        """
        Tests that when a different number of images and sheet names are
        passed, an AssertionError is correctly raised
        """
        self.args_obj.args.images = [Path(__file__).absolute()]
        self.args_obj.args.image_sheets = ['first_sheet', 'another sheet']

        with pytest.raises(AssertionError):
            self.args_obj.verify_images()

        self.args_obj.args.images = None
        self.args_obj.args.image_sheets = None


    def test_differing_images_and_image_sizes(self):
        """
        Tests that when a different number of images and image sizes are
        passed, an AssertionError is correctly raised
        """
        self.args_obj.args.images = [Path(__file__).absolute()]
        self.args_obj.args.image_sheets = ['1920:1080', '1000:500']

        with pytest.raises(AssertionError):
            self.args_obj.verify_images()

        self.args_obj.args.images = None
        self.args_obj.args.image_sheets = None


    def test_invalid_image_sizes(self):
        """
        Test that when invalid image sizes are passed, an AssertionError is
        correctly raised
        """
        self.args_obj.args.image_sizes = ['1920:1080', '500-250', 'test', '']

        with pytest.raises(AssertionError):
            self.args_obj.verify_images()

        self.args_obj.args.image_sheets = None


class TestVerifyColours():
    """
    Tests for generate_workbook.verify_colours to check for valid
    colouring expressions being given
    """
    with patch("sys.argv", []):
        args_obj = object.__new__(arguments)
        args_obj.args = args_obj.parse_args()


    def test_valid_colour_expressions(self):
        """
        Test a range of valid expressions pass the check
        """
        self.args_obj.args.colour = [
            'VF:>=0.9:green',
            'VF:>0.4:red',
            'VF:<0.9&>=0.4:orange',
            'Consequence:=synonymous_variant|=upstream_variant:green'
        ]

        self.args_obj.verify_colours()


    def test_invalid_colour_expression(self):
        """
        Test that the assertion in verify_colours() is correctly
        raised if both & and | are used in the same expression
        """
        self.args_obj.args.colour = ['VF:<0.9&>=0.4|<1:orange']

        with pytest.raises(AssertionError):
            self.args_obj.verify_colours()


class TestDxFileID():
    """
    Tests for the generate_workbook.dx_file_id method.
    Ensures that a DNAnexus file ID is correctly formatted and parsed.
    """

    # File IDs formatted as raw strings to prevent the hyphen "-" being
    # interpreted by argsparse as a new argument
    @pytest.mark.parametrize(
        "dx_id",
        ["notafileid", r"file-J0527G8487XjKgKJXjYqbgP",
         r"file-J0527G8487XjKgKJXjYqbgPFX1",
         r"file-J0527G8487XjKgKJXjYqbgPFfile-J0527G8487XjKgKJXjYqbgPF"]
    )
    def test_invalid_dx_file_id_raises_correct_error(self, dx_id):
        """
        Tests the correct error is raised when invalid DNAnexus file IDs are
        passed to the --m_codes input
        """
        args_obj = object.__new__(arguments)
        expected = "Invalid DNAnexus file ID:"

        args = ['generate_workbook.py', '--m_codes', dx_id]
        with patch("sys.argv", args):
            # Complicated to test raising of argparse.ArgumentTypeError, see
            # https://stackoverflow.com/a/49324489 for explanation
            with pytest.raises(SystemExit) as e:
                args_obj.args = args_obj.parse_args()

        assert isinstance(e.value.__context__, argparse.ArgumentError)
        assert expected in e.value.__context__.message

    def test_valid_dx_file_id_is_returned(self):
        """
        Test valid DNAnexus file IDs are able to be parsed by argparse.
        """
        args_obj = object.__new__(arguments)

        # File IDs formatted as raw strings to prevent the hyphen "-" being
        # interpreted by argparse as a new argument
        valid_dx_id = "file-J0527G8487XjKgKJXjYqbgPF"
        args = ['generate_workbook.py', '--m_codes', valid_dx_id]

        with patch("sys.argv", args):
            args_obj.args = args_obj.parse_args()

        assert args_obj.args.m_codes == "file-J0527G8487XjKgKJXjYqbgPF"


class TestVerifySortBy():
    """
    Tests for the generate_workbook.verify_sort_by method.
    Ensures that --sort_by input is correctly validated and parsed.
    """

    with patch("sys.argv", []):
        args_obj = object.__new__(arguments)
        args_obj.args = args_obj.parse_args()

    def test_sort_by_input_without_colon_raises_error(self):
        """
        Test that a sort_by string without a colon raises a ValueError.
        """
        invalid_sort = "a_column_name 'True'"
        self.args_obj.args.sort_by = invalid_sort
        with pytest.raises(ValueError, match="Expected format"):
            self.args_obj.verify_sort_by()

    @pytest.mark.parametrize(
        "sort", [["CHROM:"], [":True"]]
    )
    def test_sort_by_input_without_col_or_bool_raises_error(self, sort):
        """
        Test that sort_by values missing column name or boolean raise a
        ValueError.
        """
        self.args_obj.args.sort_by = sort
        with pytest.raises(ValueError, match="no column name or bool"):
            self.args_obj.verify_sort_by()

    @pytest.mark.parametrize(
        "sort", [["CHROM:Not_a_bool"], ["CHROM:F"], ["CHROM:0"], ["CHROM:1"]]
    )
    def test_sort_by_input_with_invalid_bool_raises_error(self, sort):
        """
        Test that sort_by values with invalid booleans raise a ValueError.
        """
        self.args_obj.args.sort_by = sort
        with pytest.raises(ValueError, match="Expected 'True'"):
            self.args_obj.verify_sort_by()

    @pytest.mark.parametrize(
        "sort, expected",
        [
            (["CHROM:True"], {"CHROM": True}),
            (["POS:False"], {"POS": False}),
            (
                ["REF:True", "POS:True", "CHROM:False"],
                {"REF": True, "POS": True, "CHROM": False}
            )
        ]
    )
    def test_valid_sort_by_input_is_returned_correctly(self, sort, expected):
        """
        Test that valid sort_by input is correctly parsed into a dictionary.
        """
        self.args_obj.args.sort_by = sort
        self.args_obj.verify_sort_by()
        assert self.args_obj.args.sort_by == expected
