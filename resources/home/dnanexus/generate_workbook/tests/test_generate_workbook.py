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
        with pytest.raises(ValueError, match="Expected format"):
            self.args_obj.verify_sort_by([invalid_sort])

    @pytest.mark.parametrize("sort", [["CHROM:"], [":True"]])
    def test_sort_by_input_without_col_or_bool_raises_error(self, sort):
        """
        Test that sort_by values missing column name or boolean raise a
        ValueError.
        """
        with pytest.raises(ValueError, match="no column name or bool"):
            self.args_obj.verify_sort_by(sort)

    @pytest.mark.parametrize(
        "sort", [["CHROM:Not_a_bool"], ["CHROM:F"], ["CHROM:0"], ["CHROM:1"]]
    )
    def test_sort_by_input_with_invalid_bool_raises_error(self, sort):
        """
        Test that sort_by values with invalid booleans raise a ValueError.
        """
        with pytest.raises(ValueError, match="Expected 'True'"):
            self.args_obj.verify_sort_by(sort)

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
        assert self.args_obj.verify_sort_by(sort) == expected
