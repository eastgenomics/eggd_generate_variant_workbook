from pathlib import Path
import pytest
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
        self.args_obj.args.colours = [
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
        self.args_obj.args.colours = ['VF:<0.9&>=0.4|<1:orange']

        with pytest.raises(AssertionError):
            self.args_obj.verify_colours()
