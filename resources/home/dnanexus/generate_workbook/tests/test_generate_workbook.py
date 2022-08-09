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
