import os
import sys
import pytest
from pathlib import Path
import openpyxl

sys.path.append(os.path.abspath(
    os.path.join(os.path.realpath(__file__), '../../')
))

from utils.excel import excel


# Fixtures
@pytest.fixture
def mocked_excel_file(mocker):
    """
    Fixture for tests where a utils.excel.excel object is required
    """
    mock_args = mocker.Mock()
    mock_args.output = Path(os.getcwd() + "/tmp.xlsx")
    mock_args.sheets = ["sheet1", "sheet2", "sheet3"]
    mock_vcfs = mocker.Mock()
    mock_additional_files = mocker.Mock()
    mock_refs = mocker.Mock()
    return excel(mock_args, mock_vcfs, mock_additional_files, mock_refs)

@pytest.fixture
def locked_worksheet(mocked_excel_file):
    """
    Fixture for tests where a locked openpyxl.Workbook object is required
    """
    wb = openpyxl.Workbook()
    ws = wb.active
    mocked_excel_file.lock_sheet(ws)
    return ws

@pytest.fixture
def unlocked_worksheet():
    """
    Fixture for tests where an unlocked openpyxl.Workbook object is required
    """
    wb = openpyxl.Workbook()
    ws = wb.active
    return ws

class TestLockSheet:
    """
    Collection of test cases for utils.excel.excel.lock_sheet
    """
    def test_attributes_are_set_correctly(self, locked_worksheet):
        """
        Tests that attributes are correctly set when worksheet is locked
        """
        expected = {
                "sheet": "1",
                "autoFilter": "0",
                "formatColumns": "0",
                "formatRows": "0",
                "formatCells": "0",
                "password": locked_worksheet.protection.password
            }
        actual = dict(locked_worksheet.protection)
        assert expected.items() <= actual.items()

class TestUnlockSpecifiedCells:
    """
    Collection of test cases for utils.excel.excel.unlock_specified_cells
    """
    def test_unlock_specified_cells(self, mocked_excel_file, locked_worksheet):
        """
        Tests that cells targeted by utils.excel.excel.unlock_specified_cells are unlocked without
        unlocking unrelated cells
        """
        target_cells = [locked_worksheet[cell] for cell in ["A1", "B2", "C3"]]
        # ensure cells are currently locked
        cells_are_locked = [cell.protection.locked for cell in target_cells]
        assert all(cells_are_locked)
        # unlock cells and ensure they are unlocked
        mocked_excel_file.unlock_specified_cells(locked_worksheet, ["A1", "B2", "C3"])
        cells_are_unlocked = [cell.protection.locked is False for cell in target_cells]
        assert all(cells_are_unlocked)
        # ensure unrelated cell remains locked
        assert locked_worksheet["D4"].protection.locked

class TestUnlockRegion:
    """
    Collection of test cases for utils.excel.excel.unlock_region
    """
    def test_unlock_region(self, mocked_excel_file, locked_worksheet):
        """
        Tests that cells targeted by utils.excel.excel.unlock_region are unlocked without
        unlocking unrelated cells
        """
        target_cells = [locked_worksheet[cell] for cell in ["A1", "A2", "B1", "B2"]]
        # ensure cells are currently locked
        cells_are_locked = [cell.protection.locked for cell in target_cells]
        assert all(cells_are_locked)
        # unlock cells and ensure they are unlocked
        mocked_excel_file.unlock_region(
                ws=locked_worksheet,
                start_row=1,
                start_col=1,
                unlock_row_num=2,
                unlock_col_num=2
                )
        cells_are_unlocked = [cell.protection.locked is False for cell in target_cells]
        assert all(cells_are_unlocked)
        # ensure unrelated cell remains unlocked
        assert locked_worksheet["B3"].protection.locked

class TestGetCellsInColumns:
    """
    Collection of test cases for utils.excel.excel.get_cells_in_columns
    """
    def test_get_cells_in_columns_returns_cells(self, mocked_excel_file, unlocked_worksheet):
        """
        Test that utils.excel.excel.get_cells_in_columns only returns the targeted cells
        """
        # set header for unlocked worksheet
        unlocked_worksheet["A1"].value = "COL_1"
        unlocked_worksheet["B1"].value = "COL_2"
        returned_cells = mocked_excel_file.get_cells_in_columns(
                sheet=unlocked_worksheet, 
                cols=["COL_1", "COL_2"], 
                num_rows=2)
        assert returned_cells == ["A2", "A3", "B2", "B3"]

    def test_get_cells_missing_column_error(self, mocked_excel_file, unlocked_worksheet):
        """
        Test that utils.excel.excel.get_cells_in_columns throws a `RunTimeError` when the target
        column doesn't exist
        """
        with pytest.raises(RuntimeError, match=".*The column.*"):
            mocked_excel_file.get_cells_in_columns(
                    sheet=unlocked_worksheet,
                    cols=["I_DONT_EXIST"],
                    num_rows=1)

class TestStoreListInSheet:
    """
    Collection of test cases for utils.excel.excel.store_list_in_sheet
    """
    def test_store_list_in_sheet(self, mocked_excel_file):
        """
        Test that excel.store_list_in_sheet stores an input list into
        the correct column
        """
        mocked_excel_file.store_list_in_sheet(
                values=["Waterloo", "Gimme Gimme Gimme", "Lay All Of Your Love On Me"],
                sheet_name="sheet1",
                col="A")
        wb = mocked_excel_file.workbook
        ws = wb.active
        assert ws["A1"].value == "Waterloo"
        assert ws["A2"].value == "Gimme Gimme Gimme"
        assert ws["A3"].value == "Lay All Of Your Love On Me"


class TestReadMCodesFile:
    """
    Collection of test cases for utils.excel.excel.read_m_codes_file
    """

    def test_success_upon_compliant_mcodes(self, tmp_path):
        """
        Test that a file containing m-codes that conform to the regex
        specification are read without throwing an error when read by
        utils.excel.excel.read_m_codes_file
        """

        compliant_content = "M1\nM2\nM3\n"
        file_path = tmp_path / "compliant_mcodes.txt"
        file_path.write_text(compliant_content, encoding="utf-8")

        result = excel.read_m_codes_file(file_path)

        assert result == ["M1", "M2", "M3"]

    def test_exception_upon_noncompliant_mcodes(self, tmp_path):
        """
        Test that a file containing m-codes that do not conform to the regex
        specification cause a `ValueError` to be thrown when it is read by
        utils.excel.excel.read_m_codes_file
        """
        noncompliant_content = "hi\nthis\nshouldnt\nwork\n"
        file_path = tmp_path / "noncompliant_mcodes.txt"
        file_path.write_text(noncompliant_content, encoding="utf-8")

        with pytest.raises(ValueError,
                           match=r".*M-codes file not formatted correctly.*"):
            excel.read_m_codes_file(file_path)


class TestStrToDropdown:
    """
    Collection of test cases for utils.excel.excel.str_to_dropdown
    """
    def test_str_to_drop_down_options_too_long(self, mocked_excel_file):
        """
        Test that utils.excel.excel.str_to_dropdown throws a `ValueError` when the input string
        is longer than 256 characters
        """
        horrid_string = "".join(["A" for i in range(300)])
        with pytest.raises(ValueError, match=r".*>256 characters long.*"):
            mocked_excel_file.str_to_drop_down(
                    dropdown_options=horrid_string,
                    prompt="A",
                    title="A",
                    sheet="A",
                    cells=["A1"])

class TestCellsFormattedAsPerc:
    """
    Collection of test cases for utils.excel.excel.cells_formatted_as_perc
    """
    def test_cells_formatted_as_perc(self, mocked_excel_file):
        """
        Test that the openpyxl.Workbook.number_format attribute is set to "0.0%" when
        utils.excel.excel.format_cells_as_percentages is executed
        """
        wb = openpyxl.Workbook()
        wb.save(mocked_excel_file.args.output)
        ws = wb.worksheets[0]
        ws["A1"] = "41.23"
        mocked_excel_file.format_cells_as_percentage(ws, cells=["A1"])
        assert ws["A1"].number_format == "0.0%"


class TestReadCancerPanelsFile:
    """
    Collection of test cases for utils.excel.excel.read_cancer_panels_file
    """
    def test_valid_file_read_and_deduplicated(self, tmp_path):
        """
        Test that a valid file is read with whitespace stripped, blank lines
        skipped and duplicate cancer type / gene pairs removed
        """
        content = (
            "cancer_type\tgene\n"
            "Colorectal\tKRAS\n"
            " Colorectal \t NRAS \n"
            "\n"
            "Colorectal\tKRAS\n"
            "Lung\tEGFR\n"
        )
        file_path = tmp_path / "panels.tsv"
        file_path.write_text(content, encoding="utf-8")

        result = excel.read_cancer_panels_file(file_path)

        assert result.values.tolist() == [
            ["Colorectal", "KRAS"],
            ["Colorectal", "NRAS"],
            ["Lung", "EGFR"]
        ]

    def test_incorrect_header_raises_error(self, tmp_path):
        """
        Test that a file without the expected header raises a ValueError
        """
        file_path = tmp_path / "panels.tsv"
        file_path.write_text("Colorectal\tKRAS\n", encoding="utf-8")

        with pytest.raises(ValueError, match=r".*Expected header.*"):
            excel.read_cancer_panels_file(file_path)

    @pytest.mark.parametrize(
        "line", ["Colorectal\t\n", "\tKRAS\n", "Colorectal\n",
                 "Colorectal\tKRAS\textra\n"]
    )
    def test_malformed_line_raises_error(self, tmp_path, line):
        """
        Test that a line with an empty value or wrong number of columns
        raises a ValueError
        """
        file_path = tmp_path / "panels.tsv"
        file_path.write_text(f"cancer_type\tgene\n{line}", encoding="utf-8")

        with pytest.raises(ValueError, match=r".*found in line 2.*"):
            excel.read_cancer_panels_file(file_path)

    def test_header_only_raises_error(self, tmp_path):
        """
        Test that a file with no cancer types / genes raises a ValueError
        """
        file_path = tmp_path / "panels.tsv"
        file_path.write_text("cancer_type\tgene\n", encoding="utf-8")

        with pytest.raises(ValueError, match=r".*No cancer types.*"):
            excel.read_cancer_panels_file(file_path)


class TestAtlasCancerPanelDropDown:
    """
    Collection of test cases for
    utils.excel.excel.atlas_cancer_panel_drop_down
    """
    def test_panel_sheets_and_drop_down_added(self, mocked_excel_file, tmp_path):
        """
        Test that hidden sheets are created holding the panels and cancer
        types, and the summary sheet has a drop-down referencing the
        cancer types
        """
        file_path = tmp_path / "panels.tsv"
        file_path.write_text(
            "cancer_type\tgene\nColorectal\tKRAS\nColorectal\tNRAS\n"
            "Lung\tEGFR\n",
            encoding="utf-8"
        )
        mocked_excel_file.args.cancer_panels = file_path
        mocked_excel_file.summary = mocked_excel_file.workbook.create_sheet(
            'summary')

        mocked_excel_file.atlas_cancer_panel_drop_down()

        wb = mocked_excel_file.workbook
        panels = wb['cancer_panels']
        types = wb['cancer_types']

        assert panels.sheet_state == 'hidden'
        assert types.sheet_state == 'hidden'
        assert [list(x) for x in panels.values] == [
            ["Colorectal", "KRAS"], ["Colorectal", "NRAS"], ["Lung", "EGFR"]
        ]
        assert [x[0] for x in types.values] == ["Colorectal", "Lung"]

        assert mocked_excel_file.summary["A8"].value == "Cancer type panel"
        validations = mocked_excel_file.summary.data_validations.dataValidation
        assert len(validations) == 1
        assert validations[0].formula1 == "='cancer_types'!A1:A2"
        assert str(validations[0].sqref) == "B8"


class TestAddInPanelColumn:
    """
    Collection of test cases for utils.excel.excel.add_in_panel_column
    """
    def test_in_panel_column_added(self, mocked_excel_file, unlocked_worksheet):
        """
        Test that the In panel column is added after the last column with a
        formula per row referencing that row's gene cell
        """
        mocked_excel_file.args.panel_gene_column = "SYMBOL"
        unlocked_worksheet.append(["CHROM", "SYMBOL", "VF"])
        unlocked_worksheet.append(["1", "KRAS", 0.5])
        unlocked_worksheet.append(["2", "EGFR", 0.4])

        mocked_excel_file.add_in_panel_column(unlocked_worksheet, 2)

        assert unlocked_worksheet["D1"].value == "In panel"
        assert unlocked_worksheet["D2"].value == (
            "=IF('summary'!$B$8=\"\",\"\",IF(COUNTIFS('cancer_panels'!$A:$A,"
            "'summary'!$B$8,'cancer_panels'!$B:$B,B2)>0,\"Yes\",\"No\"))"
        )
        assert "B3)" in unlocked_worksheet["D3"].value
        assert unlocked_worksheet["D4"].value is None

    def test_missing_gene_column_skipped(
        self, mocked_excel_file, unlocked_worksheet
    ):
        """
        Test that no column is added when the gene column is not present
        """
        mocked_excel_file.args.panel_gene_column = "SYMBOL"
        unlocked_worksheet.append(["CHROM", "POS"])
        unlocked_worksheet.append(["1", 100])

        mocked_excel_file.add_in_panel_column(unlocked_worksheet, 1)

        assert unlocked_worksheet.max_column == 2
