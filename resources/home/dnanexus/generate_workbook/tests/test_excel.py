import os
import sys
import pytest
from pathlib import Path
from dxpy.bindings.dxfile_functions import open_dxfile
import openpyxl

sys.path.append(os.path.abspath(
    os.path.join(os.path.realpath(__file__), '../../')
))

from utils.excel import excel


# fixtures
@pytest.fixture
def mocked_excel_file(mocker):
    mock_args = mocker.Mock()
    mock_args.output = Path(os.getcwd() + "tmp.xlsx")
    mock_args.sheets = ["sheet1", "sheet2", "sheet3"]
    mock_vcfs = mocker.Mock()
    mock_additional_files = mocker.Mock()
    mock_refs = mocker.Mock()
    return excel(mock_args, mock_vcfs, mock_additional_files, mock_refs)

@pytest.fixture
def locked_worksheet(mocked_excel_file):
    wb = openpyxl.Workbook()
    ws = wb.active
    mocked_excel_file.lock_sheet(ws)
    return ws

@pytest.fixture
def unlocked_worksheet():
    wb = openpyxl.Workbook()
    ws = wb.active
    return ws

# lock_sheet()
def test_attributes_are_set_correctly(locked_worksheet):
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

# unlock_specified_cells()
def test_unlock_specified_cells(mocked_excel_file, locked_worksheet):
    target_cells = [locked_worksheet[cell] for cell in ["A1", "B2", "C3"]]
    # ensure cells are currently locked
    cells_are_locked = [cell.protection.locked for cell in target_cells]
    assert all(cells_are_locked)
    # unlock cells and ensure they are unlocked
    mocked_excel_file.unlock_specified_cells(locked_worksheet, ["A1", "B2", "C3"])
    cells_are_unlocked = [cell.protection.locked is False for cell in target_cells]
    assert all(cells_are_unlocked)
    # ensure unrelated cell remains unlocked
    assert locked_worksheet["D4"].protection.locked

# unlocked_region()
def test_unlock_region(mocked_excel_file, locked_worksheet):
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

# get_cells_in_columns()
def test_get_cells_in_columns_returns_cells(mocked_excel_file, unlocked_worksheet):
    # set header for unlocked worksheet
    unlocked_worksheet["A1"].value = "COL_1"
    unlocked_worksheet["B1"].value = "COL_2"
    returned_cells = mocked_excel_file.get_cells_in_columns(
            sheet=unlocked_worksheet, 
            cols=["COL_1", "COL_2"], 
            num_rows=2)
    assert returned_cells == ["A2", "A3", "B2", "B3"]

def test_get_cells_missing_column_error(mocked_excel_file, unlocked_worksheet):
    with pytest.raises(RuntimeError, match=".*The column.*"):
        mocked_excel_file.get_cells_in_columns(
                sheet=unlocked_worksheet,
                cols=["I_DONT_EXIST"],
                num_rows=1)

# store_list_in_sheet()
def test_store_list_in_sheet(mocked_excel_file):
    mocked_excel_file.store_list_in_sheet(
            values=["Waterloo", "Gimme Gimme Gimme", "Lay All Of Your Love On Me"],
            sheet_name="sheet1",
            col="A")
    wb = mocked_excel_file.workbook
    ws = wb.active
    assert ws["A1"].value == "Waterloo"
    assert ws["A2"].value == "Gimme Gimme Gimme"
    assert ws["A3"].value == "Lay All Of Your Love On Me"

# read_m_codes_file
def test_success_upon_compliant_mcodes(mocked_excel_file, monkeypatch):
    class CompliantMFile:
        def __init__(self):
            pass

        def read(self):
            return "M123"

    def return_compliant_file(*args, **kwargs):
        return CompliantMFile() 

    monkeypatch.setattr("utils.excel.open_dxfile", return_compliant_file)
    mocked_excel_file.read_m_codes_file()

def test_exception_upon_noncompliant_mcodes(mocked_excel_file, monkeypatch):
    class NonCompliantMFile:
        def __init__(self):
            pass

        def read(self):
            return "hi\nthis\nshouldnt\nwork"

    def return_noncompliant_file(*args, **kwargs):
        return NonCompliantMFile() 

    monkeypatch.setattr("utils.excel.open_dxfile", return_noncompliant_file)
    with pytest.raises(ValueError, match=r".*M-codes file not formatted correctly.*"):
        mocked_excel_file.read_m_codes_file()

# str_to_dropdown()
def test_str_to_drop_down_options_too_long(mocked_excel_file):
    horrid_string = "".join(["A" for i in range(300)])
    with pytest.raises(ValueError, match=r".*>256 characters long.*"):
        mocked_excel_file.str_to_drop_down(
                dropdown_options=horrid_string,
                prompt="A",
                title="A",
                sheet="A",
                cells=["A1"])

# format_cells_as_percentage
def test_cells_formatted_as_perc(mocked_excel_file):
    wb = openpyxl.Workbook()
    wb.save(mocked_excel_file.args.output)
    ws = wb.worksheets[0]
    ws["A1"] = "41.23"
    mocked_excel_file.format_cells_as_percentage(ws, cells = ["A1"])
    assert ws["A1"].number_format == "0.0%"
