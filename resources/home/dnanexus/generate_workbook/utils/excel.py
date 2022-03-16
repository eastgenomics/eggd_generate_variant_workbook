from email import header
from pathlib import Path
from string import ascii_uppercase as uppercase
import sys
from timeit import default_timer as timer

import Levenshtein as levenshtein
from openpyxl import load_workbook
from openpyxl.styles import Border, DEFAULT_FONT, Font, Side
from openpyxl.styles.fills import PatternFill
import pandas as pd
from pandas.api.types import is_numeric_dtype


# openpyxl style settings
THIN = Side(border_style="thin", color="000000")
THIN_BORDER = Border(left=THIN, right=THIN, top=THIN, bottom=THIN)
DEFAULT_FONT.name = 'Calibri'


class excel():
    """
    Functions for wrangling variant data into spreadsheet formatting and
    writing output file

    Attributes
    ----------
    args : argparse.Namespace
        arguments passed from command line
    vcfs : list of pd.DataFrame
        list of dataframes formatted to write to file from vcf() methods
    refs : list
        list of reference names parsed from vcf headers
    writer : pandas.io.excel._openpyxl.OpenpyxlWriter
        writer object for writing Excel data to file
    workbook : openpyxl.workbook.workbook.Workbook
        openpyxl workbook object for interacting with per-sheet writing and
        formatting of output Excel file

    Outputs
    -------
    {args.output}.xlsx : file
        Excel file with variants written to, name passed from command line or
        inferred from input vcf name if not specified
    """
    def __init__(self, args, vcfs, refs) -> None:
        print(f"Writing to output file: {Path(args.output).absolute()}")
        self.args = args
        self.vcfs = vcfs
        self.refs = refs
        self.writer = pd.ExcelWriter(args.output, engine='openpyxl')
        self.workbook = self.writer.book


    def generate(self) -> None:
        """
        Calls all methods in excel() to generate output file
        """
        self.write_summary()
        self.write_variants()

        self.workbook.save(self.args.output)
        print('Done!')


    def write_summary(self) -> None:
        """
        Write summary sheet to excel file
        """
        print('Writing summary sheet')
        if self.args.summary == 'dias':
            # generate summary sheet in format for RD/dias
            self.summary = self.workbook.create_sheet('summary')
            self.dias_summary()


    def dias_summary(self) -> None:
        """
        Write summary sheet in format for RD group, adds the following info:
            - sample ID, panel(s), run IDs etc.
            - formatted tables for them to fill in reporting
        """
        # write titles for summary values
        self.summary.cell(1, 1).value = "Sample ID:"
        self.summary.cell(1, 5).value = "Clinical Indication(s):"
        self.summary.cell(2, 5).value = "Panel(s):"
        self.summary.cell(34, 1).value = "Workflow:"
        self.summary.cell(35, 1).value = "Workflow ID:"
        self.summary.cell(36, 1).value = "VCF Analysis:"
        self.summary.cell(37, 1).value = "VCF Analysis ID:"
        self.summary.cell(39, 1).value = "Total records:"

        # write summary values
        self.summary.cell(1, 2).value = self.args.sample
        self.summary.cell(1, 6).value = self.args.clinical_indication
        self.summary.cell(2, 6).value = self.args.panel
        self.summary.cell(34, 2).value = self.args.workflow[0]
        self.summary.cell(35, 2).value = self.args.workflow[1]
        self.summary.cell(36, 2).value = self.args.analysis[0]
        self.summary.cell(37, 2).value = self.args.analysis[1]

        # write total rows in each sheet
        count = 40
        for sheet, vcf in zip(self.args.sheets, self.vcfs):
            self.summary.cell(count, 2).value = sheet
            self.summary.cell(count, 3).value = len(vcf.index)
            self.summary[f"A{count}"].font = Font(bold=True, name=DEFAULT_FONT.name)
            count += 1


        # write args passed to script to generate report
        self.summary.cell(46, 1).value = "Filters applied:"
        if self.args.filter:
            self.summary.cell(46, 2).value = self.args.filter
        else:
            self.summary.cell(46, 2).value = "None"

        # write center reporting section tables
        self.summary.cell(9, 2).value = "Phenotype:"

        self.summary.cell(16, 2).value = "Panels"
        self.summary.cell(16, 3).value = "Excel file"
        self.summary.cell(16, 4).value = "Comments"
        self.summary.cell(16, 6).value = "Analysis by"
        self.summary.cell(16, 7).value = "Date"
        self.summary.cell(16, 8).value = "Checked by"
        self.summary.cell(16, 9).value = "Date"

        self.summary.cell(21, 2).value = "Sanger sequencing confirmation"
        self.summary.cell(22, 2).value = "Gene"
        self.summary.cell(22, 3).value = "NM_#"
        self.summary.cell(22, 4).value = "Coordinate"
        self.summary.cell(22, 5).value = "cDNA"
        self.summary.cell(22, 6).value = "Protein change"
        self.summary.cell(22, 7).value = "WS#"
        self.summary.cell(22, 8).value = "Confirmed (Y/N)"

        self.summary.cell(28, 2).value = "GEM comments summary"
        self.summary.cell(28, 4).value = "Date"

        # merge some title columns that have longer text
        self.summary.merge_cells(
            start_row=9, end_row=9, start_column=2, end_column=5)
        self.summary.merge_cells(
            start_row=21, end_row=21, start_column=2, end_column=8)
        self.summary.merge_cells(
            start_row=16, end_row=16, start_column=4, end_column=5)
        self.summary.merge_cells(
            start_row=28, end_row=28, start_column=2, end_column=3)
        self.summary.merge_cells(
            start_row=28, end_row=28, start_column=4, end_column=6)

        # set titles to bold
        title_cells = [
            "A1", "A34", "A35", "A36", "A37", "A39", "B1", "B9", "B16",
            "B21", "B22", "B28", "B34", "B35", "B36", "B37",
            "C16", "C22", "D16", "D22", "D28", "E1", "E2", "E22",
            "F1", "F2", "F16", "F22", "G16", "G22", "H16", "H22", "I16"
        ]
        for cell in title_cells:
            self.summary[cell].font = Font(bold=True, name=DEFAULT_FONT.name)

        # set column widths for readability
        self.summary.column_dimensions['A'].width = 18
        self.summary.column_dimensions['B'].width = 13
        self.summary.column_dimensions['C'].width = 13
        self.summary.column_dimensions['D'].width = 13
        self.summary.column_dimensions['E'].width = 18
        self.summary.column_dimensions['F'].width = 16
        self.summary.column_dimensions['G'].width = 16
        self.summary.column_dimensions['H'].width = 16

        # colour title cells
        blueFill = PatternFill(patternType="solid", start_color="0CABA8")

        colour_cells = [
            "B9", "B16", "B21", "B22", "B28", "C16", "C22", "D16", "D22",
            "D28", "E22", "F16", "F22", "G16", "G22", "H16", "H22", "I16"
        ]
        for cell in colour_cells:
            self.summary[cell].fill = blueFill

        # set borders around table areas
        row_ranges = [
            'B9:E9', 'B10:E10', 'B11:E11', 'B12:E12', 'B13:E13',
            'B16:I16', 'B17:I17', 'B18:I18',
            'B21:H21', 'B22:H22', 'B23:H23', 'B24:H24', 'B25:H25',
            'B28:F28', 'B29:F29', 'B30:F30', 'B31:F31', 'B32:F32'
        ]
        for row in row_ranges:
            for cells in self.summary[row]:
                for cell in cells:
                    cell.border = THIN_BORDER


    def write_variants(self) -> None:
        """
        Writes all variants from dataframe(s) to sheet(s) specified in
        self.args.sheets.

        If sheet names not specified, these will be set as "variant" where one
        dataframe is being written, or the vcf filename prefix if there are
        more than one dataframes to write
        """
        total_rows = sum([len(x) for x in self.vcfs])
        print(f"\nWriting total of {total_rows} rows to output xlsx file")
        if total_rows > 5000:
            print(
                "Writing many rows to Excel is slow, "
                "this may take a few minutes..."
            )

        with self.writer:
            # add variants
            for sheet, vcf in zip(self.args.sheets, self.vcfs):
                sheet_no = self.args.sheets.index(sheet) + 1
                print(
                    f"\nWriting {len(vcf)} rows to {sheet} sheet "
                    f"({sheet_no}/{len(self.args.sheets)})"
                )

                # timing how long it takes to write because its slow
                start = timer()
                vcf.to_excel(
                    self.writer, sheet_name=sheet,
                    index=False, float_format="%.3f"
                )

                curr_worksheet = self.writer.sheets[sheet]
                self.set_widths(curr_worksheet, vcf.columns)
                self.set_font(curr_worksheet)
                self.colour_hyperlinks(curr_worksheet)

                self.workbook.save(self.args.output)
                end = timer()
                print(
                    f"Writing to Excel for sheet {sheet} took "
                    f"{round(end - start)}s"
                )

                # read back the written sheet to check its written correctly
                self.check_written_sheets(vcf, sheet)


    def check_written_sheets(self, vcf, sheet) -> None:
        """"
        Check that the written sheet is exactly the same as the dataframe
        that was meant to be written

        Parameters
        ----------
        vcf : pd.DataFrame
            dataframe of sheet that was written to file
        sheet : str
            sheet name to read from file

        Raises
        ------
        AssertionError
            Raised when data written to sheet does not match the dataframe
        """
        print(f"\nVerifying data written to file for {sheet} sheet\n")

        # read in written sheet using openpyxl to deal with Excel oddities
        sheet_data = load_workbook(filename=self.args.output)[sheet]
        written_sheet = pd.DataFrame(
            sheet_data.values, columns=vcf.columns.tolist())
        written_sheet = written_sheet.iloc[1:]  # drop header on first row

        # openpyxl read sets NaNs to None, so match it
        vcf.fillna('None', inplace=True)

        # set all columns of both dfs to strings
        vcf = vcf.astype(str)
        written_sheet = written_sheet.astype(str).reset_index(drop=True)

        # floats with trailing zero seem inconsistent when writing to file,
        # since we've cast everything to a string strip for ease of comparing
        vcf = vcf.applymap(
            lambda x: x.replace('.0', '') if x.endswith('.0') else x
        )
        written_sheet = written_sheet.applymap(
            lambda x: x.replace('.0', '') if x.endswith('.0') else x
        )

        print("Checking")
        assert vcf.equals(written_sheet), (
            f"Written data for sheet: {sheet} does not seem to match the "
            "dataframe to be written"
        )


    def set_font(self, worksheet) -> None:
        """
        Set font to all cells in variant sheet to Calibri

        Default is Times New Roman and it is ugly

        Parameters
        ----------
        worksheet : openpyxl.Writer
            writer object for current sheet
        """
        for cells in worksheet.rows:
            for cell in cells:
                cell.font = Font(name=DEFAULT_FONT.name)


    def colour_hyperlinks(self, worksheet) -> None:
        """
        Set text colour to blue if text contains hyperlink

        Parameters
        ----------
        worksheet : openpyxl.Writer
            writer object for current sheet
        """
        for cells in worksheet.rows:
            for cell in cells:
                if 'HYPERLINK' in str(cell.value):
                    cell.font = Font(color='00007f', name=DEFAULT_FONT.name)


    def set_widths(self, worksheet, sheet_columns) -> None:
        """
        Set widths for variant sheets off common names to be more readable,
        calls get_closest_match() method to determine appropriate column
        width to use.

        Parameters
        ----------
        worksheet : openpyxl.Writer
            writer object for current sheet
        sheet_columns : list
            column names for sheet from DataFrame.columns
        """
        widths = {
            "chrom": 8,
            "pos": 12,
            "ref": 10,
            "alt": 10,
            "qual": 10,
            "af": 6,
            "dp": 10,
            'ac': 10,
            'af': 10,
            'an': 10,
            'dp': 10,
            'baseqranksum': 15,
            'clippingranksum': 16,
            "symbol": 12,
            "exon": 9,
            "variant class": 15,
            "consequence": 25,
            "hgvsc": 27,
            "hgvsp": 27,
            "gnomad": 13,
            "existing variation": 18,
            "clinvar": 10,
            "clinvar clndn": 18,
            "clinvar clinsig": 18,
            "cosmic": 15,
            "feature": 17
        }

        # generate list of 286 potential xlsx columns from A,B,C...JX,JY,JZ
        # allows for handling a lot of columns
        column_list = [
            f"{x}{y}" for x in [
                '', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J'
            ] for y in uppercase
        ]

        for idx, column in enumerate(sheet_columns):
            # loop over column names, select width by closest match from dict
            # and set width in sheet letter
            width = self.get_closest_match(column.lower(), widths)
            worksheet.column_dimensions[column_list[idx]].width = width


    def get_closest_match(self, column, widths) -> int:
        """
        Given a column name, find the closest match (if there is one) in the
        widths dict and return its width value to set. Using imprecise name
        matching as columns can differ between variant callers, annotation in
        VEP and renaming in vcf.rename_columns()

        Parameters
        ----------
        column : str
            name of column
        widths : dict
            dict of common column names and widths

        Returns
        -------
        width : int
            column width value to set
        """
        distances = {x: levenshtein.distance(column, x) for x in widths.keys()}
        closest_match = min(distances, key=distances.get)

        if distances[closest_match] <= 5:
            # close enough match to probably be correct
            width = widths[closest_match]
        else:
            # no close matches to name, use default width
            width = 13

        return width
