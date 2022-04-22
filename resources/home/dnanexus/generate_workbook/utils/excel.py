from pathlib import Path
from string import ascii_uppercase as uppercase
import sys
from timeit import default_timer as timer

import Levenshtein as levenshtein
import numpy as np
from openpyxl import load_workbook
from openpyxl.styles import Alignment, Border, DEFAULT_FONT, Font, Side
from openpyxl.styles.fills import PatternFill
import pandas as pd


# openpyxl style settings
THIN = Side(border_style="thin", color="000000")
MEDIUM = Side(border_style="medium", color="000001")
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
        if self.args.acmg:
            self.write_reporting_template()
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
        self.summary.cell(34, 1).value = "Total records:"

        # get sample name from vcf, should only be one but handle everything
        # list-wise just in case
        sample = [
            Path(x).name.replace('.vcf', '').replace('.gz', '')
            for x in self.args.vcfs
        ]
        sample = [x.split('_')[0] if '_' in x else x for x in sample]
        sample = str(sample).strip('[]').strip("'")

        # write summary values
        self.summary.cell(1, 2).value = sample
        self.summary.cell(1, 6).value = self.args.clinical_indication
        self.summary.cell(2, 6).value = self.args.panel

        # write total rows in each sheet
        count = 34

        for sheet, vcf in zip(self.args.sheets, self.vcfs):
            self.summary.cell(count, 2).value = sheet
            self.summary.cell(count, 3).value = len(vcf.index)
            self.summary[f"A{count}"].font = Font(
                bold=True, name=DEFAULT_FONT.name
            )
            count += 1

        count += 5

        # write genome reference(s) parsed from vcf header
        if self.refs:
            self.summary.cell(count, 1).value = "Reference:"
            self.summary[f"A{count}"].font = Font(
                bold=True, name=DEFAULT_FONT.name
            )
            for ref in list(set(self.refs)):
                self.summary.cell(count, 2).value = ref
                count += 1

        # write args passed to script to generate report
        self.summary.cell(count, 1).value = "Filters applied:"
        self.summary[f"A{count}"].font = Font(bold=True, name=DEFAULT_FONT.name)
        if self.args.filter:
            self.summary.cell(count, 2).value = self.args.filter
        else:
            self.summary.cell(count, 2).value = "None"

        count += 2

        self.summary.cell(count, 1).value = "Workflow:"
        self.summary.cell(count + 1, 1).value = "Workflow ID:"
        self.summary.cell(count + 2, 1).value = "Report Job ID:"
        self.summary[f"A{count}"].font = Font(bold=True, name=DEFAULT_FONT.name)
        self.summary[f"A{count + 1}"].font = Font(bold=True, name=DEFAULT_FONT.name)
        self.summary[f"A{count + 2}"].font = Font(bold=True, name=DEFAULT_FONT.name)

        self.summary.cell(count, 2).value = self.args.workflow[0]
        self.summary.cell(count + 1, 2).value = self.args.workflow[1]
        self.summary.cell(count + 2, 2).value = self.args.job_id


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
            start_row=1, end_row=1, start_column=2, end_column=4)
        self.summary.merge_cells(
            start_row=1, end_row=1, start_column=6, end_column=11)
        self.summary.merge_cells(
            start_row=2, end_row=2, start_column=6, end_column=11)
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
            "A1", "A34", "A35", "A36","A38", "B1",
            "B9", "B16", "B21", "B22", "B28", "B34", "B35", "B36", "B37",
            "C16", "C22", "D16", "D22", "D28", "E1", "E2", "E22",
            "F16", "F22", "G16", "G22", "H16", "H22", "I16"
        ]
        for cell in title_cells:
            self.summary[cell].font = Font(bold=True, name=DEFAULT_FONT.name)

        # set column widths for readability
        self.summary.column_dimensions['A'].width = 18
        self.summary.column_dimensions['B'].width = 13
        self.summary.column_dimensions['C'].width = 13
        self.summary.column_dimensions['D'].width = 13
        self.summary.column_dimensions['E'].width = 22
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


    def write_reporting_template(self) -> None:
        """
        Writes sheet to Excel file with formatting for reporting against
        ACMG criteria
        """
        report = self.workbook.create_sheet('report')

        titles = {
            "Gene": [2, 2],
            "HGVSc": [2, 3],
            "HGVSp": [2, 4],
            "Evidence": [4, 3],
            "Pathogenic": [4, 7],
            "Yes / No": [4, 8],
            "Benign": [4, 9],
            "Yes / No": [4, 10],
            "Associated disease": [5, 2],
            "Known inheritance": [6, 2],
            "Prevalence": [7, 2],
            "Estimated allele frequency": [8, 2],
            "Null variant where LOF function of disease": [9, 2],
            "Same AA change as pathogenic change,regardless\nof nucleotide": [10, 2],
            "De novo inheritance or inheritance confirmed / observed in\nhealthy adult": [11, 2],
            "In vivo / in vitro functional studies": [12, 2],
            "Prevalence in affected > controls": [13, 2],
            "In mutational hotspot, without benign variation": [14, 2],
            "Freq in controls eg ExAC, low/absent or >5%": [15, 2],
            "Confirmation of in trans/in cis with pathogenic variant": [16, 2],
            "In frame protein length change, non repeating vs. repeating": [17, 2],
            "Same AA as a different pathogenic change": [18, 2],
            "Cosegregation with disease in family, not in unnaffected": [18, 2],
            ("Missense where low rate of benign missense and common\nmechanism"
                "(Z score >3), or missense where LOF common\nmechanism"): [19, 2],
            "Multiple lines of computational evidence (Cant use with PS3)": [20, 2],
            ("Phenotype/FH specific for disease of single etiology, or\n"
                "alternative genetic cause of disease detected"): [21, 2],
            "Reputable source reports but evidence not available": [22, 2],
            "Synonymous change, no affect on splicing, not conserved": [23, 2],
            "ACMG Classification": [26, 2],
        }

        for key, val in titles.items():
            report.cell(val[0], val[1]).value = key
            report.cell(val[0], val[1]).font = Font(
                bold=True, name=DEFAULT_FONT.name
            )

        classifications = {
            "Extra": [(5, 7), (6, 7), (7, 7), (8, 7), (5, 9), (6, 9), (7, 9)],
            "PVS1": [(9, 7)],
            "PS1": [(10, 7)],
            "PS2": [(11, 7)],
            "PS3": [(12, 7)],
            "PS4": [(13, 7)],
            "PM1": [(14, 7)],
            "PM2": [(15, 7)],
            "PM3": [(16, 7)],
            "PM4": [(17, 7)],
            "PM5": [(18, 7)],
            "PM6": [(19, 7)],
            "PP1": [(20, 7)],
            "PP2": [(21, 7)],
            "PP3": [(22, 7)],
            "PP4": [(23, 7)],
            "PP5": [(24, 7)],
            "BS1": [(8, 9)],
            "BS2": [(11, 9)],
            "BS3": [(12, 9)],
            "BA1": [(15, 9)],
            "BP2": [(16, 9)],
            "BP3": [(17, 9)],
            "BS4": [(20, 9)],
            "BP1": [(21, 9)],
            "BP4": [(22, 9)],
            "BP5": [(23, 9)],
            "BP6": [(24, 9)],
            "BP7": [(25, 9)]
        }

        for key, values in classifications.items():
            for val in values:
                report.cell(val[0], val[1]).value = key

        # nice formatting of title text
        for cell in report['B']:
            breaks = str(cell.value).count("\n") + 1
            report.row_dimensions[cell.row].height = 20 * breaks
            report[f"B{cell.row}"].alignment = Alignment(
                wrapText=True, vertical="center"
            )

        # merge evidence cells
        for row in range(4, 26):
            report.merge_cells(
                start_row=row, end_row=row, start_column=3, end_column=6)

        # set appropriate widths
        report.column_dimensions['B'].width = 60
        report.column_dimensions['C'].width = 35
        report.column_dimensions['D'].width = 35
        report.column_dimensions['E'].width = 5
        report.column_dimensions['F'].width = 5
        report.column_dimensions['G'].width = 12


        # do some colouring
        colour_cells = {
            'FAC090': ['B2', 'B3', 'C2', 'C3', 'D2', 'D3'],
            '8EB4E3': ['B4', 'B5', 'B6', 'B7', 'B8', 'C4', 'G4', 'H4', 'I4', 'J4'],
            'FFFF99': ['B26', 'G5', 'G6', 'G7', 'G8', 'I5', 'I6', 'I7', 'I8'],
            'E46C0A': ['G9', 'G10', 'G11', 'G12', 'G13'],
            'FFC000': ['G14', 'G15', 'G16', 'G17', 'G18', 'G19'],
            'FFFF00': ['G20', 'G21', 'G22', 'G23', 'G24'],
            '00B0F0': ['I8', 'I11', 'I12', 'I20'],
            '92D050': ['I16', 'I17', 'I21', 'I22', 'I23', 'I24', 'I25'],
            'FF0000': ['G9'],
            'D9D9D9': [
                'I9', 'I10', 'J9', 'J10', 'I13', 'I14', 'J13', 'J14',
                'I18', 'I19', 'J18', 'J19', 'G25', 'H25'
            ]

        }
        for colour, cells in colour_cells.items():
            for cell in cells:
                report[cell].fill = PatternFill(
                    patternType="solid", start_color=colour
                )

        # add some borders
        row_ranges = {
            'horizontal': [
                'B3:D3', 'B4:J4', 'B5:J5',
                'B6:J6', 'B7:J7', 'B8:J8', 'B9:J9', 'B10:J10', 'B11:J11',
                'B12:J12', 'B13:J13', 'B14:J14', 'B15:J15', 'B16:J16',
                'B17:J17', 'B18:J18', 'B19:J19', 'B20:J20', 'B21:J21',
                'B22:J22', 'B23:J23', 'B24:J24', 'B25:J25'
            ],
            'horizontal_thick': [
                'B2:D2', 'B4:J4', 'B26:J26', 'B27:J27'
            ],
            'vertical': [
                'E2:E3', 'G4:G26', 'H4:H26', 'I4:I26', 'J4:J26'
            ],
            'vertical_thick': [
                'B2:B26', 'C2:C26', 'K4:K26', 'E2:E3'
            ]
        }

        for side, values in row_ranges.items():
            for row in values:
                for cells in report[row]:
                    for cell in cells:
                        # border style is immutable => copy current and modify
                        cell_border = cell.border.copy()
                        if side == 'horizontal':
                            cell_border.top = THIN
                        if side == 'horizontal_thick':
                            cell_border.top = MEDIUM
                        if side == 'vertical':
                            cell_border.left = THIN
                        if side == 'vertical_thick':
                            cell_border.left = MEDIUM
                        cell.border = cell_border


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

                if self.args.add_comment_column:
                    vcf = vcf.insert(loc=-1, column='comment')

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

                # set Excel types for numeric cells to suppress Excel warnings
                self.set_types(curr_worksheet)
                self.workbook.save(self.args.output)


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

        # force nan values to be None strings for consistency on comparing
        vcf = vcf.replace(r'^\s*$', np.nan, regex=True)
        vcf.fillna('None', inplace=True)
        written_sheet = written_sheet.replace(r'^\s*$', np.nan, regex=True)
        written_sheet.fillna('None', inplace=True)

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


    def set_types(self, worksheet) -> None:
        """
        Iterate over all worksheet cells and test if cell value can be numeric,
        if so sets the type to numeric to suppress 'Number stored as text'
        warning in Excel

        Parameters
        ----------
        worksheet : openpyxl.Writer
            writer object for current sheet
        """
        for cells in worksheet.rows:
            for cell in cells:
                if self.is_numeric(cell.value):
                    cell.data_type = 'n'


    def is_numeric(self, value) -> bool:
        """
        Returns true if given value is in some form numeric

        Parameters
        ----------
        value : str
            string to check if can be cast to numeric type

        Returns
        -------
        bool
            True if value can be numeric
        """
        return str(value).lstrip('-').replace('.', '').replace(
            'e-', '', 1).replace('e', '').isdigit()


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
            "gnomad_af": 20,
            "gnomad_exomes_af": 20,
            "gnomad_genomes_af": 20,
            "hgmd": 13,
            "hgmd_phen": 15,
            "spliceai_ds": 18,
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
