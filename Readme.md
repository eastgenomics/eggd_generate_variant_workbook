<!-- dx-header -->

# egg_generate_workbook (DNAnexus Platform App)

## What does this app do?

Generate an Excel workbook from VEP annotated vcf(s)

## What are typical use cases for this app?

This app may be executed as a standalone app.

## What data are required for this app to run?

**Packages**

* Packages from htslib DNAnexus asset:
  * bcftools
  * bcftools -split-vep plugin
* Python packages (specified in requirements.txt)

**File inputs (required)**:

- `--vcfs`: VEP annotated vcf(s)

**Other Inputs:**

`--assay` (string):  type of assay to generate report for, determines assay specific things in-app (currently only supports `dias`)

`--exclude` (string): Columns of VCF to exclude from output xlsx

`--include` (string): Columns of VCF to only include in output xlsx

`--reorder` (string): Order of fields from VCF in output xlsx columns, any not specified will be appended to end

`--rename` (string): = separated key value pairs of VCF fields to rename in output xlsx (e.g. `"gnomADg_AF=gnomAD_genomes_AF"`)

`--filter` (string): `bcftools filter` formatted string (see examples below)

`--keep` (bool): Determines if filtered rows from `--filter` are retained in a seperate 'filtered' tab (default: True)

`--add_name` (bool): Determines if to add samplename as first column in each sheet (default: False)

`--sheets` (list): Names to use for sheets where multiple VCFs are passed

`--output_prefix` (string): Prefix for naming output xlsx file. If not given for single file, the VCF prefix will be used. For 2+ input VCFs this must be specified.

`--merge` (bool): Determines if to merge multiple VCFs to one sheet (default: False)

`--summary` (string): If to include summary sheet, specify key of assay. Currently only supports `dias`.

`--flagstat_file` (file): Flagstat file to parse values from to display in summary.

`--panel_bed` (file): Panel bed to filter VCF with before parsing to xlsx file.

`--print_columns` (bool): Print total columns of all vcfs that will be output to the xlsx. Useful to identify what will be in the output to include/exclude.

`--print_header` (bool): Print header of first vcf, useful for inspecting all fields and types for filtering.

Example:

```bash
dx run app-eggd_vcf2xls/2.0.0 -ivcfs="file-G70BB1j45jFpjkPJ2ZB10f47" -ifilter="bcftools filter -e 'CSQ_gnomAD_AF>0.02''" --ikeep=true -irename_cols="gnomADg_AF=gnomAD_genomes_AF" -isummary="dias" -iassay="dias" -iexclude="MLEAC" -iexclude="MLEAF" -iexlcude="MQRankSum"
```

## What does this app output?

This app outputs an Excel workbook.

This is the source code for an app that runs on the DNAnexus Platform.
For more information about how to run or modify it, see
https://documentation.dnanexus.com/.

#### This app was made by EMEE GLH
