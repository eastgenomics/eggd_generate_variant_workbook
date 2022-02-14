<!-- dx-header -->

# vcf2xls_nirvana (DNAnexus Platform App)

## What does this app do?

Generate an xlsx report from a VEP annotated vcf

## What are typical use cases for this app?

This app may be executed as a standalone app.

## What data are required for this app to run?

**File inputs**:

- `--vcfs`: VEP annotated vcf(s)
- `--flagstat_file`: Flagstat file (optional)
- `--panel_bed`: flanked bed file used to filter vcf (optional)


**Other Inputs:**

`--assay` (string):  type of assay to generate report for, determines assay specific things in-app (currently only supports `dias`)

`--exclude_cols` (string): Columns of VCF to exclude from output xlsx

`--include_cols` (string): Columns of VCF to only include in output xlsx

`--reorder_cols` (string): Order of fields from VCF in output xlsx columns, any not specified will be appended to end

`--rename_cols` (string): = separated key value pairs of VCF fields to rename in output xlsx (e.g. `"gnomADg_AF=gnomAD_genomes_AF"`)

`--filter` (string): Columns on which to filter out variants. Format should be as `<column><operator><value>` e.g. `(gnomAD_AF<0.02)`. `Supported operands are >,`<,>,=,<=,== and !=.

`--keep` (bool): Determines if filtered rows from `--filter` are retained in a seperate 'filtered' tab

`--add_name` (bool): Determines if to add samplename as first column in each sheet

`--sheets` (bool): Names to use for sheets where multiple VCFs are passed

`--output_prefix` (string): Prefix for naming output xlsx file. If not given for single file, the VCF prefix will be used. For 2+ input VCFs this must be specified.

`--merge` (bool): Determines if to merge multiple VCFs to one sheet

`--summary` (string): If to include summary sheet, specify key of assay. Currently only supports `dias`.

`--flagstat_file` (file): Flagstat file to parse values from to display in summary.

`--panel_bed` (file): Panel bed to filter VCF with before parsing to xlsx file.

`--print-columns` (bool): Print total columns of all vcfs that will be output to the xlsx. Useful to identify what will be in the output to include/exclude.



Example:

```bash
dx run app-eggd_vcf2xls/2.0.0 -ivcfs="file-G70BB1j45jFpjkPJ2ZB10f47" -ifilter="gnomad_AF<0.02" -ifilter="DP>99" --ikeep=true -irename_cols="gnomADg_AF=gnomAD_genomes_AF" -isummary="dias" -iassay="dias" -iexclude="MLEAC" -iexclude="MLEAF" -iexlcude="MQRankSum"
```

## What does this app output?

This app outputs an xlsx report.

This is the source code for an app that runs on the DNAnexus Platform.
For more information about how to run or modify it, see
https://documentation.dnanexus.com/.

#### This app was made by EMEE GLH
