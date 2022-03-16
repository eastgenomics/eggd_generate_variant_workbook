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

`--exclude_columns` (`string`): Columns of VCF to exclude from output workbook.

`--include_columns` (`string`): Columns of VCF to only include in output workbook.

`--reorder_columns` (`string`): Order of fields from VCF in output xlsx columns, any not specified will be appended to end.

`--rename_columns` (`string`): `=` separated key value pairs of VCF fields to rename in output xlsx (e.g. `"gnomADg_AF=gnomAD_genomes_AF"`).

`--filter` (`string`): `bcftools filter` formatted string (see examples below).

`--types` (`string`): `=` separated key value pairs of `{field}={type}` to overwrite in VCF header (i.e `CSQ_gnomADg_AF=Float`). This is required where the field is specified for filtering, but the type is wrongly set in the header. Field types for a given vcf may be inspected with `--print_header`.

`--keep_tmp` (`bool`): Determines if to upload the intermediate bcftools split and filtered vcfs

`--keep_filtered` (`bool`): Determines if filtered rows from `--filter` are retained in a seperate 'filtered' tab (default: `True`).

`--add_samplename` (`bool`): Determines if to add sample name as first column in each sheet (default: `False`).

`--sheet_names` (`list`): Names to use for workbook sheets, these MUST be the same number as the number of vcfs passed and in the same order. If not given, and if there is 1 vcf passed the sheet will be named `variants`, else if multiple vcfs are passed the name prefix of the vcf will be used.

`--output_prefix` (`string`): Prefix for naming output xlsx file. If not given for single file, the VCF prefix will be used. For 2+ input VCFs this must be specified.

`--merge_vcfs` (`bool`): Determines if to merge multiple VCFs to one sheet (default: `False`).

`--summary` (`string`): If to include summary sheet, specify key of assay. Currently only supports `dias`.

`--panel_bed` (`file`): Panel bed to filter VCF with before parsing to xlsx file.

`--print_columns` (`bool`): Print column names of all vcfs that will be output to the xlsx. Useful to identify what will be in the output to include/exclude.

`--print_header` (`bool`): Print header of first vcf, useful for inspecting all fields and types for filtering.


**Example**:

```bash
dx run app-eggd_vcf2xls/2.0.0 \
  -ivcfs="file-G70BB1j45jFpjkPJ2ZB10f47" \
  -ifilter="bcftools filter -e 'CSQ_gnomAD_AF>0.02 | CSQ_gnomADg_AF>0.02'" \
  -itypes="CSQ_gnomADg_AF=Float" \
  -irename_cols="gnomADg_AF=gnomAD_genomes_AF" \
  -isummary="dias" \
  -iexclude="MLEAF MQRankSum MLEAC" \
```

The above will do the following:

- filter **out** variants with gnomAD exomes or genomes AF above 2%
- set the type for gnomAD genomes AF in the vcf header to `Float` to allow filtering
- rename the gnomAD genomes AF column
- add the summary sheet for Dias
- exlcude the `MLEAF`, `MGRankSum` and `MLEAC` columns from the output
- keep the filtered out variants in a separate  sheet (default behaviour)


**Notes on inputs**

- `-summary` currently only accepts `dias` as an option, which adds the summary sheet specific for samples processed through the Dias pipeline
- filtering with `-filter` accepts any valid `bcftools` command where the vcf would be passed at the end (i.e. `bcftools filter -e 'CSQ_gnomAD_AF>0.02' sample.vcf`). Commands that are chained with pipes where the vcf is passed to the first command will NOT work (i.e. `bcftools filter -e 'CSQ_gnomAD_AF>0.02' sample.vcf | bcftoools ...`). See the filtering section below for examples.
- `-include` and `-exclude` are mutually exclusive and can not be passed together
- any arguments that are strings and passing multiple should be one space seperated string (i.e. `-iexclude="MLEAF MQRankSum MLEAC"`)


**Filtering**

Filtering of the input vcf(s) is achieved by passing `bcftools` commands to the `-filter` argument. These should be passed without specifying the vcf or `-o` output argument, as these are set in the app. These require passing the exact name of the field as will be in the vcf, any fields coming from the VEP annotation will be prefixed with `CSQ_` (i.e. the `gnomAD_AF` will be `CSQ_gnomAD_AF`). The names and types set can be viewed by passing `-view_header`.

If the type is wrongly set in the vcf header (i.e. AF values set `String`), these will not be valid for filtering. Therefore, if you are specifying these fields for filtering, the `-types` argument must also be passed to specify the correct field type for the given column (i.e. `-itypes="CSQ_gnomADg_AF=Float"`).

Furthermore, the filtering allows for using the logic described in the [bcftools documentation](bcftools). No checking of what is being filtered is done (other than if bcftools returns a non-zero exit code, in which case an AssertionError will be raised). Therefore, it is up to the user to ensure the bcftools command being passed filters as intended.

Please be aware of the difference between using `-i / --include` and `-e / --exclude`, as they may not filter the same (i.e. including under 2% (`-i 'gnomAD_AF<0.02'`) != excluding over 2% (`-e gnomAD_AF>0.02`)). In this case, missing values (represented in the vcf with a `'.'`) will be filtered **out** with `--include` as they will not meet the filter. If they should be retained, `--exclude` should be used.

**Examples of filters**
```
# filtering gnomAD at 2%
-ifilter="bcftools filter -e 'CSQ_gnomAD_AF>0.02'"

# filtering gnomAD exomes and genomes at 2%
-ifilter="bcftools filter -e 'CSQ_gnomAD_AF>0.02 | CSQ_gnomADg_AF>0.02'"

# filtering out synonymous consequence
-ifilter="bcftools filter -e 'CSQ_Consequence==\"synonymous_variant\"'"
```

Useful resources:

- [bcftools expressions](bcftools-expressions)



## What does this app output?

This app outputs an Excel workbook and (optionally) the intermediate split and filtered vcfs.

This is the source code for an app that runs on the DNAnexus Platform.
For more information about how to run or modify it, see
https://documentation.dnanexus.com/.

#### This app was made by EMEE GLH

[bcftools]:https://samtools.github.io/bcftools/bcftools.html#filter
[bcftools-expressions]:https://samtools.github.io/bcftools/bcftools.html#expressions

