# eggd_generate_variant_workbook: Technical Walkthrough

*2026-02-20T13:29:44Z by Showboat 0.6.0*
<!-- showboat-id: 8b86b075-51e0-47ea-ae2c-9f282ba1aee1 -->

## Overview

`eggd_generate_variant_workbook` (v2.11.1) is a DNAnexus platform app that converts VEP-annotated VCF files into formatted Excel workbooks for clinical genomics reporting. It supports multiple diagnostic pipelines — DIAS, HELIOS, and URANUS — and produces workbooks with filtered variants, conditional formatting, hyperlinks to external databases (gnomAD, ClinVar, COSMIC, HGMD, MasterMind, DECIPHER, OncoKB, cBioPortal), and optional summary sheets tailored to each assay type.

**Primary use case**: Transform raw variant calls from NGS pipelines into professionally formatted spreadsheets ready for clinical or research reporting.

## Repository Structure

```bash
find . -not -path '*/__pycache__/*' -not -name '__pycache__' -not -path '*/.git*' -not -name '*.pyc' | sed 's|^\./||' | grep -v '^\.$' | sort | head -80
```

```output
dxapp.json
Readme.md
requirements.txt
resources
resources/home
resources/home/dnanexus
resources/home/dnanexus/generate_workbook
resources/home/dnanexus/generate_workbook/generate_workbook.py
resources/home/dnanexus/generate_workbook/tests
resources/home/dnanexus/generate_workbook/tests/__init__.py
resources/home/dnanexus/generate_workbook/tests/pytest.ini
resources/home/dnanexus/generate_workbook/tests/test_columns.py
resources/home/dnanexus/generate_workbook/tests/test_data
resources/home/dnanexus/generate_workbook/tests/test_data/column_methods_test.vcf.gz
resources/home/dnanexus/generate_workbook/tests/test_data/HD753-unittest_annotated.split.vcf
resources/home/dnanexus/generate_workbook/tests/test_data/NA12878_unittest.split.vcf
resources/home/dnanexus/generate_workbook/tests/test_data/NA12878_unittest.vcf
resources/home/dnanexus/generate_workbook/tests/test_data/oncospan_annotated.vcf.gz
resources/home/dnanexus/generate_workbook/tests/test_data/unannotated.vcf.gz
resources/home/dnanexus/generate_workbook/tests/test_excel.py
resources/home/dnanexus/generate_workbook/tests/test_filters.py
resources/home/dnanexus/generate_workbook/tests/test_generate_workbook.py
resources/home/dnanexus/generate_workbook/tests/test_utils.py
resources/home/dnanexus/generate_workbook/tests/test_vcf.py
resources/home/dnanexus/generate_workbook/utils
resources/home/dnanexus/generate_workbook/utils/columns.py
resources/home/dnanexus/generate_workbook/utils/excel.py
resources/home/dnanexus/generate_workbook/utils/filters.py
resources/home/dnanexus/generate_workbook/utils/__init__.py
resources/home/dnanexus/generate_workbook/utils/utils.py
resources/home/dnanexus/generate_workbook/utils/vcf.py
resources/home/dnanexus/packages
resources/home/dnanexus/packages/colour-0.1.5-py2.py3-none-any.whl
resources/home/dnanexus/packages/et_xmlfile-1.1.0-py3-none-any.whl
resources/home/dnanexus/packages/filetype-1.1.0-py2.py3-none-any.whl
resources/home/dnanexus/packages/lxml-4.9.3-cp38-cp38-manylinux_2_17_x86_64.manylinux2014_x86_64.manylinux_2_24_x86_64.whl
resources/home/dnanexus/packages/numpy-1.22.2-cp38-cp38-manylinux_2_17_x86_64.manylinux2014_x86_64.whl
resources/home/dnanexus/packages/openpyxl-3.1.2-py2.py3-none-any.whl
resources/home/dnanexus/packages/pandas-1.4.1-cp38-cp38-manylinux_2_17_x86_64.manylinux2014_x86_64.whl
resources/home/dnanexus/packages/Pillow-9.2.0-cp38-cp38-manylinux_2_17_x86_64.manylinux2014_x86_64.whl
resources/home/dnanexus/packages/python_Levenshtein_wheels-0.13.2-cp38-cp38-manylinux2010_x86_64.whl
resources/home/dnanexus/packages/pytz-2021.3-py2.py3-none-any.whl
resources/usr
resources/usr/local
resources/usr/local/bin
resources/usr/local/bin/jq-linux64
resources/usr/local/bin/mark-section
resources/usr/local/bin/mark-success
src
src/code.sh
technical_walkthrough.md
```

The repo follows the DNAnexus applet layout convention:

- `src/code.sh` — DNAnexus bash entry point; downloads inputs, installs packages, and invokes the Python script
- `resources/home/dnanexus/generate_workbook/` — all Python source code, deployed to the worker's home directory at runtime
  - `generate_workbook.py` — CLI argument parsing and top-level orchestration
  - `utils/vcf.py` — VCF reading, bcftools integration, and all data transformations
  - `utils/excel.py` — Excel workbook construction and formatting
  - `utils/filters.py` — bcftools-based variant filtering
  - `utils/columns.py` — INFO / FORMAT / VEP column splitting
  - `utils/utils.py` — hyperlink builders and miscellaneous helpers
- `resources/home/dnanexus/packages/` — vendored Python wheels for offline installation on DNAnexus workers
- `dxapp.json` — DNAnexus app manifest (inputs, outputs, instance type, etc.)
- `tests/` — pytest suite co-located with source

## Command-Line Interface

The Python entry point accepts arguments directly. On DNAnexus, `code.sh` constructs and passes these arguments automatically. Locally, invoke as:

    python3 generate_workbook.py --vcfs sample.vcf.gz [OPTIONS]

Key argument groups:

| Group | Arguments | Purpose |
|-------|-----------|---------|
| **Input** | `--vcfs`, `--additional_files`, `--images`, `--m_codes` | Data inputs |
| **Column selection** | `--include`, `--exclude`, `--reorder`, `--rename` | Control which columns appear and in what order |
| **Annotation columns** | `--add_comment_column`, `--add_classification_column`, `--add_allele_origin_column`, `--add_report_text_column`, `--split_hgvs`, `--add_raw_change` | Append empty or derived columns to each sheet |
| **External links** | `--additional_columns decipher oncokb cbioportal pecan` | Add clickable hyperlink columns to external databases |
| **Filtering** | `--filter`, `--types`, `--keep` | bcftools filter expression; field-type overrides; keep excluded variants in a separate tab |
| **Output** | `--output`, `--out_dir`, `--sheets`, `--merge` | Filename prefix, destination dir, sheet names, merge VCFs |
| **Formatting** | `--colour`, `--af_format`, `--freeze_column`, `--sort_by`, `--add_auto_filter` | Cell colouring, AF as decimal/percent, freeze rows/cols, sort order |
| **Summary sheet** | `--summary {dias,helios,uranus}`, `--panel`, `--clinical_indication`, `--acmg`, `--lock_sheet` | Prepend a pipeline-specific summary tab, add ACMG templates |
| **Diagnostics** | `--print_columns`, `--print_header`, `--keep_tmp` | Inspect available columns/header without writing output |

## Processing Pipeline

The end-to-end data flow from raw VCF to formatted Excel workbook:

    Input VCF(s) (VEP-annotated)
             │
             ▼
      1. Argument validation (generate_workbook.py → arguments class)
             │
             ▼
      2. VCF pre-processing (utils/vcf.py → vcf.process)
         ├─ Check VEP annotation present (check_vep_vcf)
         ├─ bcftools +split-vep  → one row per transcript
         ├─ bcftools filter      → tag excluded variants (utils/filters.py)
         ├─ Expand INFO / FORMAT columns (utils/columns.py)
         ├─ Derived columns: split_hgvs, add_raw_change, make_report_text, join_columns
         ├─ Hyperlink columns added (add_hyperlinks, add_additional_columns)
         ├─ Drop / reorder / rename columns
         └─ Sort variants
             │
             ▼
      3. Additional file ingestion (TSV/CSV/VCF; TSO500 special handling)
             │
             ▼
      4. Excel generation (utils/excel.py → excel.generate)
         ├─ Summary sheet  (dias / helios / uranus specific layout)
         ├─ ACMG reporting template sheets
         ├─ Variant sheet(s)  — one per VCF (or merged)
         ├─ Additional file sheet(s)
         ├─ Image sheet(s)
         └─ Formatting pass
               ├─ Column widths, alignment, borders
               ├─ Conditional cell colouring
               ├─ Hyperlink styling
               ├─ Data-validation dropdowns
               ├─ Sheet protection (lock_sheet)
               └─ Freeze panes
             │
             ▼
      Output: sample.xlsx  +  DNAnexus variant-count metadata

## Step 1 — VCF Input and VEP Annotation Detection

The tool requires VEP-annotated VCFs. `check_vep_vcf` confirms the presence of a `##VEP` header line and a `CSQ` INFO field before any processing begins. It also detects whether the VCF has already been split by bcftools (checking for individual CSQ sub-field headers) to avoid double-splitting.

```bash
bcftools view -h resources/home/dnanexus/generate_workbook/tests/test_data/NA12878_unittest.vcf | grep -E '^##(VEP|INFO=<ID=CSQ)' | head -5
```

```output
##VEP="v105" time="2022-05-10 17:57:27" cache="/opt/vep/.vep/homo_sapiens_refseq/105_GRCh37" ensembl=105.f357e33 ensembl-variation=105.ac8178e ensembl-funcgen=105.660df8f ensembl-io=105.2a0a40c 1000genomes="phase3" COSMIC="92" ClinVar="202012" HGMD-PUBLIC="20204" assembly="GRCh37.p13" dbSNP="154" gencode="GENCODE 19" genebuild="2011-04" gnomAD="r2.1" polyphen="2.2.2" refseq="2020-10-26 17:03:42 - GCF_000001405.25_GRCh37.p13_genomic.gff" regbuild="1.0" sift="sift5.2.2"
##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: Allele|SYMBOL|HGNC_ID|VARIANT_CLASS|Consequence|IMPACT|EXON|INTRON|Feature|HGVSc|HGVSp|HGVS_OFFSET|Existing_variation|STRAND|ClinVar|ClinVar_CLNSIG|ClinVar_CLNSIGCONF|ClinVar_CLNDN|gnomADg_AC|gnomADg_AN|gnomADg_AF|gnomADe_AC|gnomADe_AN|gnomADe_AF|gnomADe_Hom|TWE_AF|TWE_AC_Hom|TWE_AC_Het|TWE_AN|HGMD|HGMD_PHEN|HGMD_CLASS|HGMD_RANKSCORE|SpliceAI_pred_DS_AG|SpliceAI_pred_DS_AL|SpliceAI_pred_DS_DG|SpliceAI_pred_DS_DL|SpliceAI_pred_DP_AG|SpliceAI_pred_DP_AL|SpliceAI_pred_DP_DG|SpliceAI_pred_DP_DL|REVEL|Mastermind_MMID3|CADD_PHRED">
```

The `##VEP` line records the VEP version and databases used at annotation time. The `##INFO=<ID=CSQ,...>` line defines the pipe-delimited sub-fields packed into each variant's `CSQ` field. In the example above, 44 annotation fields are present — from transcript consequences (`SYMBOL`, `HGVSc`, `HGVSp`) through population frequencies (`gnomADe_AF`, `TWE_AF`) and pathogenicity predictors (`CADD_PHRED`, `REVEL`, `SpliceAI_pred_*`).

## Step 2 — CSQ Field Splitting with bcftools +split-vep

VEP writes one CSQ entry per affected transcript for each variant, all packed into a single comma-separated string. The tool calls `bcftools +split-vep` to explode this into one row per transcript, making every annotation a proper VCF INFO field:

```bash
bcftools +split-vep resources/home/dnanexus/generate_workbook/tests/test_data/NA12878_unittest.vcf -c - -d 2>/dev/null | grep -v '^#' | cut -f1-8 | head -3
```

```output
1	1271940	rs149812757	C	T	748.77	.	AC=1;AF=0.5;AN=2;BaseQRankSum=-2.623;ClippingRankSum=0;DB;DP=55;ExcessHet=3.0103;FS=0;MLEAC=1;MLEAF=0.5;MQ=60;MQRankSum=0;QD=13.87;ReadPosRankSum=1.764;SOR=0.788;CSQ=T|DVL1||SNV|intron_variant|MODIFIER||14/14|NM_004421.3|NM_004421.3:c.1640-45G>A|||rs149812757|-1|1188053|Likely_benign||not_provided|334|30884|0.0108147|1674|129544|0.0129223|18|0.016|0|16|1000|||||0.01|0.00|0.00|0.00|-15|-1|-4|38|||2.674;Allele=T;SYMBOL=DVL1;HGNC_ID=.;VARIANT_CLASS=SNV;Consequence=intron_variant;IMPACT=MODIFIER;EXON=.;INTRON=14/14;Feature=NM_004421.3;HGVSc=NM_004421.3:c.1640-45G>A;HGVSp=.;HGVS_OFFSET=.;Existing_variation=rs149812757;STRAND=-1;ClinVar=1188053;ClinVar_CLNSIG=Likely_benign;ClinVar_CLNSIGCONF=.;ClinVar_CLNDN=not_provided;gnomADg_AC=334;gnomADg_AN=30884;gnomADg_AF=0.0108147;gnomADe_AC=1674;gnomADe_AN=129544;gnomADe_AF=0.0129223;gnomADe_Hom=18;TWE_AF=0.016;TWE_AC_Hom=0;TWE_AC_Het=16;TWE_AN=1000;HGMD=.;HGMD_PHEN=.;HGMD_CLASS=.;HGMD_RANKSCORE=.;SpliceAI_pred_DS_AG=0.01;SpliceAI_pred_DS_AL=0;SpliceAI_pred_DS_DG=0;SpliceAI_pred_DS_DL=0;SpliceAI_pred_DP_AG=-15;SpliceAI_pred_DP_AL=-1;SpliceAI_pred_DP_DG=-4;SpliceAI_pred_DP_DL=38;REVEL=.;Mastermind_MMID3=.;CADD_PHRED=2.674
1	1273116	rs307371	A	G	165.8	.	AC=2;AF=1;AN=2;DB;DP=7;ExcessHet=3.0103;FS=0;MLEAC=2;MLEAF=1;MQ=60;QD=23.69;SOR=2.584;CSQ=G|DVL1||SNV|intron_variant|MODIFIER||14/14|NM_004421.3|NM_004421.3:c.1639+241T>C|||rs307371|-1|1277157|Benign||not_provided|25171|30908|0.814385|107666|124494|0.864829|46999|0.908|844|64|1000|||||0.00|0.00|0.00|0.00|28|29|-35|29|||0.824;Allele=G;SYMBOL=DVL1;HGNC_ID=.;VARIANT_CLASS=SNV;Consequence=intron_variant;IMPACT=MODIFIER;EXON=.;INTRON=14/14;Feature=NM_004421.3;HGVSc=NM_004421.3:c.1639+241T>C;HGVSp=.;HGVS_OFFSET=.;Existing_variation=rs307371;STRAND=-1;ClinVar=1277157;ClinVar_CLNSIG=Benign;ClinVar_CLNSIGCONF=.;ClinVar_CLNDN=not_provided;gnomADg_AC=25171;gnomADg_AN=30908;gnomADg_AF=0.814385;gnomADe_AC=107666;gnomADe_AN=124494;gnomADe_AF=0.864829;gnomADe_Hom=46999;TWE_AF=0.908;TWE_AC_Hom=844;TWE_AC_Het=64;TWE_AN=1000;HGMD=.;HGMD_PHEN=.;HGMD_CLASS=.;HGMD_RANKSCORE=.;SpliceAI_pred_DS_AG=0;SpliceAI_pred_DS_AL=0;SpliceAI_pred_DS_DG=0;SpliceAI_pred_DS_DL=0;SpliceAI_pred_DP_AG=28;SpliceAI_pred_DP_AL=29;SpliceAI_pred_DP_DG=-35;SpliceAI_pred_DP_DL=29;REVEL=.;Mastermind_MMID3=.;CADD_PHRED=0.824
1	1273278	rs307370	A	G	1138.77	.	AC=2;AF=1;AN=2;DB;DP=33;ExcessHet=3.0103;FS=0;MLEAC=2;MLEAF=1;MQ=60;QD=34.51;SOR=4.112;CSQ=G|DVL1||SNV|intron_variant|MODIFIER||14/14|NM_004421.3|NM_004421.3:c.1639+79T>C|||rs307370|-1|1243362|Benign||not_provided|22727|29146|0.779764|94197|132642|0.71016|34586|0.812|662|150|1000|||||0.01|0.00|0.00|0.00|0|-12|0|2||DVL1:G547int|7.352;Allele=G;SYMBOL=DVL1;HGNC_ID=.;VARIANT_CLASS=SNV;Consequence=intron_variant;IMPACT=MODIFIER;EXON=.;INTRON=14/14;Feature=NM_004421.3;HGVSc=NM_004421.3:c.1639+79T>C;HGVSp=.;HGVS_OFFSET=.;Existing_variation=rs307370;STRAND=-1;ClinVar=1243362;ClinVar_CLNSIG=Benign;ClinVar_CLNSIGCONF=.;ClinVar_CLNDN=not_provided;gnomADg_AC=22727;gnomADg_AN=29146;gnomADg_AF=0.779764;gnomADe_AC=94197;gnomADe_AN=132642;gnomADe_AF=0.71016;gnomADe_Hom=34586;TWE_AF=0.812;TWE_AC_Hom=662;TWE_AC_Het=150;TWE_AN=1000;HGMD=.;HGMD_PHEN=.;HGMD_CLASS=.;HGMD_RANKSCORE=.;SpliceAI_pred_DS_AG=0.01;SpliceAI_pred_DS_AL=0;SpliceAI_pred_DS_DG=0;SpliceAI_pred_DS_DL=0;SpliceAI_pred_DP_AG=0;SpliceAI_pred_DP_AL=-12;SpliceAI_pred_DP_DG=0;SpliceAI_pred_DP_DL=2;REVEL=.;Mastermind_MMID3=DVL1:G547int;CADD_PHRED=7.352
```

After splitting, each row's INFO field contains individual key=value pairs (e.g., `SYMBOL=DVL1`, `Consequence=intron_variant`, `gnomADe_AF=0.0129`). The internal `splitColumns` class then reads these into a pandas DataFrame — one column per key — which forms the basis for the Excel sheet.

## Step 3 — Variant Filtering (utils/filters.py)

Filtering is driven by a bcftools filter expression passed via `--filter`. The approach uses *soft-filtering*: rather than discarding variants, failing records are tagged with `EXCLUDE` in the FILTER column. The tool then separates the DataFrame into passing and excluded subsets:

- Passing variants → primary variant sheet
- Excluded variants → separate `filtered` sheet (when `--keep` is set)

This preserves all data while making the clinical view clean. A verification step confirms that passing + excluded counts equal the original total, guarding against data loss.

**Example filter expression** (exclude common variants):

    --filter "-e 'CSQ_gnomADe_AF > 0.02'"

If a field is typed incorrectly in the VCF header (e.g., `AF` stored as `String` instead of `Float`), `--types CSQ_gnomADe_AF=Float` rewrites the header type before filtering.

## Step 4 — Column Transformations

### INFO and FORMAT expansion (utils/columns.py)

`splitColumns.info()` splits the semicolon-delimited INFO column into individual columns. `splitColumns.format_fields()` joins the FORMAT and SAMPLE columns (e.g., `GT:AD:DP` / `0/1:45,12:57`) into per-sample columns like `GT`, `AD`, `DP`, appending `(FMT)` to disambiguate any name clash with an INFO field of the same name.

### Derived columns

| Option | Output column | Logic |
|--------|--------------|-------|
| `--split_hgvs` | `DNA`, `Protein` | Strips transcript prefix from HGVSc/HGVSp (e.g., `NM_004421.3:c.1640-45G>A` → `c.1640-45G>A`) |
| `--add_raw_change` | `rawChange` | Concatenates `CHROM:g.POSREF>ALT` (e.g., `1:g.1271940C>T`) |
| `--add_report_text_column` | `report` | Constructs a summary sentence: gene, consequence, exon/intron, HGVS, rsID, and AF as a percentage |
| `--join_columns A:B:sep` | `A` | Appends column B to A with separator `sep`; or prepends/appends a literal string |

### Hyperlinks (utils/utils.py → buildHyperlink)

Fields are converted to clickable Excel hyperlinks based on the reference build detected from the VCF header:

| Database | Build | Field used |
|----------|-------|-----------|
| gnomAD | GRCh38 → v4, GRCh37 → v2 | CHROM, POS, REF, ALT |
| ClinVar | either | ClinVar accession |
| COSMIC | either | COSMIC ID |
| HGMD | either | HGMD accession |
| MasterMind | either | Mastermind MMID3; requires chr→NC mapping |
| DECIPHER | GRCh38 only | CHROM, POS, REF, ALT |
| OncoKB | either | Gene SYMBOL |
| cBioPortal | either | Gene SYMBOL |
| PeCan | either | Gene SYMBOL |

## Step 5 — Excel Generation (utils/excel.py)

`excel.generate()` orchestrates workbook construction using **openpyxl**. Sheets are written in this order:

1. **Summary sheet** (optional, pipeline-specific)
2. **ACMG reporting template sheets** (optional, `--acmg N` adds N copies)
3. **Variant sheets** — one per VCF, or a single merged sheet if `--merge` is set
4. **Additional file sheets** — TSV/CSV files appended as-is
5. **Image sheets** — each image on its own sheet

### Formatting applied to variant sheets

- **Column widths**: Set to match content length; long columns capped
- **Cell types**: Numeric-looking values stored as Excel numbers (not strings), enabling proper sorting/filtering
- **Allele frequency**: Optionally formatted as percentage (`--af_format percent`)
- **Conditional cell colouring** (`--colour`): Expressions like `VF:>=0.9:green` or `VF:<0.6&>=0.3:orange` paint cells based on value thresholds; supports compound AND (`&`) or OR (`|`) conditions
- **Hyperlink styling**: Hyperlinked cells coloured blue by default
- **Freeze panes**: Keeps header row (and optionally leading columns) always visible; default `A2`
- **Sheet protection** (`--lock_sheet`): Locks all cells except designated annotation columns (comment, classification, allele origin, etc.) — up to 500 rows × 200 columns unlocked for manual entry
- **Auto-filter** (`--add_auto_filter`): Adds Excel dropdown filters to the header row
- **Data-validation dropdowns**: Certain annotation columns (e.g., allele origin) use restricted pick-lists

## Pipeline-Specific Summary Sheets

Three `--summary` modes generate a pre-populated first sheet tailored to each assay:

### DIAS (diagnostic sequencing)
- Patient ID box, referral/clinical indication, panel name
- Workflow name and DNAnexus job IDs
- Human-readable filter description
- Interpretation table with pass/fail variant counts
- Links to relevant ACMG classification templates

### HELIOS (haemato-oncology)
- Similar structure to DIAS but with haemato-oncology-specific layout
- Optimised for reporting somatic and germline variants together

### URANUS (myeloid/MDS)
- MYE workbook layout compatible with ClinVar submission format
- Additional metadata fields for myelodysplastic syndromes reporting

## DNAnexus Integration (dxapp.json and src/code.sh)

The app runs on DNAnexus using a `mem1_ssd1_v2_x4` instance (4 vCPUs, SSD storage). The manifest below shows the declared inputs and outputs.

<details>
<summary><code>dxapp.json</code> — full manifest</summary>

```bash
cat dxapp.json
```

```output
{
  "name": "eggd_generate_variant_workbook",
  "title": "eggd_generate_variant_workbook",
  "summary": "Create Excel workbook from VEP annotated vcf",
  "dxapi": "1.0.0",
  "version": "2.11.1",
  "whatsNew": "* v2.0.0 Rewrite of previous app to generate xlsx file from a VEP annotated VCF(s); * v2.0.1 Bug fix to correctly treat CHROM as string values; * v2.0.2 Bug fix for ACMG report template structure; * v2.0.3 Bug fixes for issues with hyperlinks, changed app name to eggd_generate_variant_workbook; * v2.1.0 Handle VCFs from GATK gCNV and Illumina TSO500, readability tweaks to variant sheets; * v2.1.1 Bug fix for typing of numeric values in hyperlinks; * v2.2.0 Added ability to pass in non VCF files (tsvs/csvs and images) to additional sheets, optional adding of links to DECIPHER with --decipher; * v2.3.0 Added conditional colouring of cells in variant sheets, new 'basic' summary sheet;  * v2.4.0 Added handling for duplicate annotation in VEP fields (i.e. cosmic, CGC, etc..); * v2.5.0 Better parsing of CombinedVariantOutput files as additional files; * v2.6.0 Add variant counts as DNAnexus file details to the .xlsx workbook; *v2.7.0 Handle pre-split and non VEP annotated VCFs, improvements to Dias reporting templates and Excel data validation; * v2.7.1 v.2.7.0 app was accidentally published on DNAnexus before testing; so a new version is created. Everything except version number is the same as v2.7.0; * v2.8.0 Add PID box in summary, swap BS1 and BA1 in interpret table, increase the number of unlocked rows to 500; ; * v2.8.1 Update lock_sheet function to allow formatting cols/rows; * v2.8.2 Update summary sheet titles and tables * v2.9.0 Presents Uranus VCFs; * v2.10.0 Update helios summary sheet titles; * v2.10.1 Remove COSMIC mention from report_text,* v2.10.2 Update gnomAD hyperlink to v4 for build 38 * v2.11.0 Make Uranus/MYE workbooks compatible for clinvar submission; * v2.11.1 merge hotfix to v2.10.2 with features from v2.11.0",
  "authorizedUsers": [
    "org-emee_1"
  ],
  "developers":[
    "org-emee_1"
  ],
  "inputSpec": [
    {
      "name": "vcfs",
      "label": "VEP annotated vcf(s)",
      "class": "array:file",
      "optional": false,
      "help": "",
      "group": "app"
    },
    {
      "name": "additional_files",
      "label": "additional files",
      "class": "array:file",
      "optional": true,
      "help": "additional tsv/csv file(s) to read in and add to separate sheets. Sheet names may be specified with --additional_sheet_names."
    },
    {
      "name": "images",
      "label": "images",
      "class": "array:file",
      "optional": true,
      "help": "Images to write to separate sheets. Sheet names may be specified with --image_sheet_names."
    },
    {
      "name": "image_sheet_names",
      "label": "Image sheet names",
      "class": "string",
      "optional": true,
      "help": "'Names to use for image sheets, these MUST be the same number as the number of images passed and in the same order (i.e. -iimages=graph1.png -iimages=another_image.jpeg -iimage_sheet_names='myNiceGraph someImage'). If not given, this will default image1, image2...",
      "group": "generate_workbook.py"
    },
    {
      "name": "image_sizes",
      "label": "Image sizes",
      "class": "string",
      "optional": true,
      "help": "Sizes to set for images passed with --images, formatted as colon separated width:height in px, these MUST be same number as the number of images passed and in same order (i.e. -iimages=file1 -iimages=file2 -iimage_sizes='1920:1080 1000:500')",
      "group": "generate_workbook.py"
    },
    {
      "name": "exclude_columns",
      "label": "Exclude columns",
      "class": "string",
      "optional": true,
      "help": "Columns of VCF to exclude from output xlsx",
      "group": "generate_workbook.py"
    },
    {
      "name": "include_columns",
      "label": "Include columns",
      "class": "string",
      "optional": true,
      "help": "Columns of VCF to only include in output xlsx",
      "group": "generate_workbook.py"
    },
    {
      "name": "reorder_columns",
      "label": "Reorder columns",
      "class": "string",
      "optional": true,
      "help": "Set order for columns in workbook, any not specified will be appended to the end",
      "group": "generate_workbook.py"
    },
    {
      "name": "rename_columns",
      "label": "Rename columns",
      "class": "string",
      "optional": true,
      "help": "= separated key value pairs of VCF fields to rename in output xlsx",
      "group": "generate_workbook.py"
    },
    {
      "name": "types",
      "label": "Types",
      "class": "string",
      "optional": true,
      "help": "= separated key value pairs of field=type to overwrite in VCF header (i.e CSQ_gnomADg_AF=Float)",
      "group": "generate_workbook.py"
    },
    {
      "name": "filter",
      "label": "Filter",
      "class": "string",
      "optional": true,
      "help": "Filters to apply to variants",
      "group": "generate_workbook.py"
    },
    {
      "name": "keep_filtered",
      "label": "Keep filtered",
      "class": "boolean",
      "optional": true,
      "default": true,
      "help": "Determines if to keep filtered variants in separate 'excluded' tab",
      "group": "generate_workbook.py"
    },
    {
      "name": "keep_tmp",
      "label": "Keep tmp vcfs",
      "class": "boolean",
      "default": false,
      "help": "Determines if to upload the intermediate bcftools split and filtered vcfs",
      "group": "generate_workbook.py"
    },
    {
      "name": "add_samplename_column",
      "label": "Add sample name",
      "class": "boolean",
      "optional": true,
      "default": false,
      "help": "Determines if to add samplename as first column in each sheet",
      "group": "generate_workbook.py"
    },
    {
      "name": "add_comment_column",
      "label": "Add comment column",
      "class": "boolean",
      "optional": true,
      "default": false,
      "help": "Add empty comment column to end of each sheet of variants",
      "group": "generate_workbook.py"
    },
    {
      "name": "add_classification_column",
      "label": "Add classification column",
      "class": "boolean",
      "optional": true,
      "default": false,
      "help": "Add empty classification column to end of each sheet of variants",
      "group": "generate_workbook.py"
    },
    {
      "name": "add_allele_origin_column",
      "label": "Add allele_origin column",
      "class": "boolean",
      "optional": true,
      "default": false,
      "help": "Add empty allele origin column to end of each sheet of variants",
      "group": "generate_workbook.py"
    },
    {
      "name": "add_interpreted_column",
      "label": "Add interpreted column",
      "class": "boolean",
      "optional": true,
      "default": false,
      "help": "Add empty interpreted column to end of each sheet of variants",
      "group": "generate_workbook.py"
    },
    {
      "name": "add_reported_column",
      "label": "Add reported column",
      "class": "boolean",
      "optional": true,
      "default": false,
      "help": "Add empty reported column to end of each sheet of variants",
      "group": "generate_workbook.py"
    },
    {
      "name": "add_mnv_column",
      "label": "Add MNV column",
      "class": "boolean",
      "optional": true,
      "default": false,
      "help": "Add empty MNV column to end of each sheet of variants",
      "group": "generate_workbook.py"
    },
    {
      "name": "add_report_text_column",
      "label": "Add report text column to variant sheets",
      "class": "boolean",
      "optional" : true,
      "help": "If true, a report text column will be added that contains the key variant annotation in one cell"
    },
    {
      "name": "sheet_names",
      "label": "Sheet names",
      "class": "string",
      "optional": true,
      "help": "'Names to use for workbook sheets, these MUST be the same number as the number of vcfs passed and in the same order. If not given, and if there is 1 vcf passed the sheet will be named `variants`, else if multiple vcfs are passed the name prefix of the vcf will be used",
      "group": "generate_workbook.py"
    },
    {
      "name": "additional_sheet_names",
      "label": "Additional sheet names",
      "class": "string",
      "optional": true,
      "help": "Names to use for additional file sheets, if specified these MUST be the same number as the number of files passed and in the same order (e.g. `-iadditional_files=file1 -iadditional_files=file2 -iadditional_sheet_names='name_1 name_2'`). If not given, the first 31 characters of the filename will be used",
      "group": "generate_workbook.py"
    },
    {
      "name": "output_prefix",
      "label": "Output prefix",
      "class": "string",
      "optional": true,
      "help": "Prefix for naming output xlsx file",
      "group": "generate_workbook.py"
    },
    {
      "name": "freeze_column",
      "label": "freeze column",
      "class": "string",
      "optional": true,
      "help": "Optional column / row on which to freeze Excel scrolling for variant sheets (default: A2)",
      "group": "generate_workbook.py"
    },
    {
      "name": "colour_cells",
      "label": "colour cells",
      "class": "string",
      "optional": true,
      "help": "Add conditional colouring of cells for a given column, this should be specified as column:value_range:colour, where colour is a valid hex value or colour name. See readme for futher details.",
      "group": "generate_workbook.py"
    },
    {
      "name": "merge_vcfs",
      "label": "Merge",
      "class": "boolean",
      "optional": true,
      "default": false,
      "help": "Determines if to merge multiple VCFs to one sheet",
      "group": "generate_workbook.py"
    },
    {
      "name": "summary",
      "label": "Summary",
      "class": "string",
      "optional": true,
      "help": "If to include summary sheet, specify key of assay",
      "group": "generate_workbook.py"
    },
    {
      "name": "human_filter",
      "label": "Human filter",
      "class": "string",
      "optional": true,
      "help": "String to add to summary sheet with humanly readable form of the given filter string. No checking is done of this matching the actual filter(s) used.",
      "group": "generate_workbook.py"
    },
    {
      "name": "acmg",
      "label": "ACMG",
      "class": "int",
      "optional": true,
      "default": 0,
      "help": "Determines number of extra sheet(s) with ACMG reporting criteria",
      "group": "generate_workbook.py"
    },
    {
      "name": "print_columns",
      "label": "Print columns",
      "class": "boolean",
      "optional": true,
      "help": "Print all column names of all vcfs that will be output to the xlsx. Useful to identify what will be in the output to include/exclude.",
      "group": "generate_workbook.py"
    },
    {
      "name": "print_header",
      "label": "Print header",
      "class": "boolean",
      "optional": true,
      "help": "Print header of first vcf that will be output to the xlsx. Useful to identify field types in the VCF, which can be modified with --types.",
      "group": "generate_workbook.py"
    },
    {
      "name": "panel",
      "label": "String of panel information to display in summary sheet",
      "class": "string",
      "optional": true,
      "help": "",
      "group": "generate_workbook.py"
    },
    {
      "name": "clinical_indication",
      "label": "String of clinical indication to display in summary sheet",
      "class": "string",
      "optional": true,
      "help": "",
      "group": "app"
    },
    {
      "name": "additional_columns",
      "label": "additional columns",
      "class": "string",
      "optional": true,
      "help": "List of additional columns to add with hyperlinks to external resources. See readme for details.",
      "group": "app"
    },
    {
      "name": "split_hgvs",
      "label": "split hgvs",
      "class": "boolean",
      "optional" : true,
      "help": "If true, the c. and p. changes in HGVSc and HGVSp will be split out into DNA and Protein columns, without the transcript"
    },
    {
      "name": "add_raw_change",
      "label": "add raw change",
      "class": "boolean",
      "optional" : true,
      "help": "If true, will add a column named 'rawChange' with a concatenation of columns formatted as {CHROM}:g.{POS}{REF}>{ALT}"
    },
    {
      "name": "lock_sheet",
      "label": "lock_sheet",
      "class": "boolean",
      "optional" : true,
      "help": "If true, all sheets in the variant workbook are locked for dias pipeline except specific cells"
    },
    {
      "name": "af_format",
      "label": "af_format",
      "class": "string",
      "optional" : true,
      "help": "Presents the allele frequency (AF) as a decimal (0-1) or as a percent (0-100). Default is decimal. Options are 'decimal' or 'percent'."
    },
    {
      "name": "join_columns",
      "label": "Join columns",
      "class": "string",
      "optional" : true,
      "help": "User to join two columns from VCF into a new column with a seperator of choice. The header needs to be added to the include or rename if this is used"
    },
    {
      "name": "m_codes",
      "label": "M-codes file containing all valid M-codes",
      "class": "file",
      "optional" : true,
      "help": "DNAnexus file containing all valid M-codes. M-codes should be provided one per line in a .txt file."
    },
    {
      "name": "add_auto_filter",
      "label": "Add an excel auto filter to variant sheets",
      "class": "boolean",
      "optional" : true,
      "help": "If true, an auto-filter is applied to all column headers in variant sheets"
    },
    {
      "name": "sort_by",
      "label": "Sort variant sheets by column(s)",
      "class": "string",
      "optional" : true,
      "help": "List of column name:bool pairs, for columns to sort on and bool whether to sort in ascending order"
    }
  ],
  "outputSpec": [
    {
      "name": "xlsx_report",
      "label": "Excel workbook for the given vcf(s)",
      "class": "file",
      "patterns": [
        "*.xlsx"
      ],
      "help": ""
    },
    {
      "name": "tmp_vcfs",
      "label": "tmp vcfs",
      "class": "array:file",
      "optional": true,
      "help": "intermediate split and annotated vcfs, output with --keep_tmp"
    }
  ],
  "runSpec": {
    "timeoutPolicy": {
      "*": {
        "hours": 2
      }
    },
    "interpreter": "bash",
    "release": "20.04",
    "version": "0",
    "distribution": "Ubuntu",
    "file": "src/code.sh",
    "assetDepends": [
      {
        "name": "htslib",
        "project": "project-Fkb6Gkj433GVVvj73J7x8KbV",
        "folder": "/app_assets/htslib/htslib_v1.15.0",
        "version": "1.15.0"
      }
    ]
  },
  "access": {
    "project": "CONTRIBUTE",
    "allProjects": "VIEW",
    "network": [
      "*"
    ]
  },
  "regionalOptions": {
    "aws:eu-central-1": {
      "systemRequirements": {
        "*": {
          "instanceType": "mem1_ssd1_v2_x4"
        }
      }
    }
  }
}
```

</details>

Key DNAnexus configuration points from the manifest:

- **Version**: 2.11.1 (authorised to `org-emee_1`)
- **Runtime**: Ubuntu 20.04 on `mem1_ssd1_v2_x4` (AWS eu-central-1), 2-hour timeout
- **Asset dependency**: htslib v1.15.0 — provides `bcftools` (including the `+split-vep` plugin) and `bgzip`; installed from a pre-built DNAnexus asset rather than the Ubuntu repos
- **Network access**: `["*"]` — required for pulling bcftools plugin if not bundled
- **Primary output**: `xlsx_report` (required); optional `tmp_vcfs` array when `--keep_tmp` is set

## Dependencies

### System (provided via DNAnexus asset)
| Package | Version | Role |
|---------|---------|------|
| bcftools | 1.15.0 | VCF splitting (`+split-vep`) and filtering |
| bgzip / htslib | 1.15.0 | VCF compression/indexing |

### Python (vendored as wheels in `resources/home/dnanexus/packages/`)
| Package | Version | Role |
|---------|---------|------|
| pandas | 1.4.1 | DataFrame-based VCF manipulation |
| openpyxl | 3.1.2 | Excel workbook creation |
| numpy | 1.22.2 | Numeric operations |
| colour | 0.1.5 | Colour name/hex parsing for cell colouring |
| filetype | 1.1.0 | Detect image file types |
| Pillow | 9.2.0 | Image embedding in Excel |
| python-Levenshtein | 0.13.2 | String similarity (column fuzzy matching) |
| pytz | 2021.3 | Timezone handling for summary timestamps |

Wheels are installed offline from the bundled `packages/` directory by `code.sh`, so the app has no pip internet access requirement at runtime.

## Tests

The test suite lives at `resources/home/dnanexus/generate_workbook/tests/` and is run with **pytest**. Tests are organised by module:

<details>
<summary>Full test function listing</summary>

```bash
grep -h 'def test_' resources/home/dnanexus/generate_workbook/tests/test_vcf.py resources/home/dnanexus/generate_workbook/tests/test_filters.py resources/home/dnanexus/generate_workbook/tests/test_excel.py resources/home/dnanexus/generate_workbook/tests/test_columns.py resources/home/dnanexus/generate_workbook/tests/test_generate_workbook.py resources/home/dnanexus/generate_workbook/tests/test_utils.py | sed 's/    def //' | sed 's/(self.*//' | grep -v '__'
```

```output
test_sort_vcfs
test_column_names
test_parse_reference_vep
test_parse_reference_no_vep
test_only_header_parsed
test_annotated_and_not_split
test_vcf_already_split
test_vcf_not_annotated_with_vep
test_tmp_vcf_made
test_drop_columns_exclude
test_drop_columns_include
test_reorder_columns_correct_order
test_reorder_columns_no_dropped_columns
test_non_rename_columns_unaffacted
test_renamed_correctly
test_join_columns_right_not_comma
test_join_columns_right_with_other_characters
test_type_error_raised_when_too_many_separators_provided
test_type_error_raised_when_space_separator_is_provided
test_type_error_raised_when_no_separator_is_provided
test_type_error_raised_when_new_header_equal_not_provided
test_no_fail_when_vcf_doesnt_have_column_in_vcf
test_column_creation():
test_decipher_column_added():
test_decipher_column_not_added():
test_decipher_build_37():
test_decipher_links_build_38():
test_gnomad_build_37():
test_gnomad_build_38():
test_cosmic_build_37():
test_normal_df
test_empty_df
test_missing_column
test_report_text
test_type_correctly_modified
test_header_overwritten_correctly
test_filter_with_include_eq
test_filter_with_exclude_eq
test_filter_with_exclude_gt
test_combined_exclude_float_and_string
test_combined_filter_and_recover
test_attributes_are_set_correctly
test_unlock_specified_cells
test_unlock_region
test_get_cells_in_columns_returns_cells
test_get_cells_missing_column_error
test_store_list_in_sheet
test_success_upon_compliant_mcodes
test_exception_upon_noncompliant_mcodes
test_str_to_drop_down_options_too_long
test_cells_formatted_as_perc
test_chrom
test_pos
test_id
test_ref
test_alt
test_qual
test_filter
test_parsed_correct_columns_from_info_records
test_parsed_correct_gnomAD_AF_values
test_format_sample_values_are_correct
test_parsed_correct_COSMICcMuts_values
test_parsed_correct_COSMICncMuts_values
test_invalid_image_assertion
test_differing_images_and_image_sheet_names
test_differing_images_and_image_sizes
test_invalid_image_sizes
test_valid_colour_expressions
test_invalid_colour_expression
test_sort_by_input_without_colon_raises_error
test_sort_by_input_without_col_or_bool_raises_error
test_sort_by_input_with_invalid_bool_raises_error
test_valid_sort_by_input_is_returned_correctly
def test_is_numeric():
test_comma():
test_semicolon():
test_tab():
test_space():
test_mixed():
test_tsv_suffix():
test_csv_suffix():
```

</details>

The test suite covers:

| Test file | What it covers |
|-----------|---------------|
| `test_vcf.py` | VCF parsing, reference detection, VEP annotation checks, pre-split detection, column dropping/reordering/renaming, HGVS splitting, `join_columns`, hyperlink generation (gnomAD, ClinVar, COSMIC, DECIPHER), report text construction, variant sorting |
| `test_filters.py` | bcftools filter expressions (include/exclude, GT/EQ, float/string types), combined filter + recover logic, header type modification |
| `test_excel.py` | Argument validation, cell unlocking, cell region unlocking, M-code validation, dropdown string truncation, percentage formatting |
| `test_columns.py` | INFO column parsing, gnomAD AF extraction, FORMAT/SAMPLE pairing, COSMIC annotation handling |
| `test_generate_workbook.py` | Image argument validation, colour expression validation, sort-by argument parsing |
| `test_utils.py` | `is_numeric` helper, delimiter detection (comma/semicolon/tab/space/mixed), file suffix handling |

CI runs the full suite on every push via `.github/workflows/pytest.yml`.

## Notable Implementation Details

### Duplicate VEP annotation handling
VEP can annotate variants with multiple entries from sources like COSMIC or CGC, producing duplicate transcript rows. `splitColumns.unique_vep()` removes these duplicates while preserving row order, preventing inflated variant counts in the workbook.

### TSO500 special handling
When `CombinedVariantOutput.tsv` (Illumina TSO500 pipeline) is passed as an additional file, a custom parser extracts only the relevant sections: TMB (Tumour Mutational Burden), MSI (Microsatellite Instability), and Gene Amplifications. `MetricsOutput.tsv` is parsed separately to pull sample-level QC metrics.

### Excel cell locking (`--lock_sheet`)
Sheet protection locks all cells by default, then selectively unlocks:
- Up to **500 rows × 200 columns** of annotation columns (comment, classification, allele origin, interpreted, reported, MNV)
- The auto-filter control cells if `--add_auto_filter` is also set
This pattern allows analysts to enter classifications while preventing accidental modification of variant data.

### Conditional cell colouring (`--colour`)
Expressions follow the format `column:condition:colour`:
- `VF:>=0.9:green` — cells where VF ≥ 0.9 are coloured green
- `VF:<0.6&>=0.3:orange` — AND condition (must not mix with `|` in the same expression)
- `VF:0.1|0.5:blue` — exact-value OR match

The `colour` Python library is used to validate and convert colour names to hex. openpyxl `PatternFill` applies the fill to matching cells.

### Reference build detection
`vcf.parse_reference()` scans the VCF header for `GRCh38` or `GRCh37` patterns. The detected build drives:
- gnomAD hyperlink version (v4 for GRCh38, v2 for GRCh37)
- DECIPHER column availability (GRCh38 only)
- MasterMind NC chromosome mapping

### String encoding safety
`format_strings()` URL-decodes percent-encoded values (common in VEP output) and re-encodes to UTF-8, protecting against encoding errors when openpyxl writes cells.

## Typical Usage Examples

### Minimal — single VCF to Excel

    python3 generate_workbook.py \
      --vcfs sample.vcf.gz

### DIAS diagnostic report

    python3 generate_workbook.py \
      --vcfs sample.vcf.gz \
      --filter "-e 'CSQ_gnomADe_AF > 0.02'" \
      --types CSQ_gnomADe_AF=Float \
      --keep \
      --summary dias \
      --panel "Epilepsy v2.1" \
      --clinical_indication "R59 - Childhood onset epilepsy" \
      --acmg 2 \
      --split_hgvs \
      --add_classification_column \
      --add_allele_origin_column \
      --lock_sheet \
      --colour "VF:>=0.9:green" "VF:<0.9&>=0.6:orange" \
      --output sample_report

### Merging multiple VCFs with samplename tracking

    python3 generate_workbook.py \
      --vcfs sampleA.vcf.gz sampleB.vcf.gz \
      --merge \
      --add_name \
      --sort_by CHROM:True POS:True \
      --output cohort_report

### Inspect available columns before running

    python3 generate_workbook.py --vcfs sample.vcf.gz --print_columns
    python3 generate_workbook.py --vcfs sample.vcf.gz --print_header

## Version History Summary

| Version | Notable change |
|---------|---------------|
| v2.0.0 | Full rewrite to generate xlsx from VEP-annotated VCF |
| v2.1.0 | GATK gCNV and Illumina TSO500 handling |
| v2.2.0 | Non-VCF additional files (TSV/CSV, images); DECIPHER links |
| v2.3.0 | Conditional cell colouring; basic summary sheet |
| v2.4.0 | Duplicate VEP annotation removal (COSMIC, CGC) |
| v2.5.0 | Improved CombinedVariantOutput.tsv parsing |
| v2.6.0 | Variant counts stored as DNAnexus file metadata |
| v2.7.0 | Pre-split VCF support; non-VEP VCF handling; improved DIAS templates |
| v2.8.0 | PID box in summary; reordered ACMG criteria; 500-row sheet unlock |
| v2.9.0 | URANUS pipeline support |
| v2.10.0 | HELIOS summary sheet updates |
| v2.10.2 | gnomAD hyperlinks updated to v4 for GRCh38 |
| v2.11.0 | URANUS/MYE workbooks compatible for ClinVar submission |
| v2.11.1 | Merge of v2.10.2 hotfix into v2.11.0 (current) |
