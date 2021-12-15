<!-- dx-header -->
# vcf2xls_nirvana (DNAnexus Platform App)

## What does this app do?

Generate a report from an annotated vcf

## What are typical use cases for this app?

This app may be executed as a standalone app.

## What data are required for this app to run?

Inputs:

- list_panel_names_genes: List of panels/genes to get data from, used in reanalysis (optional)
- annotations: annotations present in INFO field in annotated vcf to add to the report (optional)
- annotated_vcf: Nirvana annotated vcf
- raw_vcf: Vcf before annotation
- sample_coverage_file: Sample level coverage file
- sample_coverage_index: Index of sample level coverage file
- flagstat_file: Flagstat file
- genepanels_file: Gene panels file (in 001_Reference)
- bioinformatic_manifest: Bioinformatic manifest (in 001_Reference)
- nirvana_genes2transcripts: Genes2transcripts file (in 001_Reference)
- panel_bed: flanked bed file used to filter vcf (optional)

Example:

```bash
dx run app-vcf2xls_nirvana/1.6.0 -iannotated_vcf=${annotated_vcf} -iraw_vcf=${raw_vcf} -isample_coverage_file=${nirvana_coverage.gz} -isample_coverage_index=${nirvana_coverage.gz.tbi} -iflagstat_file=${samtools_file.flagstat} -igenepanels=001_Reference:/dynamic_files/gene_panels/${genepanels.tsv} -ibioinformatic_manifest=001_Reference:/dynamic_files/BioinformaticManifest/${bioinformatic_manifest.tsv} -inirvana_genes2transcripts=001_Reference:/dynamic_files/nirvana_genes2transcripts/${g2t.tsv} [-ilist_panel_names_genes="$reanalysis_panel_string" -iannotations="$name_of_annotations"]
```

## What does this app output?
This app outputs an xls report.

This is the source code for an app that runs on the DNAnexus Platform.
For more information about how to run or modify it, see
https://documentation.dnanexus.com/.

#### This app was made by EMEE GLH