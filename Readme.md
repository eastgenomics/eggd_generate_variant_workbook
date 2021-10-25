<!-- dx-header -->
# vcf2xls_nirvana (DNAnexus Platform App)

## What does this app do?

Generate a report from an annotated vcf

## What are typical use cases for this app?

This app may be executed as a standalone app.

## What data are required for this app to run?

Inputs:
- list_panel_names_genes: List of panels/genes to get data from, used in reanalysis (optional) 
- annotated_vcf: Nirvana annotated vcf
- raw_vcf: Vcf before annotation
- sample_coverage_file: Sample level coverage file
- sample_coverage_index: Index of sample level coverage file
- flagstat_file: Flagstat file
- genepanels_file: Gene panels file (in 001_Reference)
- bioinformatic_manifest: Bioinformatic manifest (in 001_Reference)
- exons_nirvana: Dump of exons from nirvana 2.0.10 (in 001_Reference)
- nirvana_genes2transcripts: Genes2transcripts file (in 001_Reference)
- panel_bed: flanked bed file used to filter vcf (optional)

Example:
```
dx run app-vcf2xls_nirvana/1.5.2 -iannotated_vcf=X210333_markdup_recalibrated_Haplotyper.refseq_nirvana_2010.annotated.vcf -iraw_vcf="X210333_markdup_recalibrated_Haplotyper.vcf.gz" -isample_coverage_file=X210333_markdup.nirvana_2010_5bp.gz -isample_coverage_index=X210333_markdup.nirvana_2010_5bp.gz.tbi -iflagstat_file=X210333_markdup.flagstat -igenepanels=001_Reference:/dynamic_files/gene_panels/gemini_panels_200522 -ibioinformatic_manifest=001_Reference:/dynamic_files/BioinformaticManifest/BioinformaticManifest_200819 -iexons_nirvana=001_Reference:/annotation/b37/exons_nirvana2010_no_PAR_Y.tsv -inirvana_genes2transcripts=001_Reference:/dynamic_files/nirvana_genes2transcripts/nirvana_genes2transcripts_2010_200728
```

## What does this app output?
This app outputs an xls report.

This is the source code for an app that runs on the DNAnexus Platform.
For more information about how to run or modify it, see
https://documentation.dnanexus.com/.

#### This app was made by EMEE GLH