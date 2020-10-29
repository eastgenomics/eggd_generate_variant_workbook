# dnanexus_vcf2xls_nirvana Developer Readme

## Internals:

The subroutines are:

- `check_sry()`: check the presence of the sry region
- `check_vcf_integrity_after_annotation()`: Runs vcf_integrity_check.py between the annotated vcf and the "raw" vcf
- `find_sample_name()`: Get the sample name using the vcf_file variable
- `fill_summary_sheet()`: For every transcript, fill the appropriate QC sheets
- `analyse_vcf_file()`: Parse variants in the vcf file
- `setup_workbook()`: Create the xls and setup the formats
- `gemini_af()`: Go through the gemini freq file to find the expected freq for a given variant
- `external_af()`: Get the freq in the exac, NHLBI Exome Sequencing Project, 1000 genomes project for given variant
- `write_variant()`: it puts the variants
- `low_AF_variant()`: Indicates if a variant is rare ?
- `pure_cpos()`: Gets the position in the HGVSc
- `setup_worksheets()`: Add given sheets to the workbook
- `usage()`: Prints usage
- `add_worksheet()`: Put headers and static text according to sheet name
- `gene_list_to_transcript_list()`: Return hash of PANEL2gene in that panel? So perlish
- `readin_manifest()`: Using the manifest, return gene in the panel?
- `worksheet_cell_width()`: Set width
- `worksheet_write()`: Write given data in given cell
- `worksheet_merge_cells_and_write()`: 
- `one2three()`: Convert AA into 3 characters AA
- `readin_panels_n_manifest()`: Create hashes of panels name to genes
- `parameter_panels2genes()`: Return hash of genes to transcript given genes
