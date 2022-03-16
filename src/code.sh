#!/bin/bash

set -exo pipefail

_dias_report_setup () {
    # function to handle parsing values and reading
    # manifest / g2t etc. for Dias sampels
    mark-section "Parsing values"

    # get job id creating the annotated vcf
    analysis_id=$(dx describe --json ${vcfs[0]} | jq -r '.createdBy.job')
    analysis_name=$(dx describe --json ${analysis_id} | jq -r '.executableName')

    workflow_id=$(dx describe --json ${analysis_id} | jq -r '.parentAnalysis')
    workflow_name=$(dx describe --json ${workflow_id} | jq -r '.executableName')

    project_id=$DX_PROJECT_CONTEXT_ID

    # Tiny chance of race conditions leading to two files with the same name here
    version=0
    matching_files=1
    if [ -z "$output_prefix" ]; then
        while [ $matching_files -ne 0 ]; do
            version=$((version+1))
            output_name="${sample_id}_${version}*"
            matching_files=$(dx find data --path "${project_id}":/ --name "$output_name" --brief | wc -l)
        done;
    fi
}


main() {
    echo "Value of vcf(s): ${vcfs[*]}"
    export BCFTOOLS_PLUGINS=/usr/local/libexec/bcftools/

    mark-section "Downloading inputs"
    mkdir vcfs
    dx-download-all-inputs --parallel
    find ~/in/vcfs -type f -name "*" -print0 | xargs -0 -I {} mv {} ~/vcfs

    if [ "$summary" == "dias" ]; then
        # do dias specific things
        _dias_report_setup
    fi

    mkdir -p /home/dnanexus/out/xlsx_reports && sudo chmod 757 /home/dnanexus/out/xlsx_reports

    mark-section "Installing packages"
    # install required python packages
    python3 -m pip install --no-index --no-deps packages/*

    echo "keep passed: ${keep}"
    echo "merge passed: ${merge}"

    mark-section "Generating workbook"
    # build string of input arguments
    args=""
    if [ "$clinical_indication" ]; then args+="--clinical_indication ${clinical_indication} "; fi
    if [ "$exclude_columns" ]; then args+="--exclude ${exclude_columns} "; fi
    if [ "$include_columns" ]; then args+="--include ${include_columns} "; fi
    if [ "$reorder_columns" ]; then args+="--reorder ${reorder_columns} "; fi
    if [ "$rename_columns" ]; then args+="--rename ${rename_columns} "; fi
    if [ "$print_columns" ]; then args+="--print_columns "; fi
    if [ "$print_header" ]; then args+="--print_header "; fi
    if [ "$output" ]; then args+="--sample ${output} "; fi
    if [ "$output" ]; then args+="--output ${output} "; fi
    if [ "$workflow" ]; then args+="--workflow ${workflow_name} ${workflow_id} "; fi
    if [ "$analysis" ]; then args+="--analysis ${analysis_name} ${analysis_id} "; fi
    if [ "$summary" ]; then args+="--summary ${summary} "; fi
    if [ "$sheet_names" ]; then args+="--sheets ${sheet_names} "; fi
    if [ "$filter" ]; then args+="--filter ${filter} "; fi
    if [ "$types" ]; then args+="--types ${types} "; fi
    if [ "$panel" ]; then args+="--panel ${panel} "; fi
    if [ "$add_samplename_column" == true ]; then args+="--add_name "; fi
    if [ "$merge_vcfs" == true ]; then args+="--merge "; fi
    if [ "$keep_tmp" == true ]; then args+="--keep_tmp "; fi
    if [ "$keep_filtered" == true ]; then args+="--keep "; fi

    /usr/bin/time -v python3 generate_workbook/generate_workbook.py --vcfs vcfs/* --out_dir "/home/dnanexus/out/xlsx_reports" $args

    mark-section "Uploading output"
    output_xlsx=$(dx upload /home/dnanexus/out/xlsx_reports/* --brief)
    dx-jobutil-add-output xlsx_report "$output_xlsx" --class=file

    if [ "$keep_tmp" == true ]; then
        tmp_vcfs=$(find . \( -name "*.filter.vcf.gz" -o -name "*.split.vcf.gz" \))
        uploaded_vcfs=$(dx upload "${tmp_vcfs}" --brief)
        dx-jobutil-add-output tmp_vcfs "$uploaded_vcfs" --class=array:file
    fi
}
