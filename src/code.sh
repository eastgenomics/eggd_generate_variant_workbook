#!/bin/bash

set -exo pipefail

_dias_report_setup () {
    # function to handle parsing values and reading
    # manifest / g2t etc. for Dias sampels
    mark-section "Parsing values"

    # get job id creating the annotated vcf
    vcf=$(awk -F'"' '{print $4}' <<< "${vcfs[0]}")
    analysis_id=$(dx describe --json ${vcf} | jq -r '.createdBy.job')

    if [ "$analysis_id" != "null" ]; then
        workflow_id=$(dx describe --json "${analysis_id}" | jq -r '.parentAnalysis')

        if [ "$workflow_id" != "null" ]; then
            workflow_name=$(dx describe --json "${workflow_id}" | jq -r '.executableName')
        fi
    fi

    project_id=$DX_PROJECT_CONTEXT_ID
    job_id=$DX_JOB_ID

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
    python3 -m pip install --no-index --no-deps --user packages/*

    echo "keep passed: ${keep}"
    echo "merge passed: ${merge}"

    # build string of input arguments
    mark-section "Building arguments"
    args=""
    if [ "$clinical_indication" ]; then args+="--clinical_indication ${clinical_indication} "; fi
    if [ "$exclude_columns" ]; then args+="--exclude ${exclude_columns} "; fi
    if [ "$include_columns" ]; then args+="--include ${include_columns} "; fi
    if [ "$reorder_columns" ]; then args+="--reorder ${reorder_columns} "; fi
    if [ "$rename_columns" ]; then args+="--rename ${rename_columns} "; fi
    if [ "$add_samplename_column" == true ]; then args+="--add_name "; fi
    if [ "$sheet_names" ]; then args+="--sheets ${sheet_names} "; fi
    if [ "$print_columns" ]; then args+="--print_columns "; fi
    if [ "$summary" ]; then args+="--summary ${summary} "; fi
    if [ "$keep_filtered" == true ]; then args+="--keep "; fi
    if [ "$keep_tmp" == true ]; then args+="--keep_tmp "; fi
    if [ "$print_header" ]; then args+="--print_header "; fi
    if [ "$merge_vcfs" == true ]; then args+="--merge "; fi
    if [ "$output" ]; then args+="--sample ${output} "; fi
    if [ "$output" ]; then args+="--output ${output} "; fi
    if [ "$workflow" ]; then args+="--workflow ${workflow_name} ${workflow_id} "; fi
    if [ "$job_id" ]; then args+="--job_id ${job_id} "; fi
    if [ "$types" ]; then args+="--types ${types} "; fi
    if [ "$panel" ]; then args+="--panel ${panel} "; fi

    args+="--out_dir /home/dnanexus/out/xlsx_reports "

    mark-section "Generating workbook"
    if [ "$filter" ]; then
        # adding the filter variable to the args string breaks everything as it is a single
        # string with spaces and quoting in bash is some black magic that makes no sense,
        # keeping this separate works so ¯\_(ツ)_/¯
        /usr/bin/time -v python3 generate_workbook/generate_workbook.py --vcfs vcfs/* $args --filter "${filter}"
    else
        /usr/bin/time -v python3 generate_workbook/generate_workbook.py --vcfs vcfs/* $args
    fi

    mark-section "Uploading output"
    output_xlsx=$(dx upload /home/dnanexus/out/xlsx_reports/* --brief)
    dx-jobutil-add-output xlsx_report "$output_xlsx" --class=file

    if [ "$keep_tmp" == true ]; then
        # upload intermediary vcfs
        tmp_vcfs=$(find . \( -name "*.filter.vcf.gz" -o -name "*.split.vcf.gz" \) | tr '\n' ' ')
        for file in $tmp_vcfs; do
            id=$(dx upload "${file}" --brief)
            dx-jobutil-add-output tmp_vcfs "$id" --class=array:file
        done
    fi
}
