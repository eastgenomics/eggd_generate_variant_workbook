#!/bin/bash

set -exo pipefail

_dias_report_setup () {
    # function to handle parsing values and reading
    # manifest / g2t etc. for Dias sampels
    mark-section "Parsing values"

    # get job id creating the annotated vcf
    annotation_job_id=$(dx describe --json ${vcfs[0]} | jq -r '.createdBy.job')
    annotation_job_name=$(dx describe --json ${annotation_job_id} | jq -r '.executableName')

    workflow_id=$(dx describe --json ${annotation_job_id} | jq -r '.parentAnalysis')
    workflow_name=$(dx describe --json ${workflow_id} | jq -r '.executableName')

    # Placeholder text if the workflow id is not found
    if [ -z "$annotation_job_name" ]; then
        analysis_name="No workflow id found for this report."
    fi
    if [ -z "$workflow_name" ]; then
        workflow_id="This report was probably generated for development purposes, do not use for clinical reporting"
    fi

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

    # Add text to report name if workflow id hasn't been found
    if [[ $workflow_id ]]; then
        output_prefix="${sample_id}_${version}"
    else
        output_prefix="${sample_id}_${version}_FOR_DEV_USE_ONLY"
    fi
}



main() {
    echo "Value of vcf(s): ${vcfs[*]}"
    export BCFTOOLS_PLUGINS=/usr/local/libexec/bcftools/

    mark-section "Downloading inputs"
    mkdir vcfs
    dx-download-all-inputs --parallel
    find ~/in/vcfs -type f -name "*" -print0 | xargs -0 -I {} mv {} ~/vcfs

    if [ "$assay" == "dias" ]; then
        # do dias specific things
        _dias_report_setup
    fi

    mkdir -p /home/dnanexus/out/xls_reports && sudo chmod 757 /home/dnanexus/out/xls_reports

    mark-section "Installing packages"
    # install required python packages
    python3 -m pip install --no-index --no-deps packages/*

    echo "keep passed: ${keep}"
    echo "merge passed: ${merge}"

    mark-section "Generating workbook"
    # build string of input arguments
    args=""
    if [ "$clinical_indication" ]; then args+="--clinical_indication ${clinical_indication}"; fi
    if [ "$flagstat_file" ]; then args+="--usable_reads ${nb_usable_reads} "; fi
    if [ "$exclude_cols" ]; then args+="--exclude ${exclude_cols} "; fi
    if [ "$include_cols" ]; then args+="--include ${include_cols} "; fi
    if [ "$reorder_cols" ]; then args+="--reorder ${reorder_cols} "; fi
    if [ "$rename_cols" ]; then args+="--rename ${rename_cols} "; fi
    if [ "$print_columns" ]; then args+="--print_columns "; fi
    if [ "$print_header" ]; then args+="--print_header"; fi
    if [ "$output" ]; then args+="--sample ${output} "; fi
    if [ "$output" ]; then args+="--output ${output} "; fi
    if [ "$workflow" ]; then args+="--workflow ${workflow_id} ${analysis_name} "; fi
    if [ "$analysis" ]; then args+="--analysis ${analyis} "; fi
    if [ "$summary" ]; then args+="summary ${summary} "; fi
    if [ "$sheets" ]; then args+="--sheets ${sheets} "; fi
    if [ "$filter" ]; then args+="--filter ${filter} "; fi
    if [ "$panel" ]; then args+="--panel ${panel} "; fi
    if [ "$add_name" == true ]; then args+="--add_name "; fi
    if [ "$merge" == true ]; then args+="--merge "; fi
    if [ "$keep" == true ]; then args+="--keep "; fi

    /usr/bin/time -v python3 vcf2xls/vcf2xls.py --vcfs vcfs/* --out_dir "/home/dnanexus/out/xls_reports" $args

    mark-section "Uploading output"
    output_file=$(dx upload /home/dnanexus/out/xls_reports/* --brief)
    dx-jobutil-add-output xls_report "$output_file" --class=file
}
