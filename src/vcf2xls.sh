#!/bin/bash

set -euxo pipefail

_download_files () {
    # downloads input files with xargs parallel for speed
    local dir=$1
    local file_list=$2
    xargs -n1 -P"${CPU}" dx download --output "$dir" "$file_list"
}


_dias_report_setup () {
    # function to handle parsing values and reading
    # manifest / g2t etc. for Dias sampels
    mark-section "Parsing values"
    if [ ! -z ${list_panel_names_genes+x} ]; then
        echo "Value of list_panel_names_genes: '$list_panel_names_genes'"
    fi

    single_genes=""

    # get single genes from manifest/given list
    if [ ! -z ${list_panel_names_genes+x} ]; then
        default_IFS=$IFS
        IFS=";"
        single_genes=$(for ele in $list_panel_names_genes; do if [[ $ele =~ ^_ ]]; then echo $ele | sed s/_//; fi; done)
        IFS=$default_IFS
    else
        sample_in_manifest=$(grep $sample_id BioinformaticManifest)
        single_genes=$(awk -v id=$sample_id '{if(id == $1 && $2 ~ /^_/){print $4}}' BioinformaticManifest)

        # also check if sample is present in the bioinformatic manifest
        if [[ ! $sample_in_manifest ]]; then
            echo "${sample_id} is not present in the bioinformatic manifest, exiting..."
            exit
        fi
    fi

    # Check if single genes are present in g2t
    if [[ $single_genes ]]; then
        for gene in $single_genes; do
            if ! grep -q $gene nirvana_genes2transcripts; then
                echo "${gene} not found in nirvana g2t"
                exit
            fi
        done
    fi

    # Boolean to detect if workflow id has been found
    found_workflow_id=false

    # get job id creating the annotated vcf
    annotation_job_id=$(dx describe --json ${vcfs[0]} | jq -r '.createdBy.job')
    annotation_job_name=$(dx describe --json ${annotation_job_id} | jq -r '.executableName')

    workflow_id=$(dx describe --json ${annotation_job_id} | jq -r '.parentAnalysis')
    workflow_name=$(dx describe --json ${workflow_id} | jq -r '.executableName')


    # Placeholder text if the workflow id is not found
    if [ -z $annotation_job_name ]; then
        analysis_name="No workflow id found for this report."
    fi
    if [ -z $workflow_name ]; then
        workflow_id="This report was probably generated for development purposes, do not use for clinical reporting"
    fi


    # get read stats from flagstat file
    total_nb_reads=$(grep total inputs/$flagstat_file_name | cut -d+ -f1)
    nb_duplicates_reads=$(grep duplicate inputs/$flagstat_file_name | cut -d+ -f1)
    nb_aligned_reads=$(grep "mapped (" inputs/$flagstat_file_name | cut -d+ -f1)
    nb_usable_reads=$(( $nb_aligned_reads - $nb_duplicates_reads ))


    project_id=$DX_PROJECT_CONTEXT_ID

    version=0
    matching_files=1
    if [ -z "$output_prefix" ]; then
        output_name=${sample_id}_${version}
    else
        output_name="${output_prefix}"
    fi

    # Tiny chance of race conditions leading to two files with the same name here
    while [ $matching_files -ne 0 ]; do
        version=$((version+1))
        output_name=${sample_id}_${version}*.xls
        matching_files=$(dx find data --path ${project_id}:/ --name $output_name --brief | wc -l)
    done;

    # Add text to report name if workflow id hasn't been found
    if [ $found_workflow_id = true ]; then
        output_name="${sample_id}_${version}.xls"
    else
        output_name="${sample_id}_${version}_FOR_DEV_USE_ONLY.xls"
    fi
}


_panel_filter () {
    # Filters with bedtools intersect if panel bed file given
    vcf=$1
    echo "Filtering ${vcf} agaist ${panel_bed_name}"

    # get nicely formatted name of bed file for report
    panel=$(sed 's/_\([0-9]\{2\}\)bp_b\([0-9]\{2\}\).bed//g' <<< $bed)  # remove generate_bed suffix
    panel=${bed//_/}
    panel=${bed//&&/ & }

    bedtools intersect -header -a $vcf -b inputs/$panel_bed_name > vcfs/$vcf
}


main() {
    echo "Value of vcf(s): $vcfs"
    echo "Value of flagstat_file: $flagstat_file"

    # total cores available for cmds that make use of them
    CPU=$(grep -c ^processor /proc/cpuinfo)

    mark-section "Downloading inputs"
    mkdir vcfs
    _download_files --output "vcfs/" "$vcfs"

    if [ "$assay" == "dias" ]; then
        # download and do dias specific things
        _download_files "$flagstat_file" "$genepanels_file" \
            "$bioinformatic_manifest" "$nirvana_genes2transcripts"

        _dias_report_setup
    fi

    mkdir -p /home/dnanexus/out/xls_reports

    if [ "$panel_bed" ]; then
        # filtering against bed file, move all the vcfs then add back the
        # filtered ones to vcfs dir
        dx download "$panel_bed"

        mkdir full_vcfs
        mv vcfs/* unfiltered_vcfs/

        for vcf in unfiltered_vcfs/*; do
            _panel_filter "$vcf"
        done
    fi

    # build string of input arguments
    optional_args=""
    if [ "$clinical_indication" ]; then args+="--clinical_indication ${clinical_indication}"; fi
    if [ "$flagstat_file" ]; then args+="--usable_reads ${nb_usable_reads} "; fi
    if [ "$exclude_cols" ]; then args+="--exclude ${exclude_cols} "; fi
    if [ "$include_cols" ]; then args+="--include ${include_cols} "; fi
    if [ "$reorder_cols" ]; then args+="--reorder ${reorder_cols} "; fi
    if [ "$rename_cols" ]; then args+="--rename ${rename_cols} "; fi
    if [ "$flagstat_file" ]; then args+="--reads ${total_nb_reads} "; fi
    if [ "$output_prefix" ]; then args+="--sample ${output_prefix} "; fi
    if [ "$workflow" ]; then args+="--workflow ${workflow_id} ${analysis_name} "; fi
    if [ "$analysis" ]; then args+="--analysis ${analyis} "; fi
    if [ "$summary" ]; then args+="summary ${summary} "; fi
    if [ "$sheets" ]; then args+="--sheets ${sheets} "; fi
    if [ "$panel" ]; then args+="--panel ${panel} "; fi
    if [ "$add_name" ]; then args+="--add_name "; fi
    if [ "$filter" ]; then args+="--filter "; fi
    if [ "$merge" ]; then args+="--merge "; fi
    if [ "$keep" ]; then args+="--keep "; fi

    python3 vcf2xls.py --vcfs vcfs/* --output $output_name \
        --out_dir "/home/dnanexus/out/xls_reports" $optional_args


    echo "Output name: $output_name"

    output_file=$(dx upload /home/dnanexus/out/xls_reports/* --brief)
    dx-jobutil-add-output xls_report "$output_file" --class=file
}
