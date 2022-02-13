#!/bin/bash

set -euxo pipefail

_download_files () {
    # downloads input files with xargs parallel for speed
    local file_list=$1
    xargs -n1 -P"${CPU}" "$file_list"
}

_dias_report () {
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

    # Placeholder text if the workflow id is not found
    analysis_name="No workflow id found for this report."
    workflow_id="This report was probably generated for development purposes, do not use for clinical reporting"

    # get job id creating the gnomad annotated vcf
    gnomad_annotation_job_id=$(dx describe --delim "_" $annotated_vcf_name | grep job- | cut -d_ -f2)
    # get file id of vcf annotator input raw vcf
    nirvana_annotated_vcf_id=$(dx describe --delim "_" $gnomad_annotation_job_id | grep _dest_vcf | cut -d= -f2)

    # get workflow id and analysis name of annotated vcf
    if dx describe --delim "_" $nirvana_annotated_vcf_id | grep job- ; then
        job_id=$(dx describe --delim "_" $nirvana_annotated_vcf_id | grep job- | cut -d_ -f2)
        analysis=$(dx describe --delim "_" $job_id)

        if dx describe --delim "_" $job_id | grep Root ; then
            analysis_id=$(dx describe --delim "_" $job_id | grep Root | cut -d_ -f2)
            workflow=$(dx describe --delim "_" $analysis_id)

            if dx describe --delim "_" $analysis_id | grep Workflow ; then
                workflow_id=$(dx describe --delim "_" $analysis_id | grep Workflow | cut -d_ -f2)
                analysis_name=$(dx describe --name $analysis_id)
                found_workflow_id=true
            fi
        fi
    fi

    project_id=$DX_PROJECT_CONTEXT_ID

    version=0
    matching_files=1
    if [ -z "$output_prefix" ]; then
        output_name=${sample_id}_${version}
    else
        output_name="${output_prefix}"
    fi
}

_panel_filter () {
    # Filters with bedtools intersect if panel bed file given
    
    # filter annotated vcf to include regions in the flanked panel bed using bedtools
    if [ ! -z ${panel_bed+x} ]; then
        dx download "$panel_bed" -o inputs/
        echo $panel_bed_name

        # If panel bed is provided, filter the vcf
        bedtools intersect -header -a inputs/annotated_vcf -b inputs/$panel_bed_name > inputs/sliced_annotated_vcf
    else
        # Create sliced annotated vcf to be the same as the annotated vcf if the bed is not provided
        echo "VCF not filtered as panel bed not provided"
        cp inputs/annotated_vcf inputs/sliced_annotated_vcf
    fi
}


main() {
    echo "Value of raw_vcf: '$raw_vcf'"
    echo "Value of annotated_vcf: '$annotated_vcf'"
    echo "Value of sample_coverage_file: '$sample_coverage_file'"
    echo "Value of sample_coverage_index: '$sample_coverage_index'"
    echo "Value of flagstat_file: '$flagstat_file'"

    # total cores available for cmds that make use of them
    CPU=$(grep -c ^processor /proc/cpuinfo)

    mark-section "Downloading inputs"
    _download_files "${raw_vcf} ${annotated_vcf} ${sample_coverage_file}" \
        "${sample_coverage_index} ${flagstat_file} ${genepanels_file}" \
        "${bioinformatic_manifest} ${nirvana_genes2transcripts}"


    mkdir -p /home/dnanexus/out/xls_reports

    if [ "$summary" == "dias" ]; then
        # generating report for Dias workflow, call function to get appropriate things
        # for Dias report
        _dias_report
    fi


    if [ "$panel_bed" ]; then
        # filter with bed file if given
        _panel_filter
    fi

    if [ "$flagstat_file" ]; then
        # get read stats from flagstat file
        total_nb_reads=$(grep total inputs/$flagstat_file_name | cut -d+ -f1)
        nb_duplicates_reads=$(grep duplicate inputs/$flagstat_file_name | cut -d+ -f1)
        nb_aligned_reads=$(grep "mapped (" inputs/$flagstat_file_name | cut -d+ -f1)
        nb_usable_reads=$(expr $nb_aligned_reads - $nb_duplicates_reads)
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
    if [ "$workflow" ]; then args+="--workflow ${workflow} "; fi
    if [ "$analysis" ]; then args+="--analysis ${analyis} "; fi
    if [ "$summary" ]; then args+="summary ${summary} "; fi
    if [ "$sheets" ]; then args+="--sheets ${sheets} "; fi
    if [ "$panel" ]; then args+="--panel ${panel} "; fi
    if [ "$add_name" ]; then args+="--add_name "; fi
    if [ "$filter" ]; then args+="--filter "; fi
    if [ "$merge" ]; then args+="--merge "; fi
    if [ "$keep" ]; then args+="--keep "; fi

    python3 vcf2xls.py --vcfs $annotated_vcfs --output $output_name \
        --out_dir "/home/dnanexus/out/xls_reports" $optional_args

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

    echo "Output name: $output_name"

    output_file=$(dx upload /home/dnanexus/out/xls_reports/* --brief)
    dx-jobutil-add-output xls_report "$output_file" --class=file
}
