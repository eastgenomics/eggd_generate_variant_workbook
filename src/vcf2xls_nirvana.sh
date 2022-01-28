#!/bin/bash

set -euxo pipefail

main() {
    if [ ! -z ${list_panel_names_genes+x} ]; then
        echo "Value of list_panel_names_genes: '$list_panel_names_genes'"
    fi

    if [ ! -z ${annotations+x} ]; then
        echo "Value of annotations: '$annotations'"
        parsed_annotations=""
        default_IFS=$IFS
        IFS=";"

        # parse annotations and for custom renamed annotations i.e. NEW_NAME:=NAME, keep only the new_name bit
        for ele in $annotations; do
            if [[ $ele =~ .*:=.* ]]; then
                annotation_to_add=$(echo $ele | sed s/:=.*//)
            else
                annotation_to_add=$(echo $ele)
            fi

            parsed_annotations+="$annotation_to_add;"
        done

        IFS=$default_IFS
        echo "Adjusted annotations names for vcf2xls: '$parsed_annotations'"
    fi

    echo "Value of raw_vcf: '$raw_vcf'"
    echo "Value of annotated_vcf: '$annotated_vcf'"
    echo "Value of sample_coverage_file: '$sample_coverage_file'"
    echo "Value of sample_coverage_index: '$sample_coverage_index'"
    echo "Value of flagstat_file: '$flagstat_file'"

    # Download sample inputs
    mkdir inputs

    dx download "$annotated_vcf" -o inputs/
    echo $annotated_vcf_name
    dx download "$raw_vcf" -o inputs/
    echo $raw_vcf_name
    dx download "$sample_coverage_file" -o inputs/
    echo $sample_coverage_file_name
    dx download "$sample_coverage_index" -o inputs/
    echo $sample_coverage_index_name
    dx download "$flagstat_file" -o inputs/
    echo $flagstat_file_name

    cd inputs

    if [[ $annotated_vcf_name == *.gz ]]; then
        gunzip -c $annotated_vcf_name > annotated_vcf
    else
        mv $annotated_vcf_name annotated_vcf
    fi

    cd ..

    # get sample id from vcf file name
    sample_id=$(grep -oP "^[a-zA-Z0-9]*" <<< $annotated_vcf_prefix)
    echo $sample_id

    # Download dynamic reference files
    dx download "$genepanels_file" -o genepanels
    dx download "$bioinformatic_manifest" -o BioinformaticManifest
    dx download "$nirvana_genes2transcripts" -o nirvana_genes2transcripts

    mkdir -p /home/dnanexus/out/xls_reports

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

    # filter annotated vcf to include regions in the flanked panel bed using bedtools
    dx download "project-Fkb6Gkj433GVVvj73J7x8KbV:file-G1vB3JQ433Gx5GJf1FZKYxbF" -o bedtools
    chmod a+x bedtools
    export PATH=$PATH:/home/dnanexus

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

    # Boolean to detect if workflow id has been found
    found_workflow_id=false

    # Placeholder text if the workflow id is not found
    analysis_name="No workflow id found for this report."
    workflow_id="This report was probably generated for development purposes, do not use for clinical reporting"

    # get file id of vcf annotator input raw vcf
    raw_vcf_job_id=$(dx describe --json $raw_vcf_name | jq-linux64 -r .createdBy.job)

    # get workflow id and analysis name of nirvana annotated vcf
    if dx describe --json $raw_vcf_job_id | jq-linux64 -e 'has("rootExecution")'  ; then
        analysis_id=$(dx describe --json $raw_vcf_job_id | jq-linux64 -r '.rootExecution')

        if dx describe --json $analysis_id | jq-linux64 -e 'has("executable")' ; then
            workflow_id=$(dx describe --json $analysis_id | jq-linux64 -r '.executable')
            analysis_name=$(dx describe --json $workflow_id | jq-linux64 -r '.name')
            found_workflow_id=true
        fi
    fi

    # get read stats from flagstat file
    total_nb_reads=$(grep total inputs/$flagstat_file_name | cut -d+ -f1)
    nb_duplicates_reads=$(grep duplicate inputs/$flagstat_file_name | cut -d+ -f1)
    nb_aligned_reads=$(grep "mapped (" inputs/$flagstat_file_name | cut -d+ -f1)
    nb_usable_reads=$(expr $nb_aligned_reads - $nb_duplicates_reads)

    cd packages

    tar xjf htslib-1.7.tar.bz2
    cd htslib-1.7
    sudo ./configure --prefix=/usr/bin
    sudo make
    sudo make install
    cd ..

    tar xjf samtools-1.7.tar.bz2
    cd samtools-1.7
    sudo ./configure --prefix=/usr/bin
    sudo make
    sudo make install
    cd ..

    # Compile perl packages
    for perl_package in $(ls *gz); do
        base_name=${perl_package%.*.*}
        echo "Installing: $base_name"
        tar xzf $perl_package
        cd $base_name
        perl Makefile.PL
        make
        sudo make install
        cd ..
    done

    opts="-a inputs/annotated_vcf "
    opts+="-s inputs/sliced_annotated_vcf "
    opts+="-v inputs/$raw_vcf_name "
    opts+="-c inputs/$sample_coverage_file_name "
    opts+="-u $nb_usable_reads "
    opts+="-T $total_nb_reads "
    opts+="-w \"$analysis_name\" "
    opts+="-i \"$workflow_id\""

    if [ ! -z ${list_panel_names_genes+x} ]; then
        opts+=" -p \"$list_panel_names_genes\""
    fi

    if [ ! -z ${annotations+x} ]; then
        opts+=" -f \"$parsed_annotations\""
    fi

    cd /home/dnanexus

    # Run vcf2xls
    eval "perl -I /home/dnanexus/ vcf2xls_nirvana.pl ${opts}"

    project_id=$DX_PROJECT_CONTEXT_ID

    version=0
    matching_files=1
    output_name=${sample_id}_${version}.xls

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
    
    cp /home/dnanexus/out/xls_reports/report.xls /home/dnanexus/out/xls_reports/${output_name}

    output_file=$(dx upload /home/dnanexus/out/xls_reports/${output_name} --brief)

    dx-jobutil-add-output xls_report "$output_file" --class=file
}
