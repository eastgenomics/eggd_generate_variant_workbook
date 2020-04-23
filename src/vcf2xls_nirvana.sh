#!/bin/bash
# vcf2xls_nirvana 2.0.0

set -euxo pipefail

main() {
    if [ ! -z ${list_panel_names_genes+x} ]; then
        echo "Value of list_panel_names_genes: '$list_panel_names_genes'"
    fi
    echo "Value of raw_vcf: '$raw_vcf'"
    echo "Value of annotated_vcf: '$annotated_vcf'"
    echo "Value of runfolder_coverage_file: '$runfolder_coverage_file'"
    echo "Value of runfolder_coverage_index: '$runfolder_coverage_index'"
    echo "Value of sample_coverage_file: '$sample_coverage_file'"
    echo "Value of sample_coverage_index: '$sample_coverage_index'"
    echo "Value of flagstat_file: '$flagstat_file'"

    # Compile samtools and tabix for perl script
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

    cd /home/dnanexus

    # Download sample inputs
    mkdir inputs

    sample_id=$(echo $annotated_vcf_prefix | awk -F "_" '{print $1}')
    echo $sample_id

    job_id=$(dx describe --delimiter "_" $annotated_vcf_name | grep job- | cut -d_ -f2)
    analysis_id=$(dx describe --delimiter "_" $job_id | grep Root | cut -d_ -f2)
    workflow_id=$(dx describe --delim "_" $analysis_id | grep Workflow | cut -d_ -f2)
    analysis_name=$(dx describe --name $analysis_id)

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

    total_nb_reads=$(grep total inputs/$flagstat_file_name | cut -d+ -f1)
    nb_duplicates_reads=$(grep duplicate inputs/$flagstat_file_name | cut -d+ -f1)
    nb_aligned_reads=$(grep "mapped (" inputs/$flagstat_file_name | cut -d+ -f1)
    nb_usable_reads=$(expr $nb_aligned_reads - $nb_duplicates_reads)

    # Download runfolder coverage file and assign to variable
    dx download "$runfolder_coverage_file" -o inputs/
    echo $runfolder_coverage_file_name
    dx download "$runfolder_coverage_index" -o inputs/
    echo $runfolder_coverage_index_name

    # Download reference files
    dx download "project-Fjj60Qj4yBGvQXbb5Z6FXkgF:file-Fk0X8704yBGYGJYp09xBqkK0" -o exons_nirvana203
    dx download "project-Fjj60Qj4yBGvQXbb5Z6FXkgF:file-FkPp44j4yBGbVbPz1xFZ59PB" -o BioinformaticManifest
    dx download "project-Fjj60Qj4yBGvQXbb5Z6FXkgF:file-Fk0X2Q04yBGZvF276zYpbfK1" -o nirvana_genes2transcripts
    dx download "project-Fjj60Qj4yBGvQXbb5Z6FXkgF:file-Fk0X7zj4yBGbVGkjK7j7Pkg5" -o genepanels
    dx download "project-Fjj60Qj4yBGvQXbb5Z6FXkgF:file-Fk0X7Q04yBGVP2fx7967FgvK" -o gemini_freq.vcf.gz
    dx download "project-Fjj60Qj4yBGvQXbb5Z6FXkgF:file-Fk0X2p84yBGxbvYx5ZkQX0Kz" -o esp_vcf.tab.gz
    dx download "project-Fjj60Qj4yBGvQXbb5Z6FXkgF:file-Fk0X2Xj4yBGVQg5BJ1QjVBKY" -o kg_vcf.tab.gz
    dx download "project-Fjj60Qj4yBGvQXbb5Z6FXkgF:file-Fk0X2y04yBGb2Jf851vvY3Bg" -o exac_vcf.sites.vep.vcf.gz
    dx download "project-Fjj60Qj4yBGvQXbb5Z6FXkgF:file-Fk18f1j4yBGk4P8K6vPFgKYV" -o gemini_freq.vcf.gz.tbi
    dx download "project-Fjj60Qj4yBGvQXbb5Z6FXkgF:file-Fk18f6Q4yBGbJ2BPJ5vGZ62x" -o esp_vcf.tab.gz.tbi
    dx download "project-Fjj60Qj4yBGvQXbb5Z6FXkgF:file-Fk18by04yBGv3KfKPp0v19V7" -o kg_vcf.tab.gz.tbi
    dx download "project-Fjj60Qj4yBGvQXbb5Z6FXkgF:file-Fk18bkQ4yBGb1VfQGxxVpjbY" -o exac_vcf.sites.vep.vcf.gz.tbi

    mkdir -p /home/dnanexus/out/xls_reports

    # Run vcf2xls
 
    if [ -z ${list_panel_names_genes+x} ]; then
        echo "Running: perl vcf2xls_nirvana.pl -a inputs/annotated.vcf -v inputs/raw.vcf -R inputs/runfolder_coverage.gz -C inputs/sample_coverage.gz" 
        perl vcf2xls_nirvana.pl \
            -a inputs/$annotated_vcf_name \
            -v inputs/$raw_vcf_name \
            -R inputs/$runfolder_coverage_file_name \
            -C inputs/$sample_coverage_file_name \
            -u $nb_usable_reads \
            -T $total_nb_reads \
            -w $analysis_name \
            -i $workflow_id
    else
        echo "Running: perl vcf2xls_nirvana.pl -p $list_panel_names_genes -a inputs/annotated.vcf -v inputs/raw.vcf -R inputs/runfolder_coverage.gz -C inputs/sample_coverage.gz" 
        perl vcf2xls_nirvana.pl \
            -p $list_panel_names_genes \
            -a inputs/$annotated_vcf_name \
            -v inputs/$raw_vcf_name \
            -R inputs/$runfolder_coverage_file_name \
            -C inputs/$sample_coverage_file_name \
            -u $nb_usable_reads \
            -T $total_nb_reads \
            -w $analysis_name \
            -i $workflow_id
    fi

    dx select
    0
    source ~/.dnanexus_config/unsetenv

    output_name=$sample_id.xls
    version=1

    while $(dx find data | grep -q $output_name); do
        version=$((version+1))
        output_name=${sample_id}_${version}.xls
    done

    cp /home/dnanexus/out/xls_reports/report.xls /home/dnanexus/out/xls_reports/${output_name}

    output_file=$(dx upload /home/dnanexus/out/xls_reports/${output_name} --brief)

    dx-jobutil-add-output xls_report "$output_file" --class=file
}
