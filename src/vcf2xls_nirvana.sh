#!/bin/bash
# vcf2xls_nirvana 2.0.0

set -euxo pipefail

main() {
    if [ ! -z ${list_panel_names_genes+x} ]; then
        echo "Value of list_panel_names_genes: '$list_panel_names_genes'"
    fi
    echo "Value of annotated_vcfs: '$annotated_vcfs'"
    echo "Value of raw_vcfs: '$raw_vcfs'"
    echo "Value of runfolder_coverage_file: '$runfolder_coverage_file'"
    echo "Value of coverage_files: '$coverage_files'"
    echo "Value of runfolder_coverage_index: '$runfolder_coverage_index'"
    echo "Value of coverage_indexes: '$coverage_indexes'"

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

    for i in ${!annotated_vcfs[@]}; do
        dx download "${annotated_vcfs[$i]}" -o inputs/.
    done

    for i in ${!raw_vcfs[@]}; do
        dx download "${raw_vcfs[$i]}" -o inputs/.
    done

    for i in ${!coverage_files[@]}; do
        dx download "${coverage_files[$i]}" -o inputs/.
    done

    for i in ${!coverage_indexes[@]}; do
        dx download "${coverage_indexes[$i]}" -o inputs/.
    done

    cd inputs

    # Get unique sample names for use in the perl script
    basenames=()

    for file in $(ls); do
        base=${file%%.*}

        # bash < 4.4 doesn't like empty arrays with set -u
        # thus the ugly thing
        if [[ ! ${basenames[@]+"${basenames[@]}"} =~ ${base} ]] && [[ ${base} =~ ^[XG] ]]; then
            basenames+=( "${base}" )
        fi
    done

    cd ..

    # Download runfolder coverage file and assign to variable
    dx download "$runfolder_coverage_file" -o inputs/.
    dx download "$runfolder_coverage_index" -o inputs/.

    runfolder_coverage_file=$(find inputs/ | grep -E '*refseq_nirvana_5bp.gz$')

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

    # Run vcf2xls for all samples
    for sample in ${basenames[@]}; do
        annotated_vcf=inputs/${sample}.refseq_nirvana_203.annotated.vcf
        raw_vcf=inputs/${sample}.refseq_nirvana_203.vcf
        coverage_file=inputs/${sample}.nirvana_203_5bp.gz

        if [ -z ${list_panel_names_genes+x} ]; then
            echo "Running: perl vcf2xls_nirvana.pl -a $annotated_vcf -v $raw_vcf -R $runfolder_coverage_file -C $coverage_file" 
            perl vcf2xls_nirvana.pl -a $annotated_vcf -v $raw_vcf -R $runfolder_coverage_file -C $coverage_file
        else
            echo "Running: perl vcf2xls_nirvana.pl -p $list_panel_names_genes -a $annotated_vcf -v $raw_vcf -R $runfolder_coverage_file -C $coverage_file" 
            perl vcf2xls_nirvana.pl -p $list_panel_names_genes -a $annotated_vcf -v $raw_vcf -R $runfolder_coverage_file -C $coverage_file
        fi
    done

    dx-upload-all-outputs
}
