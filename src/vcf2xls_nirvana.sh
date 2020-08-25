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

    # get sample id from vcf file name
    sample_id=$(echo $annotated_vcf_prefix | awk -F "_" '{print $1}')
    echo $sample_id

    # Boolean to detect if workflow id has been found
    found_workflow_id=false

    # Placeholder text if the workflow id is not found
    analysis_name="No workflow id found for this report."
    workflow_id="This report was probably generated for development purposes, do not use for clinical reporting"


    # Get workflow name and id
    if dx describe --delim "_" $annotated_vcf_name | grep -q job- ; then
        job_id=$(dx describe --delim "_" $annotated_vcf_name | grep job- | cut -d_ -f2)
        analysis=$(dx describe --delim "_" $job_id)

        if dx describe --delim "_" $job_id | grep -q Root ; then
            analysis_id=$(dx describe --delim "_" $job_id | grep Root | cut -d_ -f2)
            workflow=$(dx describe --delim "_" $analysis_id)

            if dx describe --delim "_" $analysis_id | grep -q Workflow ; then
                workflow_id=$(dx describe --delim "_" $analysis_id | grep Workflow | cut -d_ -f2)
                analysis_name=$(dx describe --name $analysis_id)
                found_workflow_id=true
            fi
        fi
    fi

    # get read stats from flagstat file
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
    dx download "project-Fkb6Gkj433GVVvj73J7x8KbV:file-FpQpGq8433GZXPXgB34B3kz0" -o exons_nirvana203
    dx download "project-Fkb6Gkj433GVVvj73J7x8KbV:file-FpQpV0j433GqJXGvJ30B8p2Y" -o gemini_freq.vcf.gz
    dx download "project-Fkb6Gkj433GVVvj73J7x8KbV:file-FpQpG6j433Gyp0kF6F9F69qq" -o esp_vcf.tab.gz
    dx download "project-Fkb6Gkj433GVVvj73J7x8KbV:file-FpQpFpQ433Gv30GBPqz29V0k" -o kg_vcf.tab.gz
    dx download "project-Fkb6Gkj433GVVvj73J7x8KbV:file-FpQpPkj433Gfpz8g0x2X64jQ" -o exac_vcf.sites.vep.vcf.gz
    dx download "project-Fkb6Gkj433GVVvj73J7x8KbV:file-FpQpJ5Q433Gb2V5y3fxx09p0" -o gemini_freq.vcf.gz.tbi
    dx download "project-Fkb6Gkj433GVVvj73J7x8KbV:file-FpQpGPj433GzPByY1Vpfz7bb" -o esp_vcf.tab.gz.tbi
    dx download "project-Fkb6Gkj433GVVvj73J7x8KbV:file-FpQpFyj433GQkkJzFzbFb48J" -o kg_vcf.tab.gz.tbi
    dx download "project-Fkb6Gkj433GVVvj73J7x8KbV:file-FpQpGgQ433GyvPj8Fq65F4qP" -o exac_vcf.sites.vep.vcf.gz.tbi

    # Download dynamic reference files
    dx download "$genepanels_file" -o genepanels
    dx download "$bioinformatic_manifest" -o BioinformaticManifest
    dx download "$nirvana_genes2transcripts" -o nirvana_genes2transcripts
    dx download "$exons_nirvana" -o exons_nirvana

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
            -w "$analysis_name" \
            -i "$workflow_id"
    else
        echo "Running: perl vcf2xls_nirvana.pl -p \"$list_panel_names_genes\" -a inputs/annotated.vcf -v inputs/raw.vcf -R inputs/runfolder_coverage.gz -C inputs/sample_coverage.gz" 
        perl vcf2xls_nirvana.pl \
            -p "$list_panel_names_genes" \
            -a inputs/$annotated_vcf_name \
            -v inputs/$raw_vcf_name \
            -R inputs/$runfolder_coverage_file_name \
            -C inputs/$sample_coverage_file_name \
            -u $nb_usable_reads \
            -T $total_nb_reads \
            -w "$analysis_name" \
            -i "$workflow_id"
    fi

    project_id=$DX_PROJECT_CONTEXT_ID
    #dx select $project_id
    #source ~/.dnanexus_config/unsetenv

    version=0
    matching_files=1
    output_name=${sample_id}_${version}.xls

    # Tiny chance of race conditions leading to two files with the same name here
    while [ $matching_files -ne 0 ]; do 
        version=$((version+1))
        output_name=${sample_id}_${version}.xls; 
        matching_files=$(dx find data --path ${project_id}:/ --name $output_name --brief | wc -l); 
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