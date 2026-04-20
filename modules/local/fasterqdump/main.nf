process FASTERQDUMP {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? 
        'https://depot.galaxyproject.org/singularity/sra-tools:3.1.1--h4304569_2' : 
        'biocontainers/sra-tools:3.1.1--h4304569_2' }"

    input:
    tuple val(meta), file(prefetch)

    output:
    tuple val(meta), path("*.fastq.gz") , emit: fastq
    path('versions.yml'), emit: versions

    script:
    """
    fasterq-dump \\
    --threads $task.cpus \\
    --split-files \\
    --skip-technical \\
    --outdir . \\
    ${prefetch}

    gzip *.fastq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sratools: \$(fasterq-dump --version 2>&1 | sed 's/fasterq-dump : //' | awk 'NF')
    END_VERSIONS
    """

}