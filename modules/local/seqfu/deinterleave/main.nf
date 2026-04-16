process SEQFU_DEINTERLEAVE {
    tag "$meta.id"
    label 'process_low'
    conda "${moduleDir}/environment.yml"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqfu:1.20.3--h1eb128b_2':
        'biocontainers/seqfu:1.20.3--h1eb128b_2' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('*.fq.gz'), emit: reads
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    seqfu deinterleave \\
        -o ${prefix} \\
        $args \\
        $reads
    
    gzip *.fq
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqfu: \$(seqfu version)
    END_VERSIONS
    """
}