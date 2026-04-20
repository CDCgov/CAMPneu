process SEQKIT_SAMPLE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqkit:2.10.0--h9ee0642_0':
        'biocontainers/seqkit:2.10.0--h9ee0642_0' }"

    input:
    tuple val(meta), path(fastqs), val(proportion)

    output:
    tuple val(meta), path("*_downsampled_*.fastq.gz"), emit: downsampled
    path "versions.yml"                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    """
    seqkit sample \\
        -p "${proportion}" \\
        -s 42 \\
        ${fastqs[0]} > ${prefix}_downsampled_R1.fastq
    
    seqkit sample \\
        -p "${proportion}" \\
        -s 42 \\
        ${fastqs[1]} > ${prefix}_downsampled_R2.fastq
    
    gzip *.fastq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$( seqkit version | sed 's/seqkit v//' )
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    for f in ${fastqs.join(' ')}
    do
       echo "" | gzip > "${prefix}_subset_\$(basename \$f)"
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$( seqkit version | sed 's/seqkit v//' )
    END_VERSIONS
    """
}