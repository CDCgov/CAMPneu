

process BBMAP_RANDOMREADS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bbmap:39.76--h9b5c0a0_0':
        'biocontainers/bbmap:39.76--h9b5c0a0_0' }"

    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path("*.fastq"), emit: reads
    tuple val("${task.process}"), val('bbmap'), eval("bbversion.sh"), topic: versions, emit: versions_bbmap

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    randomreads.sh \\
        ref=${input} \\
        out1=${prefix}_1.fastq \\
        out2=${prefix}_2.fastq \\
        ${args}
    
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    echo $args
    
    """
}
