process PARSEKRAKEN2 {
    tag "$meta.id"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04':
        'ubuntu:20.04' }"

    input:
    tuple val(meta), path(reads), path(tax_class)

    output:
    tuple val(meta), path('*.filtered.fastq.gz'), emit: reads
    tuple val("${task.process}"), val('parsekraken2'), eval("parsekraken2 --version"), topic: versions, emit: versions_parsekraken2

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    gunzip $reads
    grep -E 'Mycoplasmoid(es|aceae|ales)' "$tax_class" | awk '{ print 2 }' > read_ids.txt

    cat read_ids.txt | xargs -I {} grep -A 4 {} ${reads[0]} > ${prefix}_1.filtered.fastq
    cat read_ids.txt | xargs -I {} grep -A 4 {} ${reads[1]} > ${prefix}_2.filtered.fastq

    gzip ${prefix}_1.filtered.fastq
    gzip ${prefix}_2.filtered.fastq
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo $args
    
    touch ${prefix}_1.filtered.fastq
    gzip ${prefix}_1.filtered.fastq

    touch ${prefix}_2.filtered.fastq
    gzip ${prefix}_2.filtered.fastq

    """
}
