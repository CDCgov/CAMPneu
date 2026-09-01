process PARSEKRAKEN2 {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqkit:2.10.0--h9ee0642_0':
        'biocontainers/seqkit:2.10.0--h9ee0642_0' }"

    input:
    tuple val(meta), path(reads), path(tax_class)

    output:
    tuple val(meta), path('*.filtered.fastq.gz'), emit: reads
    tuple val("${task.process}"), val('seqkit'), eval("seqkit version | sed 's/seqkit v//'"), topic: versions, emit: versions_seqkit

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """

    grep -E 'Mycoplasmoid(es|aceae|ales)' "$tax_class" | awk '{ print \$2 }' > read_ids.txt

    seqkit grep -f read_ids.txt ${reads[0]} -o ${prefix}_1.filtered.fastq.gz
    seqkit grep -f read_ids.txt ${reads[1]} -o ${prefix}_2.filtered.fastq.gz

    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    
    touch ${prefix}_1.filtered.fastq
    gzip ${prefix}_1.filtered.fastq

    touch ${prefix}_2.filtered.fastq
    gzip ${prefix}_2.filtered.fastq

    """
}
