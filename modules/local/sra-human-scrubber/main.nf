process SRA_HUMAN_SCRUBBER {
    tag "$meta.id"
    label 'process_medium'

    container "library://kmorin/campneu/sra-human-scrubber:2.2.1"

    input:
    tuple val(meta), path(reads)
    path(db)

    output: 
    tuple val(meta), path('*clean.fastq.gz'), emit: reads
    tuple val(meta), path('*.removed_spots'), optional:true, emit: removed
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
 
    zcat $reads | scrub.sh \\
        -d $db \\
        -o ${prefix}_clean.fastq \\
        -p $task.cpus \\
        $args
    
    gzip ${prefix}_clean.fastq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sra-human-scrubber: \$(echo 2.2.1)
    END_VERSIONS
    """
}