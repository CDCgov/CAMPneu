process SRA_HUMAN_SCRUBBER {
    tag "$meta.id"
    label 'process_medium'

    //container "https://depot.galaxyproject.org/singularity/sra-human-scrubber:2.2.1--hdfd78af_0"
    //container "docker://ncbi/sra-human-scrubber:2.2.1"
    //container "quay.io/biocontainers/sra-human-scrubber:2.2.1--hdfd78af_0"
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0c/0c06633f1ce36c4cfb4d9e0f59cd15ce38702decda24c425ca8ed0f8f516aa29/data' :
    //     'community.wave.seqera.io/library/sra-human-scrubber:2.2.1--147fe7682c448f16' }"
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