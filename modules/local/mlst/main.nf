process MLST {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mlst:2.23.0--hdfd78af_0' :
        'biocontainers/mlst:2.23.0--hdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)
    path blastdb
    path datadir
    val db_name

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args     = task.ext.args ?: ''
    def prefix   = task.ext.prefix ?: "${meta.id}"
    def blast_db = blastdb ? "--blastdb ${blastdb}/${db_name}" : ""
    def data_dir = datadir ? "--datadir ${datadir}" : ""
    """
    mlst \\
        $args \\
        --threads $task.cpus \\
        $blast_db \\
        $data_dir \\
        $fasta \\
        > ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mlst: \$( echo \$(mlst --version 2>&1) | sed 's/mlst //' )
    END_VERSIONS
    """

}
