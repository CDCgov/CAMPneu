process FASTANI {
    tag "$meta.id"
    label 'process_medium'

    // conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastani:1.32--he1c1bb9_0' :
        'biocontainers/fastani:1.32--he1c1bb9_0' }"

    input:
    tuple val(meta), path(query)
    path reference
    path type1  //bypass rawusercontent issues
    path type2  //bypass rawusercontent issues
    val querylist
    val referencelist

    output:
    tuple val(meta), path("*.ani.out.matrix") , optional: true , emit: matrix
    tuple val(meta), path("*.ani.out")        , emit: ani
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def q      = querylist       ?  "--ql $query" : "-q $query"
    def r      = referencelist   ?  "--rl $reference" : "-r $reference"

    """
    fastANI \\
        $q \\
        $r \\
        -t $task.cpus \\
        $args \\
        -o ${prefix}.ani.out

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastani: \$(fastANI --version 2>&1 | sed 's/version//;')
    END_VERSIONS
    """
}
