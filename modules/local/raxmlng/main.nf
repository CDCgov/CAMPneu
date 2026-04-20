process RAXMLNG {
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/raxml-ng:1.2.2--h6747034_1' :
        'biocontainers/raxml-ng:1.2.2--h6747034_1' }"

    input:
    tuple val(meta), path(alignment)
    val model

    output:
    tuple val(meta), path("*.bestModel")         , emit: bestModel
    tuple val(meta), path("*.bestTree")          , emit: bestTree
    tuple val(meta), path("*.bestTreeCollapsed") , emit: bestTreeCollapsed
    tuple val(meta), path("*.raxml.log")         , emit: log
    tuple val(meta), path("*.mlTrees")           , emit: mlTrees
    tuple val(meta), path("*.rba")               , emit: rba
    tuple val(meta), path("*.startTree")         , emit: startTree
    tuple val(meta), path("*.support")           , emit: support
    tuple val(meta), path("*.bootstraps")        , emit: bootstraps
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // fix random seed for reproducibility if not specified in command line
    if (!(args ==~ /.*--seed.*/)) {args += " --seed=42"}
    """
    raxml-ng \\
        $args \\
        --msa $alignment \\
        --model $model \\
        --threads $task.cpus \\
        --prefix ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        raxmlng: \$(echo \$(raxml-ng --version 2>&1) | sed 's/^.*RAxML-NG v. //; s/released.*\$//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def touch_files = args.contains('--bootstrap') || args.contains('--bs-trees') ? "touch ${prefix}.raxml.bootstraps" : "touch ${prefix}.raxml.bestTree"
    """
    # Create stub output files
    ${touch_files}

    # Create versions.yml
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        raxmlng: \$(echo \$(raxml-ng --version 2>&1) | sed 's/^.*RAxML-NG v. //; s/released.*\$//')
    END_VERSIONS
    """
}
