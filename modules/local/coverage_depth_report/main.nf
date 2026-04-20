process REPORT_DEPTH_COV {
    tag "$meta.id"
    label 'process_single'

    container "https://depot.galaxyproject.org/singularity/centos:7.9.2009"

    input: 
    tuple val(meta), path(stats1x), path(stats10x), path(depth)

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    path  "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """

    depth_cov_report.sh $stats1x $stats10x $depth $prefix

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Awk: \$(awk --version 2>&1 | sed -n 1p | sed 's/GNU Awk //')
    END_VERSIONS
    """
}