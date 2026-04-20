process GET_MP_PERCENT {
    tag "$meta.id"
    label 'process_single'

    container "docker://roboxes/rhel8:4.3.14"

    input:
    tuple val(meta), path(report)

    output:
    tuple val(meta), path("*_Mp_percent.tsv"), emit: tsv
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """

    awk -vFS="\t" -vOFS="\t" -vsample=${prefix} 'BEGIN{ print "Sample","Percent_Mp" }{ gsub(/[ ]+/,"",\$6); if( \$6 == "Mycoplasmoidespneumoniae" ){ print sample,\$1 } }' $report > ${prefix}_Mp_percent.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
         Awk: \$(awk --version 2>&1 | sed -n 1p | awk -vFS="," '{ print \$1 }' | sed 's/GNU Awk //')
    END_VERSIONS

    """
}