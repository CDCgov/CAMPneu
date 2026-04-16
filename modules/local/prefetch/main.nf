process PREFETCH{
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? 
        'https://depot.galaxyproject.org/singularity/sra-tools:3.1.1--h4304569_2' : 
        'biocontainers/sra-tools:3.1.1--h4304569_2' }"

    input:
    tuple val(meta), val(sra_id)

    output:
    tuple val(meta), path("${sra_id[0]}/"), emit: prefetch
    path('versions.yml'), emit: versions

    script:
    """
    prefetch --verify yes ${sra_id[0]}

        cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sratools: \$(prefetch --version 2>&1 | sed 's/prefetch : //' | awk 'NF')
    END_VERSIONS
    """
}
