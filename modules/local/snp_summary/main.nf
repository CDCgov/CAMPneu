process SNP_SUMMARY {
    tag "$meta.id"
    label 'process_single'

    container "docker://roboxes/rhel8:4.3.14"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*_snps.tsv") , emit: tsv
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """

    echo "" | awk -vOFS="\t" '{ print "Sample","Pos","REF","ALT","SNP","QUAL","DEPTH","ALT_COUNT","Macrolide_Susceptibility(Sensitive/Resistant)" }' > ${prefix}_snps.tsv

    awk 'BEGIN{ found=0 }{ if (\$0 !~ /^#/){ if( \$2 != 122120 || \$5 != "T" ){ print \$0 }else{ print "" }; found += 1 } }END{ if(found == 0){ print "" }}' $vcf | \\
        awk -vOFS="\t" -vsample=${prefix} '{ if (\$0 != ""){ new_col=\$2-120056; split(\$8,arr,";");split(arr[8],dp,"=");split(arr[6],ao,"=");print sample,\$2,\$4,\$5,\$4 new_col \$5,\$6,dp[2],ao[2],"Resistant" }else{ print sample,"NA","NA","NA","NA","NA","NA","NA","Sensitive"} }' >> ${prefix}_snps.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
         Awk: \$(awk --version 2>&1 | sed -n 1p | awk -vFS="," '{ print \$1 }' | sed 's/GNU Awk //')
    END_VERSIONS

    """
}