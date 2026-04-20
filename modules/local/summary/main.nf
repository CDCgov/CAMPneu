process SUMMARY_REPORT {
    label 'process_single'

    container "docker://roboxes/rhel8:4.3.14"

    input:
    path(stats)
    path(ds_stats)
    path(mp_percent)
    path(mlst)
    path(ani)
    path(snps)
    val(depth)

    output:
    path("summary_report.out"), emit: report

    when:
    task.ext.when == null || task.ext.when

    script:
    def kraken = mp_percent       ?  true : false
    def align_stats = stats       ?  true : false
    def downsampled = ds_stats    ?  true : false
    """

    touch summary_report.out

    if ${align_stats}; then
        echo "Coverage and depth without deduping or downsampling\n" >> summary_report.out
        column -t ${stats} >> summary_report.out
        echo "---------------------------------------------------------------------------------------------------------\n" >> summary_report.out
    fi 
    if ${downsampled}; then
        echo "Coverage and depth after downsampling to ${depth}x\n" >> summary_report.out
        column -t ${ds_stats} >> summary_report.out
        echo "---------------------------------------------------------------------------------------------------------\n" >> summary_report.out
    fi
    if ${kraken}; then
        echo "Mp percentage determined by Kraken2\n" >> summary_report.out
        column -t ${mp_percent} >> summary_report.out
        echo "---------------------------------------------------------------------------------------------------------\n" >> summary_report.out
    fi
    echo "Sequence typing using MLST\n" >> summary_report.out
    column -t ${mlst} >> summary_report.out
    echo "---------------------------------------------------------------------------------------------------------\n" >> summary_report.out
    echo "P1 type determined by ANI\n" >> summary_report.out
    column -t ${ani} >> summary_report.out
    echo "---------------------------------------------------------------------------------------------------------\n" >> summary_report.out
    echo "Identification of Macrolide Resistant SNPs using Freebayes\n" >> summary_report.out
    column -t ${snps} >> summary_report.out
    echo "---------------------------------------------------------------------------------------------------------\n" >> summary_report.out

    """
    stub:
    def deduped = dedup_stats     ?  true : false
    """

    summary_report.sh $stats $dedup_stats $mp_percent $mlst $ani $snps $deduped

    """

}