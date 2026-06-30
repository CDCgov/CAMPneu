process SUMMARY_REPORT {
    label 'process_single'

    container "docker://roboxes/rhel8:4.3.14"

    input:
    path(fastp_report)
    path(stats)
    path(ds_stats)
    path(mp_percent)
    path(quast)
    path(mlst)
    path(ani)
    path(snps)
    path(amr)
    val(depth)
    val(description)
    val(version)
    val(date)

   output:
    path("summary_report.txt"), emit: report

    when:
    task.ext.when == null || task.ext.when

    script:
    def fastp = fastp_report      ?  true : false
    def kraken = mp_percent       ?  true : false
    def align_stats = stats       ?  true : false
    def downsampled = ds_stats    ?  true : false
    """

    touch summary_report.txt

    echo "CAMPneu - ${description}" >> summary_report.txt
    echo "Version: ${version}" >> summary_report.txt
    echo "Run date: ${date}\n" >> summary_report.txt

    if ${fastp}; then
        echo "Read stats before and after filtering with fastp\n" >> summary_report.txt
        column -t ${fastp_report} >> summary_report.txt
        echo "---------------------------------------------------------------------------------------------------------\n" >> summary_report.txt
    fi 
    if ${align_stats}; then
        echo "Coverage and depth without downsampling\n" >> summary_report.txt
        column -t ${stats} >> summary_report.txt
        echo "---------------------------------------------------------------------------------------------------------\n" >> summary_report.txt
    fi 
    if ${downsampled}; then
        echo "Coverage and depth after downsampling to ${depth}x\n" >> summary_report.txt
        column -t ${ds_stats} >> summary_report.txt
        echo "---------------------------------------------------------------------------------------------------------\n" >> summary_report.txt
    fi
    if ${kraken}; then
        echo "Mp percentage determined by Kraken2\n" >> summary_report.txt
        column -t ${mp_percent} >> summary_report.txt
        echo "---------------------------------------------------------------------------------------------------------\n" >> summary_report.txt
    fi
    echo "Assembly metrics reported by QUAST\n" >> summary_report.txt
    column -t ${quast} >> summary_report.txt
    echo "---------------------------------------------------------------------------------------------------------\n" >> summary_report.txt
    echo "Sequence typing using MLST\n" >> summary_report.txt
    column -t ${mlst} >> summary_report.txt
    echo "---------------------------------------------------------------------------------------------------------\n" >> summary_report.txt
    echo "P1 type determined by ANI\n" >> summary_report.txt
    column -t ${ani} >> summary_report.txt
    echo "---------------------------------------------------------------------------------------------------------\n" >> summary_report.txt
    echo "Identification of macrolide resistant SNPs using FreeBayes\n" >> summary_report.txt
    column -t ${snps} >> summary_report.txt
    echo "---------------------------------------------------------------------------------------------------------\n" >> summary_report.txt
    echo "Identification of AMR genes with AMRFinderPlus\n" >> summary_report.txt
    column -t ${amr} >> summary_report.txt
    echo "---------------------------------------------------------------------------------------------------------\n" >> summary_report.txt

    """
    stub:
    """

    touch summary_report.txt

    """

}