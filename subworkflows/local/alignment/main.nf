//
// Alignment with reference and stats
//
include { MINIMAP2_ALIGN                             } from '../../../modules/nf-core/minimap2/align/main'
include { SAMTOOLS_STATS                             } from '../../../modules/local/samtools/stats/main'
include { SAMTOOLS_STATS as SAMTOOLS_STATS10X        } from '../../../modules/local/samtools/stats/main'
include { SAMTOOLS_DEPTH                             } from '../../../modules/nf-core/samtools/depth/main'
include { REPORT_DEPTH_COV                           } from '../../../modules/local/coverage_depth_report/main'

workflow ALIGNMENT {
    take:
    reads          // channel: [ id: meta, [ reads ] ]
    report_name    // val: "*.tsv"
    ch_faidx       // channel: [ meta, index ]

    main:

    ch_versions = Channel.empty()

    //
    // MODULE: Align to reference using Minimap2
    //
    MINIMAP2_ALIGN (
        reads,
        ch_faidx,
        true,
        "bai",
        false,
        false
    )
    ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions.first())

    ch_bam = MINIMAP2_ALIGN.out.bam
    ch_bam_bai = MINIMAP2_ALIGN.out.bam
        .join(MINIMAP2_ALIGN.out.index)
    

    //
    // MODULE: Run Samtools stats
    //
    SAMTOOLS_STATS (
        ch_bam_bai,
        ch_faidx.map { meta, ref -> tuple(meta, ref, "${params.ref_name}") },
        "1X"
    )
    ch_versions = ch_versions.mix(SAMTOOLS_STATS.out.versions)

    //
    // MODULE: Run Samtools stats to get breadth of coverage with 10X depth
    //
    SAMTOOLS_STATS10X (
        ch_bam_bai,
        ch_faidx.map { meta, ref -> tuple(meta, ref, "${params.ref_name}") },
        "10X"
    )
    ch_versions = ch_versions.mix(SAMTOOLS_STATS10X.out.versions)

    //
    // MODULE: Run Samtools depth
    //
    SAMTOOLS_DEPTH (
        ch_bam,
        [[],[]]
    )
    ch_versions = ch_versions.mix(SAMTOOLS_DEPTH.out.versions)

    ch_stats_depth = SAMTOOLS_STATS.out.stats
                .join(SAMTOOLS_STATS10X.out.stats)
                .join(SAMTOOLS_DEPTH.out.tsv)
    //
    // MODULE: Parse Samtools stats & depth output, get breadth & depth of cov
    //
    REPORT_DEPTH_COV (
        ch_stats_depth
    )
    ch_versions = ch_versions.mix(REPORT_DEPTH_COV.out.versions)

    ch_stats = REPORT_DEPTH_COV.out.tsv

    // Merge depth cov reports
    ch_merge = REPORT_DEPTH_COV.out.tsv
                    .collectFile(name:"${report_name}", storeDir:"${params.outdir}/reports/", keepHeader:true){
                        meta, file -> file
                    }


    emit:

    stats               = ch_merge                                      // channel: [ report_name.tsv ]

    ch_stats            = REPORT_DEPTH_COV.out.tsv                      // channel: [ meta, tsv ]


    bam_bai             = ch_bam_bai                                    // channel: [ meta, bam, bai ]


    versions            = ch_versions                                   // channel: [ versions.yml ]

}