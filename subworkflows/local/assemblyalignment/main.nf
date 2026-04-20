//
// Generate synthetic reads, perform alignment with reference
//
include { BBMAP_RANDOMREADS                          } from '../../../modules/local/bbmap/randomreads/main'
include { MINIMAP2_ALIGN                             } from '../../../modules/nf-core/minimap2/align/main'

workflow ASSEMBLYALIGNMENT {
    take:
    ch_contigs     // channel: [ id: meta, fasta ]
    ch_faidx       // channel: [ meta, index ]

    main:

    ch_versions = Channel.empty()

    //
    // MODULE: Generate synthetic reads with randomreads
    //
    BBMAP_RANDOMREADS (
        ch_contigs
    )

    //
    // MODULE: Align to reference using Minimap2
    //
    MINIMAP2_ALIGN (
        BBMAP_RANDOMREADS.out.reads,
        ch_faidx,
        true,
        "bai",
        false,
        false
    )
    ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions.first())

    ch_bam_bai = MINIMAP2_ALIGN.out.bam
        .join(MINIMAP2_ALIGN.out.index)

    emit:

    bam_bai             = ch_bam_bai                                    // channel: [ meta, bam, bai ]


    versions            = ch_versions                                   // channel: [ versions.yml ]

}