//
// Generate synthetic reads, perform alignment with reference
//
include { QUAST                                      } from '../../../modules/local/quast/main'
include { BBMAP_RANDOMREADS                          } from '../../../modules/local/bbmap/randomreads/main'
include { MINIMAP2_ALIGN                             } from '../../../modules/nf-core/minimap2/align/main'

workflow ASSEMBLYALIGNMENT {
    take:
    ch_contigs     // channel: [ id: meta, fasta ]
    ch_faidx       // channel: [ meta, index ]

    main:

    ch_versions = Channel.empty()

    //
    // MODULE: Generate assembly quality metrics 
    //
    ch_quast = ch_contigs
                .map {
                    meta, contigs ->
                    [ meta, contigs, [] ]
                }
    QUAST(
        ch_quast,
        ['ref',"${params.reference_genome}"],
        ['gff',"${params.ref_annotation}"]
    )
    ch_versions = ch_versions.mix(QUAST.out.versions)

    // Perform QC filtering on input assemblies
    ch_assembly_qc_metrics = QUAST.out.transposed
                                .map {
                                    meta, quast ->
                                    def q_parsed = quast.splitCsv(sep:'\t', header:true)
                                    //total length in KB
                                    def length = (q_parsed['Total length (>= 0 bp)'][0].toFloat()/1000).round(2)
                                    def c_count = q_parsed['# contigs'][0].toInteger()
                                    //N50 in KB
                                    def n50 = (q_parsed['N50'][0].toFloat()/1000).round(2)
                                    def gc = q_parsed['GC (%)'][0].toFloat()
                                    return [ meta, [l:length, cc:c_count, n50:n50, gc:gc ]]
                                }
    
    ch_assembly_qc = ch_assembly_qc_metrics
                        .branch {
                            meta, s ->
                                passed: ( 700 <= s.l && s.l <= 900 ) && s.cc <= 100 && s.n50 >= 49 && ( 39 <= s.gc && s.gc <= 41 )
                                failed: ( s.l < 700 || 900 < s.l ) || s.cc > 100 || s.n50 < 49 || ( s.gc < 39 || 41 < s.gc )
                        }
    
    ch_contigs_passed = ch_contigs.join(ch_assembly_qc.passed)
                            .map {
                                meta, contigs, stats ->
                                [ meta, contigs ]
                            }

    //
    // MODULE: Generate synthetic reads with randomreads
    //
    BBMAP_RANDOMREADS (
        ch_contigs_passed
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
    passed_contigs      = ch_contigs_passed                             // channel: [ meta, contigs ]
    assembly_metrics    = ch_assembly_qc_metrics                        // channel: [ meta, metrics ]

    versions            = ch_versions                                   // channel: [ versions.yml ]

}