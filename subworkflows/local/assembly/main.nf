//
// Unicycler assembly & assembly QC
//

include { UNICYCLER              } from '../../../modules/nf-core/unicycler/main'
include { QUAST                  } from '../../../modules/local/quast/main'

workflow ASSEMBLY {
    take:
    reads             // channel: [ val(sample_name), [ reads ] ]
    ch_contigs

    main:

    ch_versions = Channel.empty()
    ch_assembly = Channel.empty()

    ch_reads = reads
                    .map {
                        meta, reads ->
                            [ meta, reads, []]
                    }

    //
    // MODULE: unicycler assembly
    //
    UNICYCLER (
        ch_reads
    )
    ch_assembly = UNICYCLER.out.scaffolds
    ch_versions = ch_versions.mix(UNICYCLER.out.versions)

    ch_quast = ch_assembly.mix(ch_contigs)
                .map {
                    meta, contigs ->
                    [ meta, contigs, [] ]
                }
    //
    // MODULE: Generate assembly quality metrics 
    //
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
                                passed: ( 700 <= s.l && s.l <= 900 ) && s.cc <= 100 && s.n50 >= 25 && ( 39 <= s.gc && s.gc <= 41 )
                                failed: ( s.l < 700 || 900 < s.l ) || s.cc > 100 || s.n50 < 25 || ( s.gc < 39 || 41 < s.gc )
                        }
    
    ch_contigs_passed = ch_assembly.join(ch_assembly_qc.passed)
                            .map {
                                meta, contigs, stats ->
                                [ meta, contigs ]
                            }

    emit:

    passed_contigs      = ch_contigs_passed                          // channel: [ meta, contigs ]
    assembly_metrics    = ch_assembly_qc_metrics                     // channel: [ meta, metrics ]

    quast_results       = QUAST.out.results                          // channel: [ meta, results ]

    versions            = ch_versions                                // channel: [ versions.yml ]

}