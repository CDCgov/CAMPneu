//
// Unicycler assembly & assembly QC
//

include { UNICYCLER              } from '../../../modules/nf-core/unicycler/main'
include { QUAST                  } from '../../../modules/local/quast/main'

workflow ASSEMBLY {
    take:
    reads             // channel: [ val(sample_name), [ reads ] ]

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


    ch_quast = ch_assembly.join(reads)
    //
    // MODULE: Generate assembly quality metrics 
    //
    QUAST(
        ch_quast,
        ['ref',"${params.reference_genome}"],
        ['gff',"${params.ref_annotation}"]
    )
    ch_versions = ch_versions.mix(QUAST.out.versions)

    emit:

    contigs             = ch_assembly                                // channel:  [ id: meta, [ contigs ] ] 

    quast_results       = QUAST.out.results                          // channel: [ id: meta, [ results ] ]

    versions            = ch_versions                                // channel: [ versions.yml ]

}