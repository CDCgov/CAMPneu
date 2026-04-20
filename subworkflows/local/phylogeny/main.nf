//
// Phylogenetic analysis and tree generation
//

include { SNIPPY               } from '../../../modules/local/snippy/main'
include { SNIPPY_CORE          } from '../../../modules/local/snippy_core/main'
include { SNIPPY_CLEAN         } from '../../../modules/local/snippy_clean/main'
include { REMOVE_REF           } from '../../../modules/local/remove_ref/main'
include { SNPDISTS             } from '../../../modules/nf-core/snpdists/main'
include { GUBBINS              } from '../../../modules/nf-core/gubbins/main'
include { RAXMLNG              } from '../../../modules/local/raxmlng/main'

workflow PHYLOGENY {
    take:
    ch_contigs
    //samples        // channel: [ sample ]

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    //
    // Run Snippy with input reference
    //
    SNIPPY (
        ch_contigs,
        "${params.reference_genome}"
    )
    ch_multiqc_files = ch_multiqc_files.mix(SNIPPY.out.txt.collect{it[1]})
    ch_versions = ch_versions.mix(SNIPPY.out.versions)

    ch_vcf = SNIPPY.out.vcf
                .map {
                    meta, vcf ->
                    vcf
                }
                .collect()
    ch_aligned_fa = SNIPPY.out.aligned_fa
                        .map {
                            meta, aligned_fa ->
                            aligned_fa
                        }
                        .collect()


    //
    // Align core SNPs
    //
    SNIPPY_CORE (
        ch_vcf,
        ch_aligned_fa,
        "${params.reference_genome}"
    )
    ch_versions = ch_versions.mix(SNIPPY_CORE.out.versions)

    //
    // Cleanup core SNP alignment 
    //
    ch_clean = Channel.empty()
    SNIPPY_CLEAN (
        SNIPPY_CORE.out.full_aln
    )
    ch_clean = SNIPPY_CLEAN.out.clean_full_aln
    ch_versions = ch_versions.mix(SNIPPY_CLEAN.out.versions)

    //
    // Remove reference from Snippy core alignment
    //
    REMOVE_REF (
        ch_clean
    )
    ch_clean = REMOVE_REF.out.no_ref_aln
                .map {
                        aln ->
                        [ [ id:"core_aln" ], aln ]
                    }
    ch_versions = ch_versions.mix(REMOVE_REF.out.versions)

    //
    // MODULE: Compute SNP distances
    //
    SNPDISTS (
        ch_clean
    )
    ch_versions = ch_versions.mix(SNPDISTS.out.versions)

    if (params.gubbins){
        //
        // Mark recombination regions and contruct phylogeny on mutations outside of recombination regions
        //
        GUBBINS (
            ch_clean
        )
        ch_versions = ch_versions.mix(GUBBINS.out.versions)
        ch_clean = GUBBINS.out.phylip
    }

    //
    // Phylogenetic analysis with RAxML-NG (output Newick tree)
    //
    RAXMLNG (
        ch_clean,
        "GTR+G"
    )
    ch_versions = ch_versions.mix(RAXMLNG.out.versions)

    emit:

    versions             = ch_versions                             // channel: [ versions.yml ]

    multiqc              = ch_multiqc_files
}