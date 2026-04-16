//
// Check for AMR markers
//
include { SAMTOOLS_FAIDX                       } from '../../../modules/nf-core/samtools/faidx/main'
include { FREEBAYES                            } from '../../../modules/local/freebayes/main'
include { SNP_SUMMARY                          } from '../../../modules/local/snp_summary/main'
include { AMRFINDERPLUS_RUN                    } from '../../../modules/local/amrfinderplus/run/main'

workflow AMR {
    take:
    ch_bam_bai        // channel: [ meta, bam, bai ]
    ch_contigs        // channel: [ meta, contigs ]

    main:

    ch_versions = Channel.empty()
    
    //
    // MODULE: Index reference genome with faidx
    //
    ch_ref = Channel.of([[id: "reference"],"${params.reference_genome}"]).collect()
    SAMTOOLS_FAIDX (
        ch_ref,
        [[],[]],
        false
    )
    //ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)
    

    ch_freebayes = ch_bam_bai
                    .map {
                        meta, bam, bai ->
                        [ meta, bam, bai, [], [], "${params.target_bed}" ]
                    }
    //
    // MODULE: Variant calling with FreeBayes
    //
    FREEBAYES (
        ch_freebayes,
        ch_ref.collect(),
        SAMTOOLS_FAIDX.out.fai.collect(),
        [[],[]],
        [[],[]],
        [[],[]]
    )
    ch_versions = ch_versions.mix(FREEBAYES.out.versions)

    //
    // MODULE: Get AMR SNP summary 
    //
    SNP_SUMMARY(
        FREEBAYES.out.vcf
    )
    ch_versions = ch_versions.mix(SNP_SUMMARY.out.versions)

    // Merge SNP reports
    ch_snp_report = SNP_SUMMARY.out.tsv
                        .collectFile(name:'SNP_report.tsv', storeDir:"${params.outdir}/reports/", keepHeader:true){
                            meta, file -> file
                        }

    //
    // MODULE: AMR Gene indentification
    //
    AMRFINDERPLUS_RUN (
        ch_contigs,
        "${params.amrfinderplus_db}"
    )
    ch_versions = ch_versions.mix(AMRFINDERPLUS_RUN.out.versions)

    emit:
    snp_report          = ch_snp_report                              // channel: [ 'SNP_report.tsv' ]
    versions            = ch_versions                                // channel: [ versions.yml ]

}
