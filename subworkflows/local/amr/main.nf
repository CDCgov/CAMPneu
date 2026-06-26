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

    // Merge AMRFinderPlus reports
    ch_amr_report = AMRFINDERPLUS_RUN.out.report
                        .map {
                            meta, tsv ->
                            def file = tsv
                                        .splitCsv( sep:"\t" )
                            return [ meta, [ symbol:file[6], name:file[7], type:file[9], class:file[11], cov:file[16], ident:file[17], acc:file[19] ] ]
                        }
                        .collectFile(name:"AMRFinderPlus_report.tsv", seed: 'Sample\tElement_symbol\tElement_name\tType\tClass\t%_Coverage_of_reference\t%_Identity_to_reference\tClosest_reference_accession',storeDir:"${params.outdir}/reports/", cache:false, newLine:true){
                            meta, amr ->
                            [ 'AMRFinderPlus_report.tsv', meta.id + '\t' + amr.symbol + '\t'+ amr.name + '\t' + amr.type + '\t' + amr.class + '\t' + amr.cov + '\t' + amr.ident + '\t' + amr.acc ]
                        }
                        // .collectFile(name:'AMRFinderPlus_report.tsv', storeDir:"${params.outdir}/reports/", keepHeader:true){
                        //     meta, file -> file
                        // }

    emit:
    snp_report          = ch_snp_report                              // channel: [ SNP_report.tsv ]
    amr_report          = ch_amr_report                              // channel: [ AMRFinderPlus_report.tsv ]
    versions            = ch_versions                                // channel: [ versions.yml ]

}
