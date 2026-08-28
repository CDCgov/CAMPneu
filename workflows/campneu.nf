/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { PREFETCH                    } from '../modules/local/prefetch/main'
include { FASTERQDUMP                 } from '../modules/local/fasterqdump/main'
include { PREPROCESSING               } from '../subworkflows/local/preprocess'
include { PREPROCESSING as PREPROCESS } from '../subworkflows/local/preprocess'
include { KRAKEN2_KRAKEN2             } from '../modules/nf-core/kraken2/kraken2/main'
include { GET_MP_PERCENT              } from '../modules/local/mp_percent/main'
include { ASSEMBLY                    } from '../subworkflows/local/assembly'
include { TYPING                      } from '../subworkflows/local/typing'
include { MINIMAP2_INDEX              } from '../modules/nf-core/minimap2/index/main'
include { ASSEMBLYALIGNMENT           } from '../subworkflows/local/assemblyalignment'
include { MINIMAP2_ALIGN              } from '../modules/nf-core/minimap2/align/main'
include { AMR                         } from '../subworkflows/local/amr'
include { SUMMARY_REPORT              } from '../modules/local/summary/main'
include { PHYLOGENY                   } from '../subworkflows/local/phylogeny'
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap            } from 'plugin/nf-schema'
include { paramsSummaryMultiqc        } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML      } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText      } from '../subworkflows/local/utils_nfcore_campneu_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow CAMPNEU {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    main:

    ch_versions = channel.empty()
    ch_multiqc_files = channel.empty()

    ch_samplesheet
        .filter { meta, it ->
            meta.dtype == 'FASTQ'
        }
        .set { ch_samplesheet_fastq }

    ch_samplesheet
        .filter { meta, it ->
            meta.dtype == 'FASTA'
        }
        .set { ch_samplesheet_fasta }
    
    ch_samplesheet
        .filter { meta, it ->
            meta.dtype == 'MIX'
        }
        .map {
            meta, it ->
            [ meta, [ it[0][0], it[0][1] ], [ it[0][2] ]]
        }
        .set { ch_samplesheet_mix }
    
    ch_samplesheet
        .filter { meta, it ->
            meta.dtype == 'SRA'
        }
        .set { ch_samplesheet_sra } // e.g. [[id:SRR11725329], [SRR11725329]]
    
    ch_only_fastq = ch_samplesheet_mix
                        .map {
                            meta, reads, fasta ->
                            [ meta, reads ]
                        }

    //
    // MODULE: Create Minimap2 index of reference genome
    //
    ch_ref = Channel.of([[id: "reference"],"${params.reference_genome}"]).collect()
    MINIMAP2_INDEX (
        ch_ref
    )
    ch_versions = ch_versions.mix(MINIMAP2_INDEX.out.versions)

    if (!params.phylogeny) {
        //
        // MODULE: Prefetch SRA samples
        //
        PREFETCH (
            ch_samplesheet_sra
        )
        ch_versions = ch_versions.mix(PREFETCH.out.versions)

        //
        // MODULE: Download SRA files
        //
        FASTERQDUMP (
            PREFETCH.out.prefetch
        )
        ch_versions = ch_versions.mix(FASTERQDUMP.out.versions)
        
        ch_samplesheet_fastq = ch_samplesheet_fastq.mix(FASTERQDUMP.out.fastq)
        //
        // SUBWORKFLOW: Align samples with only assemblies to the reference
        //
        ASSEMBLYALIGNMENT (
            ch_samplesheet_fasta,
            MINIMAP2_INDEX.out.index
        )
        ch_versions = ch_versions.mix(ASSEMBLYALIGNMENT.out.versions)

        //
        // SUBWORKFLOW: Preprocess FASTQ files
        //
        PREPROCESSING (
            ch_samplesheet_fastq,
            MINIMAP2_INDEX.out.index
        )
        ch_versions = ch_versions.mix(PREPROCESSING.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(PREPROCESSING.out.json.collect{it[1]})

        PREPROCESS (
            ch_only_fastq,
            MINIMAP2_INDEX.out.index
        )
        ch_versions = ch_versions.mix(PREPROCESS.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(PREPROCESS.out.json.collect{it[1]})


        // Merge fastp reports
        ch_fastp_report = PREPROCESSING.out.fastp_reports.mix(PREPROCESS.out.fastp_reports)
                        .map {
                            meta, tsv ->
                            def file = tsv
                                        .splitCsv( header:true, sep:"\t" )
                            return [ meta, [ before_reads:file.before_filtering_total_reads[0], before_q30:file.before_filtering_q30_rate[0], after_reads:file.after_filtering_total_reads[0], after_q30:file.after_filtering_q30_rate[0] ] ]
                        }
                        .collectFile(name:"fastp_report.tsv", seed: 'sample\tbefore_filtering_total_reads\tbefore_filtering_q30_rate\tafter_filtering_total_reads\tafter_filtering_q30_rate',storeDir:"${params.outdir}/reports/", cache:false, newLine:true){
                            meta, fastp ->
                            [ 'fastp_report.tsv', meta.id + '\t' + fastp.before_reads + '\t'+ fastp.before_q30 + '\t' + fastp.after_reads + '\t' + fastp.after_q30 ]
                        }

        ch_reads = PREPROCESSING.out.reads.mix(PREPROCESS.out.reads)

        // //
        // // MODULE: Classify reads with Kraken2
        // //
        // KRAKEN2_KRAKEN2 (
        //     ch_reads,
        //     "${params.kraken2db}",
        //     false,
        //     false
        // )
        // ch_versions = ch_versions.mix(KRAKEN2_KRAKEN2.out.versions)
        // ch_multiqc_files = ch_multiqc_files.mix(KRAKEN2_KRAKEN2.out.report.collect{it[1]})

        // //
        // // MODULE: Determine Mp percent from Kraken2 read classification
        // //
        // GET_MP_PERCENT (
        //     KRAKEN2_KRAKEN2.out.report
        // )
        // ch_versions = ch_versions.mix(GET_MP_PERCENT.out.versions)

        // // check if Mp percent is >90% to pass/fail samples
        // ch_checkmp = GET_MP_PERCENT.out.tsv
        //                 .map {
        //                     meta, tsv ->
        //                     def file = tsv
        //                                 .splitCsv( header:true, sep:"\t" )
        //                     def percent = file.Percent_Mp[0] as Float
        //                     return [ meta, percent ]
        //                 }
        //                 .branch {
        //                     meta, percent ->
        //                     pass: percent >= 90
        //                     fail: percent < 90
        //                 }

        // //Merge Mp percent reports
        // ch_percent_mp = GET_MP_PERCENT.out.tsv
        //                         .collectFile(name:'Mp_report.tsv', storeDir:"${params.outdir}/reports/", keepHeader:true, cache:false){
        //                             meta, file -> file
        //                         }
        //                         .ifEmpty([])
        
        // ch_reads = ch_reads.join(ch_checkmp.pass)
        //                             .map {
        //                                 meta, reads, percent ->
        //                                 [ meta, reads ]
        //                             }
        //only filtering this channel by avg depth (checked in preprocess) & Mp percent, not checking if corresponding assemblies passed
        // ch_bam_bai = PREPROCESSING.out.bam_bai.mix(PREPROCESS.out.bam_bai).join(ch_checkmp.pass)
        //                             .map {
        //                                 meta, bam, bai, percent ->
        //                                 [ meta, bam, bai ]
        //                             }
        //                             .mix(ASSEMBLYALIGNMENT.out.bam_bai)

        // ch_reads_to_assemble = PREPROCESSING.out.reads.join(ch_checkmp.pass)
        //                             .map{
        //                                 meta, reads, percent ->
        //                                 [ meta, reads ]
        //                             }
        
        //assemblies that were input alongside their readset 
        ch_only_fasta = ch_samplesheet_mix.join(PREPROCESS.out.checkmp)
                            .map {
                                meta, reads, fasta, percent ->
                                [ meta, fasta ]
                            }
        //
        // SUBWORKFLOW: Assembly & assembly QC
        //
        ASSEMBLY (
            ch_reads_to_assemble,
            ch_only_fasta
        )
        ch_multiqc_files = ch_multiqc_files.mix(ASSEMBLY.out.quast_results.collect{it[1]})
        ch_versions = ch_versions.mix(ASSEMBLY.out.versions)

        //Merge assembly QC reports
        ch_assembly_qc = ASSEMBLYALIGNMENT.out.assembly_metrics.mix(ASSEMBLY.out.assembly_metrics)
                                .collectFile(name:'quast_report.tsv', seed: 'sample\ttotal_length(KB)\tcontig_count\tN50(KB)\tGC_content(%)', storeDir:"${params.outdir}/reports/", cache:false, newLine:true){
                                    meta, metrics -> 
                                    [ 'quast_report.tsv', meta.id + '\t' + metrics.l + '\t'+ metrics.cc + '\t' + metrics.n50 + '\t' + metrics.gc ]
                                }

        
        ch_assembly = ASSEMBLYALIGNMENT.out.passed_contigs.mix(ASSEMBLY.out.passed_contigs)

        //
        // SUBWORKFLOW: P1 typing, sequence typing
        //
        TYPING (
            ch_assembly
        )
        ch_versions = ch_versions.mix(TYPING.out.versions)

        //
        // SUBWORKFLOW: AMR characterization
        //
        AMR (
            ch_bam_bai,
            ch_assembly
        )
        ch_versions = ch_versions.mix(AMR.out.versions)

        ch_stats = PREPROCESSING.out.stats.mix(PREPROCESS.out.stats)
                        .ifEmpty([])
        ch_ds_stats = PREPROCESSING.out.ds_stats.mix(PREPROCESS.out.ds_stats)
                        .ifEmpty([])

        def run_date = workflow.start.format("yyyy-MM-dd")
        //
        // Report results from run
        //
        SUMMARY_REPORT(
            ch_fastp_report,
            ch_stats,
            ch_ds_stats,
            ch_percent_mp,
            ch_assembly_qc,
            TYPING.out.mlst_report,
            TYPING.out.ani_report,
            AMR.out.snp_report,
            AMR.out.amr_report,
            "${params.depth}",
            "${workflow.manifest.description}",
            "${workflow.manifest.version}",
            "${run_date}"
        )
    } 

    if (params.phylogeny){ 
        //
        // SUBWORKFLOW: phylogenetic analysis
        //
        PHYLOGENY (
            ch_samplesheet_fasta
        )
        ch_multiqc_files = ch_multiqc_files.mix(PHYLOGENY.out.multiqc.collect{it[1]})
        ch_versions = ch_versions.mix(PHYLOGENY.out.versions)

    }

    
    //
    // Collate and save software versions
    //
    def topic_versions = Channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'campneu_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        channel.fromPath(params.multiqc_config, checkIfExists: true) :
        channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
