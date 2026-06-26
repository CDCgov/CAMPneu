//
// sequence typing with MLST, P1 typing with FastANI
//
include { MLST            } from '../../../modules/local/mlst/main'
include { FASTANI         } from '../../../modules/local/fastani/main'

workflow TYPING {
    take:
    ch_contigs        // channel: [ val(meta), contigs ]

    main:

    ch_versions = Channel.empty()

    //
    // Sequence typing with mlst
    //
    MLST (
        ch_contigs,
        "${params.blastdb}",
        "${params.datadir}",
        "${params.db_name}"
    )
    ch_versions = ch_versions.mix(MLST.out.versions)

    // Merge mlst reports
    ch_merge = MLST.out.tsv
                    .collectFile(name:'mlst_report.tsv', storeDir:"${params.outdir}/reports/", keepHeader:true){
                        meta, file -> file
                    }

    //
    // 1 vs typ1 & type2 reference to determine P1 type with FastANI
    //
    FASTANI (
        ch_contigs,
        "${params.ref_list}",
        "${params.reference_genome}",
        "${params.type2}",
        false,
        true
    )
    ch_versions = ch_versions.mix(FASTANI.out.versions)

    // Get P1 types for references
    ch_ref_type = Channel.fromPath("${params.ref_type}")
                                    .splitCsv( header: true )
                                    .map {
                                        row ->
                                        [ row.ref, row.type ]
                                    }
    
    ch_ref_type = ch_ref_type
                    .map {
                        name, type ->
                        [ [ id:"reference" ], name, type ]
                    }


    // Get top hit name and ANI for each sample
    ch_tophit = FASTANI.out.ani
                    .map {
                        meta, ani_file ->
                        def tophit = ani_file
                                        .splitCsv( sep:"\t" )
                                        .max { a, b -> a[2] <=> b[2] }
                        def name   = tophit[1].tokenize('/').last().split('_')[0..1].join('_')
                        return [ meta, name, tophit[2] as Float ]
                    }

    // Determine P1 type based on top ANI hit
    ch_match_ref_type = ch_tophit
                            .combine(ch_ref_type, by:[1])
                            .map {
                                name, meta, avgnid, ref, type ->
                                [ meta, [ ani:avgnid, hit:ref, p1type:type ] ]
                            }


    // Send all ANI tophits to a file
    ch_ani = ch_match_ref_type
                .collectFile( name:'fastANI_tophits.tsv', newLine:true, storeDir:"${params.outdir}/reports/" ){
                    meta, tophit ->
                    [ 'fastANI_tophits.tsv', meta.id + '\t' + tophit.p1type + '\t' + tophit.ani ]
                }
    

    emit:
    mlst_report         = ch_merge                                   // channel: [ mlst_report.tsv ]

    ani_tophit          = ch_tophit                                  // channel: [ meta, hitname, ani ]
    ani_report          = ch_ani                                     // channel: [ fastANI_tophits.tsv ]

    versions            = ch_versions                                // channel: [ versions.yml ]

}
