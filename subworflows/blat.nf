include    {  BUILD_TEST_BLAT_DB     }               from              "../modules/build_test_blat_db.nf"
include    {  BLAT_SOFT_CLIPPED_SEQ  }               from              "../modules/blat_soft_clipped_seq.nf"
include    {  BLAT_SUMMARY           }               from              "../modules/blat_summary.nf"


workflow BLAT {
    take:
    referenceFastaForIndex
    clipped_seq             // [(val(pfix), val(prefix), path(clipped_seq))]
    bkp_info                // [(val(pfix), path(bkp_info))]

    main:

    if (!params.skipBlat && params.blatdb){
        Channel.fromPath( params.blatdb )
                .ifEmpty { exit 1, "BLAT database not found: ${params.blatdb}" }
                .set { blatDatabase }
    } 
    
    if (!params.skipBlat && !params.blatdb && (workflow.profile =~ /test/)){
        BUILD_TEST_BLAT_DB(referenceFastaForIndex)
        blatDatabase = BUILD_TEST_BLAT_DB.out
    }

    if(params.skipBlat){
        blatDatabase = Channel.empty()
    }


    BLAT_SOFT_CLIPPED_SEQ(
        clipped_seq.combine(blatDatabase)                       // [(val(pfix), val(prefix), path(clipped_seq), path(blatdb))]
    )

    blat_parser_script = Channel.fromPath("https://raw.githubusercontent.com/Truongphikt/nf-VIF/refs/heads/master/src/blatParser.py")
    BLAT_SUMMARY(
        BLAT_SOFT_CLIPPED_SEQ.out.combine(bkp_info, by: 0)      // [(val(pfix), val(prefix), path(blat_res), path(bkp_info), path(blat_parser_script)]
                             .combine(blat_parser_script)
    )

    emit:
    ttd                 =   BLAT_SUMMARY.out.ttd                // ([val(prefix), path(bkptable_filtered)])
}