include    {  BUILD_TEST_BLAT_DB  }               from              "../modules/build_test_blat_db.nf"
// include  {  BLAT_SOFT_CLIPPED_SEQ  }            from                "../modules/blat_soft_clipped_seq.nf"
// include  {  BLAT_SUMMARY           }            from                "../modules/blat_summary.nf"


workflow BLAT {
    take:
    referenceFastaForIndex

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


    // BLAT_SOFT_CLIPPED_SEQ()
    // BLAT_SUMMARY()

    // emit:
}