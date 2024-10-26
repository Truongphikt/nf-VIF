include {  MAKE_BOWTIE2_INDEX      }        from        "../modules/make_bowtie2_index.nf"
include {  MAKE_BOWTIE2_INDEX_HPV  }        from        "../modules/make_bowtie2_index_hpv.nf"
include {  MAKE_BOWTIE2_INDEX_CTRL }        from        "../modules/make_bowtie2_index_ctrl.nf"


workflow PREPROCESSING{
    take:
    bwt2RefIndex
    referenceFastaForIndex

    main:
    MAKE_BOWTIE2_INDEX(referenceFastaForIndex)

    // if ( (!params.bwt2IndexHpv | !params.bwt2IndexHpvSplit) && params.fastaHpv ){
    //     MAKE_BOWTIE2_INDEX_HPV()
    // }

    // MAKE_BOWTIE2_INDEX_CTRL()

    // emit:


}

