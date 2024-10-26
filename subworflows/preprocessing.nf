include {  MAKE_BOWTIE2_INDEX      }        from        "../modules/make_bowtie2_index.nf"
include {  MAKE_BOWTIE2_INDEX_HPV  }        from        "../modules/make_bowtie2_index_hpv.nf"
include {  MAKE_BOWTIE2_INDEX_CTRL }        from        "../modules/make_bowtie2_index_ctrl.nf"


workflow PREPROCESSING{
    take:
    bwt2RefIndex
    referenceFastaForIndex
    hpvFastaForIndex

    main:
    MAKE_BOWTIE2_INDEX(referenceFastaForIndex)
    MAKE_BOWTIE2_INDEX_HPV(hpvFastaForIndex)

    // MAKE_BOWTIE2_INDEX_CTRL()

    // emit:


}

