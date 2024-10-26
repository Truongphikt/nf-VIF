include {  MAKE_BOWTIE2_INDEX      }        from        "../modules/make_bowtie2_index.nf"
include {  MAKE_BOWTIE2_INDEX_HPV  }        from        "../modules/make_bowtie2_index_hpv.nf"
include {  MAKE_BOWTIE2_INDEX_CTRL }        from        "../modules/make_bowtie2_index_ctrl.nf"


workflow PREPROCESSING{
    take:
    bwt2RefIndex
    referenceFastaForIndex
    hpvFastaForIndex
    chFastaCtrl

    main:
    MAKE_BOWTIE2_INDEX(referenceFastaForIndex)
    MAKE_BOWTIE2_INDEX_HPV(hpvFastaForIndex)
    MAKE_BOWTIE2_INDEX_CTRL(chFastaCtrl)

    emit:
    bwt2_index_ctrl      = MAKE_BOWTIE2_INDEX_CTRL.out                  // ([path(bwt2_index_ctrl_folder)])
    bwt2_index_hpv       = MAKE_BOWTIE2_INDEX_HPV.out.main              // ([path(bwt2_index)])
    bwt2_index_hpv_split = MAKE_BOWTIE2_INDEX_HPV.out.split             // ([path(bwt2_index_split)])
}

