include {  PREPROCESSING                 }               from        "../../subworflows/preprocessing.nf"
include {  QC                            }               from        "../../subworflows/qc.nf"
include {  MAPPING                       }               from        "../../subworflows/mapping.nf"
// include {  LOCAL_MAPPING                 }               from        "../../subworflows/local_mapping.nf"
// include {  EXTRACT_BREAKPOINTS_SEQUENCE  }               from        "../../modules/extract_breakpoints_sequence.nf"
// include {  BLAT                          }               from        "../../subworflows/blat.nf"
// include {  MULTIQC_PROCESSING            }               from        "../../subworflows/multiqc_processing.nf"

workflow NF_VIF{
    take:
    bwt2RefIndex
    referenceFastaForIndex
    hpvFastaForIndex
    chFastaCtrl
    readsTrimgalore                        // ([val(name), listpath(fastq_file)])
    hpv_bwt2_base


    main:
    /****************************************************
    * PRE-PROCESSING
    */

    PREPROCESSING(
        bwt2RefIndex,
        referenceFastaForIndex,
        hpvFastaForIndex,
        chFastaCtrl
    )

    /**************************** Main worflow ******************************


    /*
    * Quality control
    */
    QC(readsTrimgalore)

    /*
    * Mapping
    */

    MAPPING(
        QC.out.trim_fastq,                          // ([val(name), listpath(trimmed_fastq)])
        PREPROCESSING.out.bwt2_index_ctrl,          // ([path(bwt2_index_ctrl_folder)])
        PREPROCESSING.out.bwt2_index_hpv,           // ([path(bwt2_index)])
        hpv_bwt2_base
    )

    // // Filter - removes all samples for which the genotype has not been detected
    // skippedNogeno = []
    // def checkGenotypes(geno) {
    // def nbGeno = 0;
    // geno.eachLine { nbGeno++; }
    // samplename = geno.getBaseName() - '_HPVgenotyping.filered'
    // if(nbGeno < 1 ){
    //     log.info "#################### NO HPV GENOTYPE DETECTED! IGNORING FOR FURTHER DOWNSTREAM ANALYSIS! ($samplename)"
    //     skippedNogeno << samplename
    //     return false
    // } else {
    //     return true
    // }
    // }


    // selHpvGeno
    //         .filter { geno -> checkGenotypes(geno) }
    //         .into { hpvGenoFilter; hpvGenoMqcConfig }

    // /*
    // * Local Mapping for selected genotypes
    // */

    // LOCAL_MAPPING()

    // /*
    // * Breakpoint detection
    // */

    // EXTRACT_BREAKPOINTS_SEQUENCE()

    // /*
    // * Blat
    // */  

    // if (!params.skipBlat){
    //     BLAT()
    // }else{
    //     ttd = Channel.from(false)
    // }


    // /*
    // /* MultiQC PROCESSING
    // */

    // MULTIQC_PROCESSING()

}