include {  PREPROCESSING                 }               from        "../../subworflows/preprocessing.nf"
include {  QC                            }               from        "../../subworflows/qc.nf"
include {  MAPPING                       }               from        "../../subworflows/mapping.nf"
include {  LOCAL_MAPPING                 }               from        "../../subworflows/local_mapping.nf"
include {  EXTRACT_BREAKPOINTS_SEQUENCE  }               from        "../../modules/extract_breakpoints_sequence.nf"
include {  BLAT                          }               from        "../../subworflows/blat.nf"
include {  MULTIQC_PROCESSING            }               from        "../../subworflows/multiqc_processing.nf"

workflow NF_VIF{
    take:
    bwt2RefIndex
    referenceFastaForIndex
    hpvFastaForIndex
    chFastaCtrl
    readsTrimgalore                        // ([val(prefix), listpath(fastq_file)])
    hpv_bwt2_base
    vif_ob


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
        QC.out.trim_fastq,                          // ([val(prefix), listpath(trimmed_fastq)])
        PREPROCESSING.out.bwt2_index_ctrl,          // ([path(bwt2_index_ctrl_folder)])
        PREPROCESSING.out.bwt2_index_hpv,           // ([path(bwt2_index)])
        hpv_bwt2_base
    )

    filtered_sel_hpv_geno = MAPPING.out.sel_hpv_geno
                                   .filter { geno -> vif_ob.checkGenotypes(geno) }                  // ([path(sel_hpv_geno)])
            // .into { hpvGenoFilter; hpvGenoMqcConfig }

    /*
    * Local Mapping for selected genotypes
    */

    LOCAL_MAPPING(
        QC.out.trim_fastq,                          // ([val(prefix), listpath(trimmed_fastq)])
        filtered_sel_hpv_geno,                      // ([path(sel_hpv_geno)])
        PREPROCESSING.out.bwt2_index_hpv_split      // ([path(bwt2_index_split)])
    )

    /*
    * Breakpoint detection
    */
    extract_softclipped_code = Channel.fromPath("https://raw.githubusercontent.com/Truongphikt/nf-VIF/refs/heads/master/src/extractSoftclipped.py")
    EXTRACT_BREAKPOINTS_SEQUENCE(
        LOCAL_MAPPING.out.hpv_soft_bam                           // [(val(prefix), val(hpv), path(local_bam))]
                     .combine(extract_softclipped_code)          // [(val(prefix), val(hpv), path(local_bam), path(extract_softclipped_code))]
    )

    /*
    * Blat
    */  

    BLAT(
        referenceFastaForIndex,
        EXTRACT_BREAKPOINTS_SEQUENCE.out.clipped_seq,            // [(val(pfix), val(prefix), path(clipped_seq))]
        EXTRACT_BREAKPOINTS_SEQUENCE.out.bkp_info                // [(val(pfix), path(bkp_info))]
    )
    // ttd = Channel.from(false)


    /*
    /* MultiQC PROCESSING
    */

    MULTIQC_PROCESSING()

}