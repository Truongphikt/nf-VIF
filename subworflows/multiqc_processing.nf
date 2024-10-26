include {  GET_SOFTWARE_VERSIONS       }                 from                "../modules/get_software_versions.nf"
include {  WORKFLOW_SUMMARY_MQC        }                 from                "../modules/workflow_summary_mqc.nf"
include {  MAKE_HPV_CONFIG_PER_SAMPLE  }                 from                "../modules/make_hpv_config_per_sample.nf"
include {  MULTIQC                     }                 from                "../modules/multiqc.nf"
include {  MAKE_HPV_CONFIG             }                 from                "../modules/make_hpv_config.nf"
include {  MULTIQC_ALL_SAMPLES         }                 from                "../modules/multiqc_all_samples.nf"


workflow MULTIQC_PROCESSING {
    take:

    main:
    GET_SOFTWARE_VERSIONS()
    WORKFLOW_SUMMARY_MQC()


    if (params.splitReport){

    MAKE_HPV_CONFIG_PER_SAMPLE()

    mqcHpvConf
            .join(fastqcResults, remainder: true)
            .join(trimgaloreResults, remainder: true)
            .join(ctrlStats)
        .join(hpvBowtie2Log)
        .join(hpvGenoMqc)
            .join(hpvCovStats.groupTuple(), remainder: true)
            .join(hpvBwCov.groupTuple(), remainder: true)
            .join(bkpPos.groupTuple(), remainder: true)
            .join(hpvGenoStats, remainder: true)
            .join(mqcGenepos, remainder: true)
        .join(ttd.groupTuple(), remainder: true)
        .dump(tag: "join")
            .set{chHpvReport}

    MULTIQC()
    }else{
    MAKE_HPV_CONFIG()
    MULTIQC_ALL_SAMPLES
    }


    emit:


}