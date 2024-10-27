include {  GET_SOFTWARE_VERSIONS       }                 from                "../modules/get_software_versions.nf"
include {  WORKFLOW_SUMMARY_MQC        }                 from                "../modules/workflow_summary_mqc.nf"
// include {  MAKE_HPV_CONFIG_PER_SAMPLE  }                 from                "../modules/make_hpv_config_per_sample.nf"
// include {  MULTIQC                     }                 from                "../modules/multiqc.nf"
include {  MAKE_HPV_CONFIG             }                 from                "../modules/make_hpv_config.nf"
// include {  MULTIQC_ALL_SAMPLES         }                 from                "../modules/multiqc_all_samples.nf"


workflow MULTIQC_PROCESSING {
    take:
    vif_ob
    filtered_sel_hpv_geno                   // ([path(sel_hpv_geno)])
    ch_hpv_genes_coord

    main:
    scrape_software_versions_script = Channel.fromPath("https://raw.githubusercontent.com/Truongphikt/nf-VIF/refs/heads/master/src/scrape_software_versions.py")
    GET_SOFTWARE_VERSIONS(scrape_software_versions_script)
    
    mqc_yaml_content = params.skipMultiqc? Channel.empty() : Channel.of(vif_ob.getYamlContent())
    WORKFLOW_SUMMARY_MQC(mqc_yaml_content)


    // if (params.splitReport){

    // MAKE_HPV_CONFIG_PER_SAMPLE()

    // mqcHpvConf
    //         .join(fastqcResults, remainder: true)
    //         .join(trimgaloreResults, remainder: true)
    //         .join(ctrlStats)
    //     .join(hpvBowtie2Log)
    //     .join(hpvGenoMqc)
    //         .join(hpvCovStats.groupTuple(), remainder: true)
    //         .join(hpvBwCov.groupTuple(), remainder: true)
    //         .join(bkpPos.groupTuple(), remainder: true)
    //         .join(hpvGenoStats, remainder: true)
    //         .join(mqcGenepos, remainder: true)
    //     .join(ttd.groupTuple(), remainder: true)
    //     .dump(tag: "join")
    //         .set{chHpvReport}

    // MULTIQC()
    // }else{
    scrape_mqc_config_script = Channel.fromPath("https://raw.githubusercontent.com/Truongphikt/nf-VIF/refs/heads/master/src/scrape_mqc_config.py")
    gene_tracks_script       = Channel.fromPath("https://raw.githubusercontent.com/Truongphikt/nf-VIF/refs/heads/master/src/gene_tracks.sh")
    MAKE_HPV_CONFIG(
        filtered_sel_hpv_geno.combine(ch_hpv_genes_coord)
                             .combine(scrape_mqc_config_script)
                             .combine(gene_tracks_script)
    )
    // MULTIQC_ALL_SAMPLES()
    // }


    // emit:


}