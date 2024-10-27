include {  GET_SOFTWARE_VERSIONS       }                 from                "../modules/get_software_versions.nf"
include {  WORKFLOW_SUMMARY_MQC        }                 from                "../modules/workflow_summary_mqc.nf"
include {  MAKE_HPV_CONFIG_PER_SAMPLE  }                 from                "../modules/make_hpv_config_per_sample.nf"
include {  MULTIQC                     }                 from                "../modules/multiqc.nf"
include {  MAKE_HPV_CONFIG             }                 from                "../modules/make_hpv_config.nf"
include {  MULTIQC_ALL_SAMPLES         }                 from                "../modules/multiqc_all_samples.nf"


workflow MULTIQC_PROCESSING {
    take:
    vif_ob
    filtered_sel_hpv_geno                    // ([path(sel_hpv_geno)])
    ch_hpv_genes_coord
    chSplan
    chMultiqcConfig
    fastqc_results                           // ([val(sname), path(fastqc_rs)])
    trimming_report                          // ([val(sname), path(trimming_report)])
    ctrl_stats                               // ([val(prefix), path(ctrl_stats)])
    hpv_bowtie2_log                          // ([val(prefix), path(hpv_bowtie2_log)])
    hpv_geno_mqc                             // ([val(prefix), path(sel_hpv_geno)])
    hpv_geno_stats                           // ([val(prefix), path(hpv_geno_stats)])
    hpv_cov_stats                            // [(val(prefix), path(coverage_stats))]
    hpv_bw_cov                               // [(val(prefix), path(covmatrix_mqc))]
    bkp_pos                                  // [(val(prefix), path(3prime_bkp_mqc)]
    ttd                                      // ([val(prefix), path(bkptable_filtered)])
    customRunName

    main:
    scrape_software_versions_script = Channel.fromPath("https://raw.githubusercontent.com/Truongphikt/nf-VIF/refs/heads/master/src/scrape_software_versions.py")
    GET_SOFTWARE_VERSIONS(scrape_software_versions_script)
    
    mqc_yaml_content = params.skipMultiqc? Channel.empty() : Channel.of(vif_ob.getYamlContent())
    WORKFLOW_SUMMARY_MQC(mqc_yaml_content)


    scrape_mqc_config_script = Channel.fromPath("https://raw.githubusercontent.com/Truongphikt/nf-VIF/refs/heads/master/src/scrape_mqc_config.py")
    gene_tracks_script       = Channel.fromPath("https://raw.githubusercontent.com/Truongphikt/nf-VIF/refs/heads/master/src/gene_tracks.sh")
    stats2_multiqc_script    = Channel.fromPath("https://raw.githubusercontent.com/Truongphikt/nf-VIF/refs/heads/master/src/stats2multiqc.sh")
    mqc_header_script        = Channel.fromPath("https://raw.githubusercontent.com/Truongphikt/nf-VIF/refs/heads/master/src/mqc_header.py")


    // if (params.splitReport){

    MAKE_HPV_CONFIG_PER_SAMPLE(
        filtered_sel_hpv_geno.combine(ch_hpv_genes_coord)
                             .combine(scrape_mqc_config_script)
                             .combine(gene_tracks_script)
    )

    MAKE_HPV_CONFIG_PER_SAMPLE.out.mqc_hpv_conf
                                .join(fastqc_results, remainder: true)
                                .join(trimming_report, remainder: true)
                                .join(ctrl_stats)
                                .join(hpv_bowtie2_log)
                                .join(hpv_geno_mqc)
                                .join(hpv_cov_stats.groupTuple(), remainder: true)
                                .join(hpv_bw_cov.groupTuple(), remainder: true)
                                .join(bkp_pos.groupTuple(), remainder: true)
                                .join(hpv_geno_stats, remainder: true)
                                .join(MAKE_HPV_CONFIG_PER_SAMPLE.out.mqc_genepos, remainder: true)
                                .join(ttd.groupTuple(), remainder: true)
                                .dump(tag: "join")
                                .set{ch_hpv_report}

    MULTIQC(
        chSplan.first().combine(chMultiqcConfig.first())
                        .combine(ch_hpv_report)
                        .combine(GET_SOFTWARE_VERSIONS.out)
                        .combine(WORKFLOW_SUMMARY_MQC.out)
                        .combine(customRunName)
                        .combine(stats2_multiqc_script)
                        .combine(mqc_header_script)
    )
    // }else{
    
    MAKE_HPV_CONFIG(
        filtered_sel_hpv_geno.combine(ch_hpv_genes_coord)
                             .combine(scrape_mqc_config_script)
                             .combine(gene_tracks_script)
    )
    qc_all_sample_ch = chSplan.combine(chMultiqcConfig)
    // MULTIQC_ALL_SAMPLES()
    // }


    // emit:


}