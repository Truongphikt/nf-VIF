process MULTIQC_ALL_SAMPLES {
    container = "phinguyen2000/multiqc:c42a7c6"

    when:
    !params.skipMultiqc

    input:
    tuple path(splan), path(multiqcConfig), path(hpvConfig),
          path('fastqc/*'), path('trimming/*'),
          path('hpv/*'), path('hpv/*'), path('hpv/*'),
          path('hpv/*'), path('hpv/*'), path('hpv/*'),
          path('hpv/*'), path('hpv/*'),
          path('ctrl/*'), 
          path('software_versions/*'), path('workflow_summary/*'),
          val(customRunName),
          path(stats2_multiqc_script), path(mqc_header_script)

    // file splan from chSplan.first()
    // file multiqcConfig from chMultiqcConfig.first()
    // file hpvConfig from mqcHpvConf.collect().ifEmpty([])
    // file('fastqc/*') from fastqcResults.map{items->items[1]}.collect().ifEmpty([])
    // file('trimming/*') from trimgaloreResults.map{items->items[1]}.collect().ifEmpty([])
    // file ('hpv/*') from hpvGenoStats.map{items->items[1]}.collect()
    // file ('hpv/*') from hpvGenoMqc.map{items->items[1]}.collect()
    // file ('hpv/*') from hpvCovStats.map{items->items[1]}.collect()
    // file ('hpv/*') from hpvBwCov.map{items->items[1]}.collect()
    // file ('hpv/*') from bkpPos.map{items->items[1]}.collect()
    // file ('hpv/*') from mqcGenepos.map{items->items[1]}.collect()
    // file ('hpv/*') from ttd.map{items->items[1]}.collect().ifEmpty([])
    // file ('hpv/*') from hpvBowtie2Log.map{items->items[1]}.collect().ifEmpty([])
    // file ('ctrl/*') from ctrlStats.map{items->items[1]}.collect()

    // file ('software_versions/*') from software_versions_yaml.collect()
    // file ('workflow_summary/*') from workflow_summary_yaml.collect()

    output:
    path("*multiqc_report.html"),       emit: multiqc_report
    path("*_data"),                     emit: data

    // file "*multiqc_report.html" into multiqcReport
    // file "*_data"

    script:
    rtitle = customRunName ? "--title \"$customRunName\"" : ''
    rfilename = customRunName ? "--filename " + customRunName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''
    metadataOpts = params.metadata ? "--metadata ${metadata}" : ""
    splanOpts = params.samplePlan ? "--splan ${params.samplePlan}" : ""
    """	
    bash $stats2_multiqc_script ${splan}
    python $mqc_header_script --name "nf-VIF" --version ${workflow.manifest.version} ${metadataOpts} ${splanOpts} > multiqc-config-header.yaml         
    multiqc . -f $rtitle $rfilename -c $multiqcConfig -c $hpvConfig -c multiqc-config-header.yaml -m fastqc -m custom_content
    """
}