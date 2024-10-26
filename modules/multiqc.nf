process MULTIQC {
    container "phinguyen2000/multiqc:48b6b9d"

    when:
    !params.skipMultiqc

    input:
    file splan from chSplan.first()
    file multiqcConfig from chMultiqcConfig.first()
    set val(prefix), file('mqc/hpv_config.mqc'), file('fastqc/*'), file('trimming/*'), file('ctrl/*'), 
    file('hpv/*'), file('hpv/*'), file('hpv/*'), file('hpv/*'), file('hpv/*'), file('hpv/*'), file('hpv/*'), file('hpv/*') from chHpvReport.dump(tag: "mqc") 
    file ('software_versions/*') from software_versions_yaml.collect()
    file ('workflow_summary/*') from workflow_summary_yaml.collect()

    output:
    file splan
    file "*.html" into multiqc_report
    file "*_data"

    script:
    rtitle = customRunName ? "--title \"$customRunName\"" : ''
    rfilename = customRunName ? "--filename " + customRunName.replaceAll('\\W','_').replaceAll('_+','_') + "_" + prefix + "_multiqc_report" : ''
    metadataOpts = params.metadata ? "--metadata ${metadata}" : ""
    splanOpts = params.samplePlan ? "--splan ${params.samplePlan}" : ""
    """
    awk -F"," -v sname=${prefix} '\$1==sname{print}' ${splan} > splan_${prefix}.csv
    stats2multiqc.sh splan_${prefix}.csv	
    mqc_header.py --name "${prefix} - nf-VIF" --version ${workflow.manifest.version} ${metadataOpts} ${splanOpts} > multiqc-config-header.yaml
    multiqc . -f $rtitle -n ${prefix}_nfvif_report.html -c $multiqcConfig -c 'mqc/hpv_config.mqc' -c multiqc-config-header.yaml -m fastqc -m custom_content
    """
}