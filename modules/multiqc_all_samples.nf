process MULTIQC_ALL_SAMPLES {
    publishDir "${params.outdir}/MultiQC/", mode: 'copy'

    when:
    !params.skipMultiqc

    input:
    file splan from chSplan.first()
    file multiqcConfig from chMultiqcConfig.first()
    file hpvConfig from mqcHpvConf.collect().ifEmpty([])
    file('fastqc/*') from fastqcResults.map{items->items[1]}.collect().ifEmpty([])
    file('trimming/*') from trimgaloreResults.map{items->items[1]}.collect().ifEmpty([])
    file ('hpv/*') from hpvGenoStats.map{items->items[1]}.collect()
    file ('hpv/*') from hpvGenoMqc.map{items->items[1]}.collect()
    file ('hpv/*') from hpvCovStats.map{items->items[1]}.collect()
    file ('hpv/*') from hpvBwCov.map{items->items[1]}.collect()
    file ('hpv/*') from bkpPos.map{items->items[1]}.collect()
    file ('hpv/*') from mqcGenepos.map{items->items[1]}.collect()
    file ('hpv/*') from ttd.map{items->items[1]}.collect().ifEmpty([])
    file ('hpv/*') from hpvBowtie2Log.map{items->items[1]}.collect().ifEmpty([])
    file ('ctrl/*') from ctrlStats.map{items->items[1]}.collect()

    file ('software_versions/*') from software_versions_yaml.collect()
    file ('workflow_summary/*') from workflow_summary_yaml.collect()

    output:
    file splan
    file "*multiqc_report.html" into multiqcReport
    file "*_data"

    script:
    rtitle = customRunName ? "--title \"$customRunName\"" : ''
    rfilename = customRunName ? "--filename " + customRunName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''
    metadataOpts = params.metadata ? "--metadata ${metadata}" : ""
    splanOpts = params.samplePlan ? "--splan ${params.samplePlan}" : ""
    """	
    stats2multiqc.sh ${splan}
    mqc_header.py --name "nf-VIF" --version ${workflow.manifest.version} ${metadataOpts} ${splanOpts} > multiqc-config-header.yaml         
    multiqc . -f $rtitle $rfilename -c $multiqcConfig -c $hpvConfig -c multiqc-config-header.yaml -m fastqc -m custom_content
    """
}