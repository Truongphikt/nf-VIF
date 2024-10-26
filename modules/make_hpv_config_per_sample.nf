process MAKE_HPV_CONFIG_PER_SAMPLE {        
    when:
    !params.skipMultiqc

    input:
    file(geno) from hpvGenoMqcConfig
    file genes from chHpvGenesCoord.collect()

    output:
    set val(prefix),file('*conf.mqc') into mqcHpvConf
    set val(prefix),file('*bkp.mqc') optional true into mqcGenepos

    script:
    prefix = geno.toString() - "_HPVgenotyping.filtered"
    """
    awk -F, '{print \$2}' ${geno} | sort -u > allgenotypes_unique.txt
    scrape_mqc_config.py  allgenotypes_unique.txt > ${prefix}_hpv_conf.mqc
    gene_tracks.sh allgenotypes_unique.txt ${genes} ${prefix} 
    """
}