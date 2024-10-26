process MAKE_HPV_CONFIG {

    when:
    !params.skipMultiqc

    input:
    file(geno) from hpvGenoMqcConfig.collect()
    file genes from chHpvGenesCoord.collect()

    output:
    file('*conf.mqc') into mqcHpvConf
    file('*bkp.mqc') into mqcGenepos

    script:
    prefix = geno.toString() - "_HPVgenotyping.filtered"
    """
    awk -F, '{print \$2}' ${geno} | sort -u > allgenotypes_unique.txt
    scrape_mqc_config.py  allgenotypes_unique.txt > hpv_conf.mqc
    gene_tracks.sh allgenotypes_unique.txt ${genes} 'genes' 
    """
}