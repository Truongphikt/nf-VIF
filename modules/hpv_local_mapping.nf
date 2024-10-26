process HPV_LOCAL_MAPPING {
  
  input:
  file index from bwt2IndexHpvSplit.first()
  set val(prefix), val(hpv), file(reads) from hpvGenoFilter
    .splitCsv(header: ["sample", "hpv"])
    .map{
      [ it["sample"], it["hpv"] ]
    }
    .combine(readsSplitmap, by: 0)
    .dump(tag: "hpvloc")

  output:
  set val(prefix), file("*.bam") into hpvLocalBam, hpvCovBam, hpvSoftBam

  script: 
  if ( params.singleEnd ){
  """
  bowtie2 --rg-id BMG --rg SM:${prefix} \\
          --local --very-sensitive-local --no-unal \\
          -p ${task.cpus} \\
          -x ${index}/${hpv} \\
          -U ${reads} > ${prefix}-${hpv}.bam 2> ${prefix}-${hpv}_bowtie2.log
  """
  }else{
  """
  bowtie2 --rg-id BMG --rg SM:${prefix} \\
          --local --very-sensitive-local --no-unal \\
          -p ${task.cpus} \\
          -x ${index}/${hpv} \\
          -1 ${reads[0]} -2 ${reads[1]} > ${prefix}-${hpv}.bam 2> ${prefix}-${hpv}_bowtie2.log
  """
  }
}