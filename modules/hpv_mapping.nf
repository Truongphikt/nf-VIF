process HPV_MAPPING {
  tag "$prefix"

  input:
  set val(prefix), file(reads) from readsHpvmap
  file index from bwt2IndexHpv.collect()

  output:
  set val(prefix), file("${prefix}_hpvs.bam") into hpvBam
  set val(prefix), file("*bowtie2.log") into hpvBowtie2Log

  script:
  if ( params.singleEnd ){
  """
  bowtie2 --rg-id BMG --rg SM:${prefix} \\
          --very-sensitive --no-unal \\
          -p ${task.cpus} \\
          -x ${index}/${hpvBwt2Base} \\
          -U ${reads} > ${prefix}_hpvs.bam 2> ${prefix}_hpvs_bowtie2.log
  """
  }else{
  """
  bowtie2 --rg-id BMG --rg SM:${prefix} \\
          --very-sensitive --no-unal \\
          -p ${task.cpus} \\
          -x ${index}/${hpvBwt2Base} \\
          -1 ${reads[0]} -2 ${reads[1]} > ${prefix}_hpvs.bam 2> ${prefix}_hpvs_bowtie2.log 
  """
  }
}