process CTRL_MAPPING {
  tag "$prefix"
  publishDir "${params.outdir}/ctrlMapping/", mode: 'copy',
      saveAs: {filename ->
          if (filename.endsWith(".log")) "logs/$filename"
          else if (params.saveAlignedIntermediates) filename
	  else null}

  input:
  set val(prefix), file(reads) from readsCtrl
  file index from  bwt2IndexCtrl.collect()

  output:
  set val(prefix), file("${prefix}_ctrl.bam") into ctrlBam
  set val(prefix), file("*ctrl_bowtie2.log") into ctrlBowtie2Log

  script:
  if ( params.singleEnd ){
  """
  bowtie2 --rg-id BMG --rg SM:${prefix} \\
          --very-sensitive \\
          -p ${task.cpus} \\
          -x ${index}/ctrlRegions \\
          -U ${reads} > ${prefix}_ctrl.bam 2> ${prefix}_ctrl_bowtie2.log
  """
  }else{
  """
  bowtie2 --rg-id BMG --rg SM:${prefix} \\
          --very-sensitive \\
          -p ${task.cpus} \\
          -x ${index}/ctrlRegions \\
          -1 ${reads[0]} -2 ${reads[1]} > ${prefix}_ctrl.bam 2> ${prefix}_ctrl_bowtie2.log
  """
  }
}