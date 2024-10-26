process SELECT_GENOTYPES{
  publishDir "${params.outdir}/hpvMapping/allref", mode: 'copy'

  input:
  set val(prefix), file(bam) from hpvBam 

  output:
  set val(prefix), file("${prefix}_HPVgenotyping.stats") into hpvGenoStats
  set val(prefix), file("${prefix}_HPVgenotyping.filtered") into hpvGenoMqc
  file("${prefix}_HPVgenotyping.filtered") into selHpvGeno
  set val(prefix), file('*fsorted_hpvs.{bam,bam.bai}') into hpvsFiltBams

  script:
  """
  samtools view -h -q ${params.minMapq} ${bam} | samtools sort -@  ${task.cpus} -o ${prefix}_fsorted_hpvs.bam -
  samtools index ${prefix}_fsorted_hpvs.bam
  samtools idxstats ${prefix}_fsorted_hpvs.bam | cut -f1,3 | sort -k2,2nr > ${prefix}_HPVgenotyping.counts
  nbreads=\$(samtools view -c ${bam})
  awk -v tot=\$nbreads '{printf("%s\\t%.02f\\n", \$1, \$2/tot)}' ${prefix}_HPVgenotyping.counts > ${prefix}_HPVgenotyping.freq
  awk -v minFreq=${params.minFreqGeno} -v sname=${prefix} '\$2>=minFreq{print sname","\$1}' ${prefix}_HPVgenotyping.freq > ${prefix}_HPVgenotyping.filtered
  awk '\$2>=10{print}\$2<10{others+=\$2}END{print "Others\t"others}' ${prefix}_HPVgenotyping.counts | grep -v "*" > ${prefix}_HPVgenotyping.stats
  """
}