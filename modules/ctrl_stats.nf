process CTRL_STATS {
  tag "$prefix"
  publishDir "${params.outdir}/ctrlMapping/", mode: 'copy'

  input:
  set val(prefix), file(bam) from ctrlBam

  output:
  set val(prefix), file('*fsorted_ctrl.{bam,bam.bai}') into ctrlFiltBams
  set val(prefix), file("*ctrl.stats") into ctrlStats

  script:
  peStatus = params.singleEnd ? "0" : "1"
  """
  nbreads=\$(samtools view -c ${bam})
  samtools view -h -q 20 ${bam} | samtools sort -@  ${task.cpus} - > ${prefix}_fsorted_ctrl.bam
  samtools index ${prefix}_fsorted_ctrl.bam
  samtools idxstats ${prefix}_fsorted_ctrl.bam | cut -f1,3 | sort -k2,2nr > ${prefix}_ctrl.stats
  awk -v isPe=${peStatus} -v tot=\$nbreads '\$1!="*"{s=s+\$2} \$1=="*"{\$1="unmapped"; \$2=tot-s} isPe==1{\$2=\$2/2} {printf "%s\\t%.0f\\n",\$1,\$2}' ${prefix}_ctrl.stats > ${prefix}_ctrl_final.stats
  mv ${prefix}_ctrl_final.stats ${prefix}_ctrl.stats
  """
}