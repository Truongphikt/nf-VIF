process HPV_COVERAGE {
  
  input:
  set val(prefix), file(bam) from hpvCovBam

  output:
  set val(prefix), file("*covmatrix.mqc") into hpvBwCov
  set val(prefix), file('*sorted.{bam,bam.bai}') into hpvSortedBams

  script:
  pfix= bam.toString() - ~/(_sorted)?(.bam)?$/
  normOpts = params.splitReport ? "" : "--normalizeUsing CPM"
  """
  samtools sort -@ ${task.cpus} -o ${pfix}_sorted.bam ${bam}
  samtools index ${pfix}_sorted.bam
  bamCoverage -b ${pfix}_sorted.bam --binSize 50 ${normOpts} --outFileFormat bedgraph -o ${pfix}.bedgraph
  awk -F"\t" '{OFS="\t"; print \$2+25,\$4}' ${pfix}.bedgraph > ${pfix}_covmatrix.mqc
  """
}