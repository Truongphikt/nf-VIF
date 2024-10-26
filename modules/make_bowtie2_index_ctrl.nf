process MAKE_BOWTIE2_INDEX_CTRL {
  
  input:
  file fasta from chFastaCtrl

  output:
  file "bowtie2IndexCtrl" into bwt2IndexCtrl

  script:
  """
  mkdir bowtie2IndexCtrl
  bowtie2-build ${fasta} bowtie2IndexCtrl/ctrlRegions
  """
}