process MAKE_BOWTIE2_INDEX_CTRL {
  
  input:
  path fasta_ctrl
  // file fasta from chFastaCtrl

  output:
  path("bowtie2IndexCtrl/")
  // file "bowtie2IndexCtrl" into bwt2IndexCtrl

  script:
  """
  mkdir bowtie2IndexCtrl
  bowtie2-build ${fasta_ctrl} bowtie2IndexCtrl/ctrlRegions
  """
}