process MAKE_BOWTIE2_INDEX_CTRL {
  publishDir path: { params.saveReference ? "${params.outdir}/references" : params.outdir },
            saveAs: { params.saveReference ? it : null }, mode: 'copy'

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