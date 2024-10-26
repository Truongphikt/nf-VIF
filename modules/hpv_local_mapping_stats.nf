process HPV_LOCAL_MAPPING_STATS {
  publishDir "${params.outdir}/hpvMapping/pergenotype", mode: 'copy'
  
  input:
  set val(prefix), file(bam) from hpvLocalBam

  output:
  set val(prefix), file("*_coverage.stats") into hpvCovStats

  script:
  pfix= bam.toString() - ~/(.bam)?$/
  """
  genomeCoverageBed -d -ibam ${bam} > ${pfix}_coverage.out
  nbaln=\$(samtools flagstat ${bam} | grep 'mapped (' | awk '{print \$1}')
  echo -e "ID,sample,HPVsubtype,mappedReads,minCov,maxCov,meanCov" > ${pfix}_coverage.stats
  awk -v id=${prefix} -v nbaln=\$nbaln 'NR==1{hpv=\$1;mn=mx=\$3}{total+=\$3}(\$3>mx){mx=\$3}(\$3<mn){mn=\$3} END{OFS=","; print id"_"hpv,id,hpv,nbaln,mn,mx,total/NR}' ${pfix}_coverage.out >> ${pfix}_coverage.stats
  """
}