process MAKE_BOWTIE2_INDEX_HPV {
    publishDir path: { params.saveReference ? "${params.outdir}/references" : params.outdir },
            saveAs: { params.saveReference ? it : null }, mode: 'copy'

    input:
    file fasta from hpvFastaForIndex

    output:
    file "bowtie2IndexHpv" into bwt2IndexHpv
    file "bowtie2IndexHpvSplit" into bwt2IndexHpvSplit

    script:
    hpvBwt2Base = fasta.toString() - ~/(\.fa)?(\.fasta)?(\.fas)?$/
    """
    mkdir bowtie2IndexHpv
    bowtie2-build ${fasta} bowtie2IndexHpv/${hpvBwt2Base}

    mkdir fastaSplit && cd fastaSplit
    cat ../$fasta | awk '{ if (substr(\$0, 1, 1)==">") {filename=(substr(\$0,2) ".fa")} print \$0 > filename }'
    cd .. && ls fastaSplit/* > listoffasta.txt

    mkdir bowtie2IndexHpvSplit
    while read ff; do
    base="\$(basename \"\${ff}\" | sed -e 's/.fa//')"
    bowtie2-build "\${ff}" bowtie2IndexHpvSplit/"\${base}"
    done < listoffasta.txt
    """
}