process MAKE_BOWTIE2_INDEX {
    tag "$refBwt2Base"
    publishDir path: { params.saveReference ? "${params.outdir}/references" : params.outdir },
            saveAs: { params.saveReference ? it : null }, mode: 'copy'
    
    input:
    file fasta from referenceFastaForIndex

    output:
    file "bowtie2Index" into bwt2IndexEnd2end
    file "bowtie2Index" into bwt2IndexTrim

    script:
    refBwt2Base = fasta.toString() - ~/(\.fa)?(\.fasta)?(\.fas)?$/
    """
    mkdir bowtie2Index
    bowtie2-build ${fasta} bowtie2Index/${refBwt2Base}
    """
}  