process MAKE_BOWTIE2_INDEX {
    tag "$refBwt2Base"
    
    input:
    path(reference_fasta)
    //path(fasta) from referenceFastaForIndex

    output:
    path("bowtie2Index/")
    // file "bowtie2Index" into bwt2IndexEnd2end
    // file "bowtie2Index" into bwt2IndexTrim

    script:
    refBwt2Base = reference_fasta.toString() - ~/(\.fa)?(\.fasta)?(\.fas)?$/
    """
    mkdir bowtie2Index
    bowtie2-build ${reference_fasta} bowtie2Index/${refBwt2Base}
    """
}  