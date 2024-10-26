process BLAT_SOFT_CLIPPED_SEQ {
    
    when:
    !params.skipBlat

    input:
    file(blatdb) from blatDatabase.collect()
    set val(pfix), val(sname), file(fasta) from clippedSeq

    output:
    set val(pfix), val(sname), file("*.tsv") into blatRes

    script:
    """
    blat ${blatdb} ${fasta} ${pfix}.tsv -noHead -minScore=25 -minIdentity=90
    """
}