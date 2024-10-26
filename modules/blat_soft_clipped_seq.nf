process BLAT_SOFT_CLIPPED_SEQ {

    container "phinguyen2000/blat:255336f"
    
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