process BLAT_SUMMARY {
    publishDir "${params.outdir}/hpvMapping/blat", mode: 'copy'

    input:
    set val(pfix), val(sname), file(psl), file(csv) from blatRes.join(bkpInfo).dump(tag:"blat")

    output:
    set val(sname), file("*_table_filtered.csv") into ttdHQ
    set val(sname), file("*_table.csv") into tableHQ
    set val(sname), file("*_bkptable_filtered.csv") into ttd
    set val(sname), file("*_bkptable.csv") into table

    script:
    """
    blatParser.py -f ${psl} -b ${csv} --sname ${sname}
    """
}