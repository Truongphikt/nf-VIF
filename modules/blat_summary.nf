process BLAT_SUMMARY {

    container "phinguyen2000/pandas:813ad74"
    
    cpus   { check_max( 1, 'cpus' ) }
    memory { check_max( 18.GB * task.attempt, 'memory' ) }
    time   { check_max( 12.h * task.attempt, 'time' ) }

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