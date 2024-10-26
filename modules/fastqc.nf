process FASTQC {
    tag "$name"

    container "biocontainers/fastqc:v0.11.8dfsg-2-deb_cv1"
    
    cpus     { check_max( 1, 'cpus' ) }
    memory   { check_max( 10.GB * task.attempt, 'memory' ) }
    time     { check_max( 12.h * task.attempt, 'time' ) } 
    errorStrategy  { task.exitStatus in [143,137] ? 'retry' : 'ignore' }
   
    when:
    !params.skipFastqc

    input:
    set val(name), file(reads) from readsFastqc

    output:
    set val(prefix), file("${prefix}*.{zip,html}") into fastqcResults

    script:
    prefix = reads[0].toString() - ~/(_1)?(_2)?(_R1)?(_R2)?(.R1)?(.R2)?(_val_1)?(_val_2)?(_trimmed)?(\.fq)?(\.fastq)?(\.gz)?$/
    """
    fastqc -q $reads
    """
}