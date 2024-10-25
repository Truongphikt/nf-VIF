process FASTQC {
    tag "$name"

    container "phinguyen2000/fastqc:7df4ab6"
    
    cpus          { 3 * task.attempt      }
    memory        { 5.GB * task.attempt   }
    time          { 12.h * task.attempt   } 
    errorStrategy { task.exitStatus in [143,137] ? 'retry' : 'ignore' }

    input:
    tuple val(name), path(trimmed_fastq)
    // set val(name), file(reads) from readsFastqc

    output:
    tuple val(prefix), path("${prefix}*.{zip,html}")
    // set val(prefix), file("${prefix}*.{zip,html}") into fastqcResults

    script:
    prefix = trimmed_fastq[0].toString() - ~/(_1)?(_2)?(_R1)?(_R2)?(.R1)?(.R2)?(_val_1)?(_val_2)?(_trimmed)?(\.fq)?(\.fastq)?(\.gz)?$/
    """
    fastqc -q $trimmed_fastq
    """
}