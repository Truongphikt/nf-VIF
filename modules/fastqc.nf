process FASTQC {
    tag "$name"
    publishDir "${params.outdir}/fastqc", mode: 'copy',
        saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}
   
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