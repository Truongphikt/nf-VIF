process TRIMGALORE {
    tag "$name" 

    container "phinguyen2000/trim-galore:231627f"

    input:
    set val(name), file(reads) from readsTrimgalore

    output:
    set val(name), file("*fq.gz") into readsHpvmap, readsSplitmap, readsCtrl, readsFastqc
    set val(prefix), file("*trimming_report.txt") into trimgaloreResults

    script:
    prefix = reads[0].toString() - ~/(_1)?(_2)?(_R1)?(_R2)?(.R1)?(.R2)?(_val_1)?(_val_2)?(\.fq)?(\.fastq)?(\.gz)?$/
    if (params.singleEnd) {
    """
    trim_galore --trim-n --quality 20 --length 20\
            --gzip $reads --basename ${prefix} --cores ${task.cpus}
    """
    }else {
    """
    trim_galore --trim-n --quality 20  --length 20 \
            --paired --gzip $reads --basename ${prefix} --cores ${task.cpus}
    mv ${prefix}_R1_val_1.fq.gz ${prefix}_R1_trimmed.fq.gz
    mv ${prefix}_R2_val_2.fq.gz ${prefix}_R2_trimmed.fq.gz
    mv ${reads[0]}_trimming_report.txt ${prefix}_R1_trimming_report.txt
    mv ${reads[1]}_trimming_report.txt ${prefix}_R2_trimming_report.txt
    """
    }
}