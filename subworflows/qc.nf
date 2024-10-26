include {  TRIMGALORE  }                  from            "../modules/trim_galore.nf"
include {  FASTQC      }                  from            "../modules/fastqc.nf"

workflow QC {
    take:
    readsTrimgalore                        // ([val(name), listpath(fastq_file)])

    main:
    TRIMGALORE(readsTrimgalore)
    // if (!params.skipTrimming){
    // TRIMGALORE(readsTrimgalore)
    // }else{
    // readsTrimgalore.into{readsHpvmap; readsSplitmap; readsCtrl; readsFastqc}
    // trimgaloreResults = Channel.from(false)
    // }

    in_fastqc = params.skipFastqc ? Channel.empty(): TRIMGALORE.out.trim_fastq // ([val(name), listpath(trimmed_fastq)])
    FASTQC(in_fastqc)

    emit:
    trim_fastq = TRIMGALORE.out.trim_fastq                  // ([val(name), listpath(trimmed_fastq)])
    trimming_report = TRIMGALORE.out.trimming_report        // ([val(name), path(trimming_report)])
}