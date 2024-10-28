include {  TRIMGALORE  }                  from            "../modules/trim_galore.nf"
include {  FASTQC      }                  from            "../modules/fastqc.nf"

workflow QC {
    take:
    readsTrimgalore                        // ([val(sname), listpath(fastq_file)])

    main:
    if (!params.skipTrimming){
        TRIMGALORE(readsTrimgalore)
        trim_fastq = TRIMGALORE.out.trim_fastq
        trimming_report = TRIMGALORE.out.trimming_report
    }else{
        trim_fastq = readsTrimgalore
        trimming_report = Channel.empty()
    }
 
    FASTQC(
        TRIMGALORE.out.trim_fastq           // ([val(name), listpath(trimmed_fastq)])
    )

    emit:
    trim_fastq = trim_fastq                  // ([val(sname), listpath(trimmed_fastq)])
    trimming_report = trimming_report        // ([val(sname), path(trimming_report)])
    fastqc_results = FASTQC.out              // ([val(sname), path(fastqc_rs)])
}