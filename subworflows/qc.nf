include {  TRIMGALORE  }                  from            "../modules/trim_galore.nf"
include {  FASTQC      }                  from            "../modules/fastqc.nf"

workflow QC {
    take:

    main:
    if (!params.skipTrimming){
    TRIMGALORE()
    }else{
    readsTrimgalore.into{readsHpvmap; readsSplitmap; readsCtrl; readsFastqc}
    trimgaloreResults = Channel.from(false)
    }

    FASTQC()

    emit:
}