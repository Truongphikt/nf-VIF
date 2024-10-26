include {  CTRL_MAPPING      }              from            "../modules/ctrl_mapping.nf"
// include {  CTRL_STATS        }              from            "../modules/ctrl_stats.nf"
// include {  HPV_MAPPING       }              from            "../modules/hpv_mapping.nf"
// include {  SELECT_GENOTYPES  }              from            "../modules/select_genotypes.nf"

workflow MAPPING {
    take:
    trim_fastq                          // ([val(name), listpath(trimmed_fastq)])
    bwt2_index_ctrl                     // ([path(bwt2_index_ctrl_folder)])

    main:
    /*
    * Mapping on control regions
    */

    CTRL_MAPPING(
        trim_fastq.combine(bwt2_index_ctrl)       // ([val(name), listpath(trimmed_fastq), path(bwt2_index_ctrl_folder)])
    )
    // CTRL_STATS()

    // /*
    // * HPV mapping and genotyping
    // */ 

    // HPV_MAPPING()
    // SELECT_GENOTYPES()

    // emit:
}