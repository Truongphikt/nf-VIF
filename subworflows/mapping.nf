include {  CTRL_MAPPING      }              from            "../modules/ctrl_mapping.nf"
include {  CTRL_STATS        }              from            "../modules/ctrl_stats.nf"
include {  HPV_MAPPING       }              from            "../modules/hpv_mapping.nf"
include {  SELECT_GENOTYPES  }              from            "../modules/select_genotypes.nf"

worflow MAPPING {
    take:

    main:
    /*
    * Mapping on control regions
    */

    CTRL_MAPPING()
    CTRL_STATS()

    /*
    * HPV mapping and genotyping
    */ 

    HPV_MAPPING()
    SELECT_GENOTYPES()

    emit:
}