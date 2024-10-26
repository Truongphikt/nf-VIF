include {  HPV_LOCAL_MAPPING        }             from            "../modules/hpv_local_mapping.nf"
include {  HPV_LOCAL_MAPPING_STATS  }             from            "../modules/hpv_local_mapping_stats.nf"
include {  HPV_COVERAGE             }             from            "../modules/hpv_coverage.nf"


workflow LOCAL_MAPPING {
    take:

    main:
    HPV_LOCAL_MAPPING()
    HPV_LOCAL_MAPPING_STATS()
    HPV_COVERAGE()

    emit:

}