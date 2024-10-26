#!/usr/bin/env nextflow

/*
Copyright Institut Curie 2019-2020
This software is a computer program whose purpose is to analyze high-throughput sequencing data.
You can use, modify and/ or redistribute the software under the terms of license (see the LICENSE file for more details).
The software is distributed in the hope that it will be useful, but "AS IS" WITHOUT ANY WARRANTY OF ANY KIND. 
Users are therefore encouraged to test the software's suitability as regards their requirements in conditions enabling the security of their systems and/or data. 
The fact that you are presently reading this means that you have had knowledge of the license and that you accept its terms.

This script is based on the nf-core guidelines. See https://nf-co.re/ for more information
*/


/*
========================================================================================
                         HPV DETECTION PIPELINE
========================================================================================
 HPV Analysis Pipeline.
 #### Homepage / Documentation
 https://gitlab.curie.fr/data-analysis/illumina-hpv
----------------------------------------------------------------------------------------


/*
 * SET UP CONFIGURATION VARIABLES
 */

vif_ob = VIF(workflow, log, params, "$USER")

// Show help emssage
if (params.help){
    vif_ob.helpMessage()
    exit 0
}

// Configure reference genomes
// Reference index path configuration

params.bwt2Index = params.genome ? params.genomes[ params.genome ].bowtie2 ?: false : false
params.fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false
params.blatdb =  params.genome ? params.genomes[ params.genome ].blatdb ?: false : false

params.bwt2IndexHpv = params.fastaHpv ? false : params.genomes['HPV'].bowtie2 ?: false
params.bwt2IndexHpvSplit = params.fastaHpv ? false : params.genomes['HPV'].bowtie2Split ?: false
params.fastaHpv = params.fastaHpv ?: params.genomes['HPV'].fasta ?: false

params.genesHpv = params.genomes['HPV'].genes ?: false
params.fastaCtrl = params.genomes['HPV'].ctrlCapture ?: false



// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
customRunName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  customRunName = workflow.runName
}

// Stage config files
chMultiqcConfig = Channel.fromPath(params.multiqcConfig)
chOutputDocs = Channel.fromPath("$baseDir/docs/output.md")
chFastaCtrl = Channel.fromPath(params.fastaCtrl)
chHpvGenesCoord = Channel.fromPath(params.genesHpv)

/*
 * CHANNELS
 */

/*
 * Create a channel for input read files
 */

if(params.samplePlan){
   if(params.singleEnd){
      Channel
         .from(file("${params.samplePlan}"))
         .splitCsv(header: false)
         .map{ row -> [ row[0], [file(row[2])]] }
         .set {readsTrimgalore}
   }else{
      Channel
         .from(file("${params.samplePlan}"))
         .splitCsv(header: false)
         .map{ row -> [ row[0], [file(row[2]), file(row[3])]] }
         .set {readsTrimgalore}
   }
   params.reads=false
}
else if(params.readPaths){
    if(params.singleEnd){
        Channel
            .from(params.readPaths)
            .map { row -> [ row[0], [file(row[1][0])]] }
            .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied" }
            .set {readsTrimgalore}
    } else {
        Channel
            .from(params.readPaths)
            .map { row -> [ row[0], [file(row[1][0]), file(row[1][1])]] }
            .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied" }
            .set {readsTrimgalore}
    }
} else {
    Channel
        .fromFilePairs( params.reads, size: params.singleEnd ? 1 : 2 )
        .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nNB: Path requires at least one * wildcard!\nIf this is single-end data, please specify --singleEnd on the command line." }
        .set {readsTrimgalore}
}

/*
 * Make sample plan if not available
 */

if (params.samplePlan){
  chSplan = Channel.fromPath(params.samplePlan)
}else if(params.readPaths){
  if (params.singleEnd){
    Channel
       .from(params.readPaths)
       .collectFile() {
         item -> ["samplePlan.csv", item[0] + ',' + item[0] + ',' + item[1][0] + '\n']
        }
       .set{ chSplan }
  }else{
     Channel
       .from(params.readPaths)
       .collectFile() {
         item -> ["samplePlan.csv", item[0] + ',' + item[0] + ',' + item[1][0] + ',' + item[1][1] + '\n']
        }
       .set{ chSplan }
  }
}else{
  if (params.singleEnd){
    Channel
       .fromFilePairs( params.reads, size: 1 )
       .collectFile() {
          item -> ["samplePlan.csv", item[0] + ',' + item[0] + ',' + item[1][0] + '\n']
       }     
       .set { chSplan }
  }else{
    Channel
       .fromFilePairs( params.reads, size: 2 )
       .collectFile() {
          item -> ["samplePlan.csv", item[0] + ',' + item[0] + ',' + item[1][0] + ',' + item[1][1] + '\n']
       }     
       .set { chSplan }
   }
}

/*
 * Other input channels
 */

// Reference genome

if ( params.bwt2Index ){
   lastPath = params.bwt2Index.lastIndexOf(File.separator)
   refBwt2Dir =  params.bwt2Index.substring(0,lastPath+1)
   refBwt2Base = params.bwt2Index.substring(lastPath+1)

   Channel.fromPath( refBwt2Dir , checkIfExists: true)
      .ifEmpty { exit 1, "Genome index: Provided index not found: ${params.bwt2Index}" }
      .set { bwt2RefIndex }

}
else if ( params.fasta ) {
   lastPath = params.fasta.lastIndexOf(File.separator)
   refBwt2Base = params.fasta.substring(lastPath+1)

   Channel.fromPath( params.fasta )
        .ifEmpty { exit 1, "Genome index: Fasta file not found: ${params.fasta}" }
        .set { referenceFastaForIndex }
}
else {
   exit 1, "No reference genome specified!"
}

if ( params.blatdb ){
   Channel.fromPath( params.blatdb )
        .ifEmpty { exit 1, "BLAT database not found: ${params.blatdb}" }
        .set { blatDatabase }
}

//HPV genome

if ( params.bwt2IndexHpv && params.bwt2IndexHpvSplit ){

   lastPath = params.bwt2IndexHpv.lastIndexOf(File.separator)
   hpvBwt2Dir =  params.bwt2IndexHpv.substring(0,lastPath+1)
   hpvBwt2Base = params.bwt2IndexHpv.substring(lastPath+1)

   Channel.fromPath( hpvBwt2Dir , checkIfExists: true)
      .ifEmpty { exit 1, "HPV index: Provided index not found: ${params.bwt2IndexHpv}" }
      .set { bwt2IndexHpv }

   lastPath = params.bwt2IndexHpvSplit.lastIndexOf(File.separator)
   hpvSplitBwt2Dir =  params.bwt2IndexHpvSplit.substring(0,lastPath+1)
 
   Channel.fromPath( hpvSplitBwt2Dir , checkIfExists: true)
      .ifEmpty { exit 1, "HPV index per strain: Provided index not found: ${params.bwt2IndexHpvSplit}" }
      .set { bwt2IndexHpvSplit }
}
else if ( params.fastaHpv ){
   lastPath = params.fastaHpv.lastIndexOf(File.separator)
   hpvBwt2Base = params.fastaHpv.substring(lastPath+1) - ~/(\.fa)?(\.fasta)?(\.fas)?$/

   Channel.fromPath( params.fastaHpv )
        .ifEmpty { exit 1, "HPV index: Fasta file not found: ${params.fastaHpv}" }
        .set { hpvFastaForIndex }
}
else{
   exit 1, "No HPV genome specified!"
}

// Header log info
vif_ob.headerInfo()


/****************************************************
 * PRE-PROCESSING
 */

if ( !params.bwt2Index && params.fasta ){
  MAKE_BOWTIE2_INDEX() 
}

if ( (!params.bwt2IndexHpv | !params.bwt2IndexHpvSplit) && params.fastaHpv ){
  MAKE_BOWTIE2_INDEX_HPV()
}

MAKE_BOWTIE2_INDEX_CTRL()



/****************************************************
 * Main worflow
 */


/*
 * Reads Trimming
 */
if (!params.skipTrimming){
   TRIMGALORE()
}else{
   readsTrimgalore.into{readsHpvmap; readsSplitmap; readsCtrl; readsFastqc}
   trimgaloreResults = Channel.from(false)
}


/*
 * FastQC
 */
FASTQC()

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


// Filter - removes all samples for which the genotype has not been detected
skippedNogeno = []
def checkGenotypes(geno) {
  def nbGeno = 0;
  geno.eachLine { nbGeno++; }
  samplename = geno.getBaseName() - '_HPVgenotyping.filered'
  if(nbGeno < 1 ){
      log.info "#################### NO HPV GENOTYPE DETECTED! IGNORING FOR FURTHER DOWNSTREAM ANALYSIS! ($samplename)"
      skippedNogeno << samplename
      return false
  } else {
      return true
  }
}


selHpvGeno
        .filter { geno -> checkGenotypes(geno) }
    	.into { hpvGenoFilter; hpvGenoMqcConfig }

/*
 * Local Mapping for selected genotypes
 */

HPV_LOCAL_MAPPING()
HPV_LOCAL_MAPPING_STATS()
HPV_COVERAGE()


/*
 * Breakpoint detection
 */

EXTRACT_BREAKPOINTS_SEQUENCE()

/*
 * Blat
 */  

if (!params.skipBlat){
   BLAT_SOFT_CLIPPED_SEQ()
   BLAT_SUMMARY()
}else{
   ttd = Channel.from(false)
}


/*
/* MultiQC
*/


GET_SOFTWARE_VERSIONS()
WORKFLOW_SUMMARY_MQC()


if (params.splitReport){

   MAKE_HPV_CONFIG_PER_SAMPLE()

   mqcHpvConf
        .join(fastqcResults, remainder: true)
        .join(trimgaloreResults, remainder: true)
        .join(ctrlStats)
	.join(hpvBowtie2Log)
	.join(hpvGenoMqc)
        .join(hpvCovStats.groupTuple(), remainder: true)
        .join(hpvBwCov.groupTuple(), remainder: true)
        .join(bkpPos.groupTuple(), remainder: true)
        .join(hpvGenoStats, remainder: true)
        .join(mqcGenepos, remainder: true)
	.join(ttd.groupTuple(), remainder: true)
	.dump(tag: "join")
        .set{chHpvReport}

  MULTIQC()
}else{
  MAKE_HPV_CONFIG()
  MULTIQC_ALL_SAMPLES
}


