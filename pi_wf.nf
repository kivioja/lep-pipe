#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * estimate nucleotide diversity pi
 */



 // default 
params.javaheapsize = "64G"

/*
  include processes 
*/
include { generateRefContigChunks; callVarsAndParents } from './modules/callVarsToLgs.nf'
include { calcNucDiversity; collectNucDiversities } from './modules/calcFeatures.nf'
  
// default file names
// - pi calculated from the parents but children aid in parent variant calling
params.mapping_bam2ind_file = params.species + "_bam2indname.csv"
params.bam2indpath = params.sampledir + "/" + params.mapping_bam2ind_file
params.pedfile =  params.species + "_ped.tsv"
params.pedpath = params.sampledir + "/" + params.pedfile

workflow {

  // pi estimation needs also non-variant positions, thus
  // variant calling files will be large, thus do by contig chunk
  // TODO: remove the large chunk files
    
  // divide the contings into junks to be processed in parallel
  contigchunks_ch = generateRefContigChunks()

  // call variants
  bam2ind_ch = Channel.fromPath(params.bam2indpath)
  ped_ch = Channel.fromPath(params.pedpath)
  parentcallchunks_ch = callVarsAndParents(contigchunks_ch, bam2ind_ch, ped_ch)

  // TODO: ADD filtering! (but do not put to output dir)

  // estimate pi
  agp_ch = Channel.fromPath(params.agpfile)
  bedcallchunks_ch = calcNucDiversity(parentcallchunks_ch, bam2ind_ch, agp_ch)

  intervals_ch = Channel.fromPath(params.genomeintervals)
  pi_ch = collectNucDiversities(bedcallchunks_ch.collect(), intervals_ch)
   
}
