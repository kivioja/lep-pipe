#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * various quality control processes to check for example the relatedness of individuals
 */

 // default 
params.javaheapsize = "64G"

/*
  include processes 
*/
include { callVarsSubset; calcRelatedness; calcParentMendelErrors } from './modules/calcRelatedness.nf'

params.runrelatedness = true;
params.runmendelerrors = true;

// default file names
params.mapping_bam2ind_file = params.species + "_bam2indname.csv"
params.bam2indpath = params.sampledir + "/" + params.mapping_bam2ind_file
params.parent_file = params.species + "_parents.txt"
params.parentpath = params.sampledir + "/" + params.parent_file

workflow {

    bam2ind_ch = Channel.fromPath(params.bam2indpath)

    post_ch = callVarsSubset(bam2ind_ch)

    if (params.runrelatedness == true) {
        calcRelatedness(post_ch)
    }

    if (params.runmendelerrors == true) {
        parents_ch = Channel.fromPath(params.parentpath)
        calcParentMendelErrors(post_ch, parents_ch)
    }

}
