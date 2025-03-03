nextflow.enable.dsl=2

/*
 * Postprocessing of result files
 *  calculate sequence features like parents variant positions
 */

 // default 
params.javaheapsize = "64G"

params.halfsibs = 0
params.LepMap_numLowerCoverage = 0.3
params.sampledir = ""

params.makewindows = true
params.getinterference = false
params.runvars = false
params.runoldcounts = false

/*
  include processes 
*/
include { generateRefContigChunks; callVarsAndParents } from './modules/callVarsToLgs.nf'
include { makeGenomeWindows; countWindowSeqFeatures; calcWindowFeatures; getWindowInterference; callToBedAllParentVars; collectParentVariantWindowCounts } from './modules/calcFeatures.nf'
  
// default file names
// - pi calculated from the parents but children aid in parent variant calling
params.mapping_bam2ind_file = params.species + "_bam2indname.csv"
params.bam2indpath = params.sampledir + "/" + params.mapping_bam2ind_file
params.pedfile =  params.species + "_ped.tsv"
params.pedpath = params.sampledir + "/" + params.pedfile
params.agpfile = params.mapdir + "/" + "densities" + "/" + params.species + "_LAfixed.agp"
params.fixedgenome = params.mapdir + "/" + "densities" + "/" + params.species + "_LAfixed.fa"
params.estfiletemplate = params.mapdir + "/densities/" + "estimates"

workflow {

  params.windowsize = 8

  bedoutfile = params.species + "_marey_" + params.windowsize + "_windows.bed"
  recoutfile = params.species + "_marey_" + params.windowsize + "_windows.tsv"

  if (params.makewindows == true) {
    windows_ch = makeGenomeWindows(params.estfiletemplate, bedoutfile, recoutfile, params.windowsize)
    windowsfile = windows_ch[0]
  } else {
    windowsfile = params.genomewindowfile
  }

  // basic sequence features like GC-content
  featoutcountfile =  params.species + "_seq_feat_counts_" + params.windowsize + "_windows.tsv"
  // bending scores
  bendscorefile =  params.species + "_bending_scores_" + params.windowsize + "_windows.tsv"

  // old, only basic 
  if (params.runoldcounts == true) {
    feats_ch = countWindowSeqFeatures(windowsfile, params.fixedgenome, featoutcountfile)
  }
  // new, bending too
  feats_ch = calcWindowFeatures(windowsfile, params.fixedgenome, 
                                featoutcountfile, bendscorefile)

  // run window interference script for each chromosome estimate file
  if (params.getinterference == true) {
    def chroms = 1..params.chromnum
    interference_ch = getWindowInterference(params.estfiletemplate, chroms, params.windowsize)
  }

  if (params.runvars == true) {
    // pi estimation needs also non-variant positions, thus    
    // divide the contings into junks to be processed in parallel
    contigchunks_ch = generateRefContigChunks()

    // call variants
    bam2ind_ch = Channel.fromPath(params.bam2indpath)
    ped_ch = Channel.fromPath(params.pedpath)
    agp_ch = Channel.fromPath(params.agpfile)
    parentcallchunks_ch = callToBedAllParentVars(contigchunks_ch, bam2ind_ch, ped_ch, agp_ch)

    varcountoutfile = params.species + "_parent_var_counts_" + params.windowsize + "_windows.tsv"
    cout_ch = collectParentVariantWindowCounts(parentcallchunks_ch.collect(), windowsfile, varcountoutfile)
  }
}

