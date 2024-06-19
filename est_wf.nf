nextflow.enable.dsl=2

/*
  Do only the Lep-Anchor fit and estimation
  - fit step function and estimate: linkage map given as input or, for example different parameters than standard
  - evaluate in physical order instead of fitting step function
  - estimate gamma model: fit given as input
 */

 // default values, put before includes
params.javaheapsize = "64G"

// defaults, Lep-MAP part
params.snpfile = params.mapdir + "/" + 'snps.txt'
params.callsfile =  params.mapdir + "/" + 'data.call.filt.gz'
params.finalmapfile = params.mapdir + "/" + 'map_js.txt'

params.scale="M/N 2"
params.proximityscale=100

// defaults, correspond to Lep-Anchor + fit output
params.inputdir = params.mapdir + "/" + "densities"
params.mapforlafile = params.inputdir + "/" + "map_for_la.liftover.txt"
params.mapforstepfile = params.inputdir + "/" + "map_for_step.liftover.txt" 
// fit and estimate output and gamma model input directory
params.recdensitydir = params.inputdir
// gamma model outputdir, change
params.gammamodeldir = params.inputdir

params.gridsize = 1024
params.numsamples = 5

params.minsupport = 2
params.minlength = 10000 

params.endsupport = params.minsupport
params.bandwidth = 0

params.evaluateorder = false
params.evaluatereforder = false
params.newfit = false
params.reffit = false
params.gammamodel = 0
params.nofemalerecomb = false

params.lg2refchrom = ""
params.refseq2refchrom = ""

/*
  include processes 
*/
include { fitStepAndEstimateRecomb; fitStepReferenceAndEstimateRecomb; evaluateOrderAndEstimateRecomb; evaluateReferenceOrderAndEstimateRecomb; estimateGammaModel } from './modules/anchorAndEstimateRecomb.nf'
  

workflow {

  def chroms = 1..params.chromnum  

  
  // Evaluate map in the physical order, uses Lep-Anchor output as default
  if (params.evaluateorder == true) {

    snps_ch = Channel.fromPath(params.snpfile)
    calls_ch =  Channel.fromPath(params.callsfile)
    map_ch = Channel.fromPath(params.finalmapfile)

    orderevals_ch = evaluateOrderAndEstimateRecomb(snps_ch, calls_ch, map_ch, chroms)
  }

  // Evaluate map in physical order defined by the reference, use its chromosome names
  if (params.evaluatereforder == true) {

    // read the mapping from Lep-Map linkage group numbers to reference chromosome numbers
    // note: each does not handle tuples, so creating one combined input channel
    lg2chrom_ch =  Channel.fromPath(params.lg2refchrom)
                          .splitCsv(sep: "\t")
                          .map { row-> tuple(row[0], row[1]) }
   
    snps_ch = Channel.fromPath(params.snpfile)
    calls_ch =  Channel.fromPath(params.callsfile)
    map_ch = Channel.fromPath(params.finalmapfile)

    ch2 = lg2chrom_ch.combine(snps_ch)
    ch3 = ch2.combine(calls_ch)
    combined_ch = ch3.combine(map_ch)

    //combined_ch.view()
    orderevals_ch = evaluateReferenceOrderAndEstimateRecomb(combined_ch)
  }


  // Fit step function and estimate, uses linkage map names 
  if (params.newfit == true) {

    mapforla_ch = Channel.fromPath(params.mapforlafile)
    mapforstep_ch = Channel.fromPath(params.mapforstepfile)   
 
    estimate_chs = fitStepAndEstimateRecomb(mapforla_ch, mapforstep_ch, chroms)
  }


  // Fit step function and estimate, use reference genome names
  if (params.reffit == true) {

    // read the mapping from Lep-Map linkage group numbers to reference chromosome numbers
    // note: each does not handle tuples, so creating one combined input channel
    lg2chrom_ch =  Channel.fromPath(params.lg2refchrom)
                          .splitCsv(sep: "\t")
                          .map { row-> tuple(row[0], row[1]) }
  
    mapforla_ch = Channel.fromPath(params.mapforlafile)
    mapforstep_ch = Channel.fromPath(params.mapforstepfile)   

    ch2 = lg2chrom_ch.combine(mapforla_ch)
    combined_ch = ch2.combine(mapforstep_ch)

    estimate_chs = fitStepReferenceAndEstimateRecomb(combined_ch)
  }


  if (params.gammamodel != 0) {

    fitpath = params.recdensitydir + "/" + "fit*.txt"
    fits_ch = Channel.fromPath(fitpath).collect()

    gamma_chs = estimateGammaModel(chroms, fits_ch)
  }

   
}

