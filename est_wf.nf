nextflow.enable.dsl=2

/*
  Do only the Lep-Anchor fit and estimation
  - fit step function and estimate: linkage map given as input or, for example different parameters than standard
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
params.newfit = false
params.gammamodel = 0


/*
  include processes 
*/
include { fitStepAndEstimateRecomb; evaluateOrderAndEstimateRecomb; estimateGammaModel } from './modules/anchorAndEstimateRecomb.nf'
  

workflow {

  def chroms = 1..params.chromnum  

  if (params.evaluateorder == true) {

    snps_ch = Channel.fromPath(params.snpfile)
    calls_ch =  Channel.fromPath(params.callsfile)
    map_ch = Channel.fromPath(params.finalmapfile)

    orderevals_ch = evaluateOrderAndEstimateRecomb(snps_ch, calls_ch, map_ch, chroms)
  }


  if (params.newfit == true) {

    mapforla_ch = Channel.fromPath(params.mapforlafile)
    mapforstep_ch = Channel.fromPath(params.mapforstepfile)   
 
    estimate_chs = fitStepAndEstimateRecomb(mapforla_ch, mapforstep_ch, chroms)
  }


  if (params.gammamodel != 0) {

    fitpath = params.recdensitydir + "/" + "fit*.txt"
    fits_ch = Channel.fromPath(fitpath).collect()

    gamma_chs = estimateGammaModel(chroms, fits_ch)
  }

   
}

