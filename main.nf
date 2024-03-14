#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 MAIN workflow for generating linkage maps
 assumes that bams already exist (see seq_fw to create those)
*/

// set default parameter values
// (needs to be before includes: https://training.nextflow.io/basic_training/modules/#parameter-scopes)
params.javaheapsize = "64G"
setParamDefaults()

//include processes 
include { generateRefContigChunks; callVarsAndParents; collectCalls; filterCalls;  separateChromosomes; genSeqNameAndSNPfiles } from './modules/callVarsToLgs.nf'
include { splitChromosomes } from './modules/refineLgs.nf'
include { removeSmallAndJoinSingles; orderMarkers; storeFinalOrdering } from './modules/orderLgsMarkers'
include { runAnchorWrapper; fixContigs; liftOverMaps; fitStepAndEstimateRecomb } from './modules/anchorAndEstimateRecomb.nf'

// which steps to run
// can start from either from beginning or from runorder (if first 2 steps already done)
params.runvars2lgs = true
params.runsplitlgs = true
params.runorder = true
params.runanchor = true

// store the parameter values used
paramValFile = file("RUN_PARAMS_" + params.species + ".tsv")
paramValFile.write "species\tparam\tvalue\n"
params.each { paramValFile.append("$params.species\t$it.key\t$it.value\n") }

/*
 MAIN workflow
 */
workflow {

    if (params.runvars2lgs == true) {
        var_chs = vars2lgs_wf()
    }
    if (params.runsplitlgs == true) {
        refinedmap_ch = split_lgs_wf(var_chs[0], var_chs[1])
    }
    if (params.runorder == true) {
        // if previous step not run assume that the files are already available
        if (params.runsplitlgs == false) {
            filteredparentcalls_ch = Channel.fromPath(params.mapdir + "/" + "data.call.filt.gz")  
            refinedmap_ch = Channel.fromPath(params.mapdir + "/" + params.refinedmap)
        } else {
            filteredparentcalls_ch = var_chs[0]
        }
        finalorders_chs = order_wf(filteredparentcalls_ch, refinedmap_ch)
    }
    if (params.runanchor == true) {
        anchor_wf(finalorders_chs[0], finalorders_chs[1])
    }

}



/*
 * Call variants, parents and then group to linkage groups
 */
workflow vars2lgs_wf {

    main:
        // divide the contings into junks to be processed in parallel
        contigchunks_ch = generateRefContigChunks()

        // Do variant and parent calling as well as filtering in one process 
        bam2ind_ch = Channel.fromPath(params.bam2indpath)
        ped_ch = Channel.fromPath(params.pedpath)
        parentcallchunks_ch = callVarsAndParents(contigchunks_ch, bam2ind_ch, ped_ch)

        // collect calls to one 
        allparentcalls_ch = collectCalls(parentcallchunks_ch.first(), parentcallchunks_ch.collect())

        // filter variant calls
        filteredparentcalls_ch = filterCalls(allparentcalls_ch)

        // Run linkage grouping with different lod values
        def lodlimits = params.minlod..params.maxlod
        seplgs_ch = separateChromosomes(filteredparentcalls_ch, lodlimits)

        genSeqNameAndSNPfiles(filteredparentcalls_ch)

    emit:
        filteredparentcalls_ch
        seplgs_ch
}


/*
 * after choosing the LOD limit, split chromosomes if needed
 */
workflow split_lgs_wf {
    take:
        filteredparentcalls_ch
        seplgs_ch

    main:    
        initmap_ch = seplgs_ch
                        .first { it.name == "map_" + params.LepMap_lodLimit + ".txt" }
        refinedmap_ch = splitChromosomes(filteredparentcalls_ch, initmap_ch)
    
    emit:
        refinedmap_ch
}


/*
 * split LG:s if needed and order chromosomes 
 */
workflow order_wf {
    take:
        filteredparentcalls_ch 
        refinedmap_ch

    main:    
        jsmap_ch = removeSmallAndJoinSingles(filteredparentcalls_ch, refinedmap_ch)

        // order markers, if several attempts for each LG, choose the best
        def chroms = 1..params.chromnum
        def attempts = 1..3

        orderattempts_ch = orderMarkers(filteredparentcalls_ch, jsmap_ch[0], chroms, attempts)
        finalorders_chs = storeFinalOrdering(filteredparentcalls_ch, orderattempts_ch[0].collect(), orderattempts_ch[1].collect())
    
    emit:
        finalorders_chs[2]
        finalorders_chs[3]
}


/*
 *  anchor and fix genome, estimate recombination rates
 */
workflow anchor_wf {
    take:
        mapforla_ch
        mapforstep_ch

    main:
        def chroms = 1..params.chromnum

        // to get soft link "bin" directory pointing to the java classes for anchor wrapper 
        anchorbin_ch = Channel.fromPath(params.lepanchordir + "/" + "bin")
        anchoroutput_chs = runAnchorWrapper(anchorbin_ch, mapforla_ch)

        // fix contigs
        fixedcontigs_chs = fixContigs(mapforla_ch, anchoroutput_chs[3], chroms)

        // lift over to LG coordinates
        liftedmap_chs = liftOverMaps(fixedcontigs_chs[0].collect(), mapforla_ch, mapforstep_ch)

        // fit step function and estimate recombination rates
        estimate_chs = fitStepAndEstimateRecomb(liftedmap_chs[1], liftedmap_chs[2], chroms)
    
    emit:
        estimate_chs[0] 
}




/*
 * set the default parameter values  
 */
def setParamDefaults() {

    params.numthreads = 4
    params.jvmoptions = "-XX:+PerfDisableSharedMem"

    // default file names
    params.mapping_bam2ind_file = params.species + "_bam2indname.csv"
    params.bam2indpath = params.sampledir + "/" + params.mapping_bam2ind_file
    params.pedfile =  params.species + "_ped.tsv"
    params.pedpath = params.sampledir + "/" + params.pedfile

    params.mappingfile = 'mapping.txt'

    params.minlod = 10
    params.maxlod = 30

    params.refinedmap = "map_refined.txt"

    //
    // callVarsToLgs.nf
    //
    params.bamtemplate = '/something/*.bam'
    params.bamdir = 'somepath/'
    params.mapdir = 'somepath/'

    params.sampleibd = false
    params.LepMap_numLowerCoverage = 0.3
    params.LepMap_minAlleleFreq = 0.1
    params.LepMap_lod3Mode = 1
    params.contigchunksize = 10
    params.samplepairsfraction = 1
    params.xlimit = "inf"
    params.zlimit = "inf"
    params.halfsibs = 0
    params.removeNonInformative = 1
    params.informativemask = "0123"
    params.datatolerance = 0.001 

    //
    // refineLgs.nf
    // - no splitting by default
    params.splitchromosomes = "()"
    params.splitlimits = "()"

    //
    // orderLgsMarkers  
    //
    params.joinsingles = true
    params.scale="M/N 2"
    params.proximityscale=100

    // if specific LOD limit not given for joining singles, 
    // use the same as for separating chrom
    params.joinsingleslimit = params.LepMap_lodLimit

    // set true for species like some butteflies such that there is no 
    // recombination in females
    params.nofemalerecomb = false

    //
    // anchorAndEstimateRecomb
    //
    params.fixcontigs = false
    params.LepAnchor_minImprovement = 1 

    params.makechain = false
    params.makepaf = false

    params.anchordir = params.mapdir + "/" + "anchoring"
    params.anchorchaindir = params.mapdir + "/" + "anchoring_chained"
    params.anchorpafdir = params.mapdir + "/" + "anchoring_with_paf"
    params.recdensitydir = params.mapdir + "/" + "densities"

    lepanchorbinpath_ch = Channel.fromPath(params.lepanchordir + "/" + "bin")

    // defaults for estimation 
    params.gridsize = 1024
    params.numsamples = 5

    params.minsupport = 2
    params.minlength = 10000 

    params.endsupport = params.minsupport
    params.bandwidth = 0

}
