#!/usr/bin/env nextflow
nextflow.enable.dsl=2


// default 
params.javaheapsize = "64G"

/*
  include processes 
*/
include { generateRefContigChunks; callVarsAndParents; collectCalls; filterCalls;  separateChromosomes; genSeqNameAndSNPfiles } from './modules/callVarsToLgs.nf'
include { splitChromosomes } from './modules/refineLgs.nf'
include { removeSmallAndJoinSingles; orderMarkers; chooseHighestLikelihoodOrdering; filterEndMarkers; reorderMarkersWithoutEndFiltered; storeFinalOrdering } from './modules/orderLgsMarkers'
include { runAnchorWrapper; fixContigs; liftOverMaps; fitStepAndEstimateRecomb } from './modules/anchorAndEstimateRecomb.nf'

// default file names
params.mapping_bam2ind_file = params.species + "_bam2indname.csv"
params.bam2indpath = params.sampledir + "/" + params.mapping_bam2ind_file
params.pedfile =  params.species + "_ped.tsv"
params.pedpath = params.sampledir + "/" + params.pedfile

params.mappingfile = 'mapping.txt'

params.minlod = 10
params.maxlod = 30

params.refinedmap = "map_refined.txt"


  
params.runvars2lgs = true
params.runsplitlgs = true
params.runorder = true
params.runanchor = true

/*
 MAIN workflow
 */
workflow {
    if (params.runvars2lgs == true) {
        var_chs = vars2lgs_wf()
    }
    if (params.runsplitlgs == true) {
        refinedmap_ch = split_lgs_wf(var_chs[0], var_chs[1])
    }ß
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
 * split LG:s if needed order, filter ends, reorder chromosomes 
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
        firstorders_ch = chooseHighestLikelihoodOrdering(orderattempts_ch[0].collect(), orderattempts_ch[1].collect())

        // filter ends 
        toberemoved_ch = filterEndMarkers(filteredparentcalls_ch, jsmap_ch[0], firstorders_ch[1].collect())

        // reorder without filtered end markers 
        reorderattempts_ch = reorderMarkersWithoutEndFiltered(filteredparentcalls_ch, jsmap_ch[0], toberemoved_ch, chroms, attempts)
        finalorders_chs = storeFinalOrdering(filteredparentcalls_ch, reorderattempts_ch[0].collect(), reorderattempts_ch[1].collect())
    
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
