#!/usr/bin/env nextflow
nextflow.enable.dsl=2


//defaults
params.numthreads = 4
params.LepMap_lod3Mode = 3
params.joinsingles = true

params.scale="M/N 2"

// if specific LOD limit not given for joining singles, 
// use the same as for separating chrh
params.joinsingleslimit = params.LepMap_lodLimit


/*
 * first remove small ones and join singles to the chosen map - also reorders 
 */
process removeSmallAndJoinSingles {

  publishDir params.mapdir, mode: 'copy', overwrite: true

  input:
  path callsfile
  path refinedmap

  output:
  path "map_js.txt"
  path "map_js_lgs_only.txt"

  // note: the joinSingles2All output file can contain variable number columns 
  script:
  if (params.joinsingles == true)
    """
    zcat $callsfile | java -Xmx${params.javaheapsize} -cp ${params.lepmapdir} SeparateChromosomes2 data=- \
      map=${refinedmap} \
      distortionLod=1 lodLimit=${params.LepMap_lodLimit} \
      lod3Mode=${params.LepMap_lod3Mode} \
      sizeLimit=${params.LepMap_sizeLimit} >map_nosmall.txt
    zcat $callsfile | java -Xmx${params.javaheapsize} -cp ${params.lepmapdir} JoinSingles2All data=- \
      map=map_nosmall.txt distortionLod=1 lodLimit=${params.joinsingleslimit} \
      iterate=1 lodDifference=2 \
      numThreads=${params.numthreads} >map_js.txt
    cut -f 1 map_js.txt >map_js_lgs_only.txt 
    """
  else 
    """
    zcat $callsfile | java -Xmx${params.javaheapsize} -cp ${params.lepmapdir} SeparateChromosomes2 data=- \
      map=${refinedmap} \
      distortionLod=1 lodLimit=${params.LepMap_lodLimit} \
      lod3Mode=${params.LepMap_lod3Mode} \
      sizeLimit=${params.LepMap_sizeLimit} >map_js.txt
    cut -f 1 map_js.txt >map_js_lgs_only.txt 
    """

}



/*
 * order each chromosome independently several times 
 */
process orderMarkers {
 
  input:
  path callsfile
  path finalmap 
  each chrom 
  each attempt 

  output:
  path "order_${chrom}_${attempt}.txt"
  path "order_${chrom}_${attempt}.int"

  script:
  """
  zcat $callsfile | java -Xmx${params.javaheapsize} -cp ${params.lepmapdir} OrderMarkers2 data=- \
    useMorgan=1 \
    scale=${params.scale} \
    map=${finalmap} chromosome=$chrom calculateIntervals=order_${chrom}_${attempt}.int \
    numThreads=${params.numthreads} >order_${chrom}_${attempt}.txt
  """
}



/*
 * pick the highhest likelihood chromosome orderings
 */
process chooseHighestLikelihoodOrdering {
  publishDir params.mapdir + "/" + "initorders", mode: 'copy', overwrite: true 

  input:
  path 'order_attempt*.txt' 
  path 'order_attempt*.int' 

  output: 
  path 'order*.txt'
  path 'order*.int'

  script:
  """
  ls -1 order_attempt*.txt >files_list.txt
  collect_highest_likelihood_orderings.pl < files_list.txt
  """
}



/*
 * Filter end markers
 */
process filterEndMarkers { 

  publishDir params.mapdir, mode: 'copy', overwrite: true 

  input:
  path callsfile
  path finalmap
  path 'order*.int'

  output:
  path 'removed_by_endfiltering.txt'

  shell:
  """
  zcat data.call.filt.gz | tail -n +7 | cut -f 1,2 >snps.txt
  for X in {1..!{params.chromnum}}; do \
    awk -vn=\$X '(NR==FNR){map[NR-1]=\$0}(NR!=FNR){\$1=map[\$1] "\t"n;print}' \
      snps.txt order\$X.int; done >map_all.txt
  gawk -f !{params.scriptdir}/filterEnds.awk map_all.txt >removed_by_endfiltering.txt 
  """

}


/*
 * Rerun the ordering
 */
process reorderMarkersWithoutEndFiltered {

  input:
  path callsfile
  path finalmap 
  path 'removed_by_endfiltering.txt'
  each chrom 
  each attempt 

  output:
  path "order_${chrom}_${attempt}.txt"
  path "order_${chrom}_${attempt}.int"

  script:
  """
  zcat $callsfile | java -Xmx${params.javaheapsize} -cp ${params.lepmapdir} OrderMarkers2 data=- \
    removeSNPs=file:removed_by_endfiltering.txt \
    useMorgan=1 \
    scale=${params.scale} \
    map=${finalmap} chromosome=$chrom calculateIntervals=order_${chrom}_${attempt}.int \
    numThreads=${params.numthreads} >order_${chrom}_${attempt}.txt
  """
 
}


/*
 * pick the final highhest likelihood chromosome orderings and store
 * - collect also to the whole-genome map files, needed if Lep-Anchor is run 
 */
process storeFinalOrdering {
  publishDir params.mapdir + "/" + "orders", mode: 'copy', overwrite: true 

  input:
  path 'data.call.filt.gz'
  path 'order_attempt*.txt' 
  path 'order_attempt*.int' 

  output: 
  path 'order*.txt'
  path 'order*.int'
  path 'map_for_la.txt'
  path 'map_for_step.txt'

  // shell needed to differentiate between bash/awk and groovy variables 
  shell:
  """
  ls -1 order_attempt*.txt >files_list.txt
  collect_highest_likelihood_orderings.pl < files_list.txt
  zcat data.call.filt.gz | tail -n +7 | cut -f 1,2 >snps.txt
  cat !{params.refgenome} | awk -f !{params.lepanchordir}/contigLength.awk >contigs.length
  for X in {1..!{params.chromnum}}; do \
    awk -vn=\$X '(NR==FNR){map[NR-1]=\$0}(NR!=FNR){\$1=map[\$1] "\t"n;print}' \
      snps.txt order\$X.int; done >map_for_la.txt
  for X in {1..!{params.chromnum}}; do \
    awk -vFS="\t" -vOFS="\t" -vn=\$X '(NR==FNR){map[NR-1]=\$0}(NR!=FNR && \$1 ~ /^[^#]/){\$1=map[\$1] "\t"n;print}' \
      snps.txt order\$X.txt; done >map_for_step.txt      
  """
}

