#!/usr/bin/env nextflow
nextflow.enable.dsl=2



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
    zcat $callsfile | java ${params.jvmoptions} -Xmx${params.javaheapsize} -cp ${params.lepmapdir} SeparateChromosomes2 data=- \
      map=${refinedmap} \
      distortionLod=1 lodLimit=${params.LepMap_lodLimit} \
      lod3Mode=${params.LepMap_lod3Mode} \
      sizeLimit=${params.LepMap_sizeLimit} >map_nosmall.txt
    zcat $callsfile | java ${params.jvmoptions} -Xmx${params.javaheapsize} -cp ${params.lepmapdir} JoinSingles2All data=- \
      map=map_nosmall.txt distortionLod=1 lodLimit=${params.joinsingleslimit} \
      iterate=1 lodDifference=2 \
      numThreads=${params.numthreads} >map_js.txt
    cut -f 1 map_js.txt >map_js_lgs_only.txt 
    """
  else 
    """
    zcat $callsfile | java ${params.jvmoptions} -Xmx${params.javaheapsize} -cp ${params.lepmapdir} SeparateChromosomes2 data=- \
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
  if (params.nofemalerecomb == false) {
    """
    zcat $callsfile | java ${params.jvmoptions} -Xmx${params.javaheapsize} -cp ${params.lepmapdir} OrderMarkers2 data=- \
      useMorgan=1 \
      scale=${params.scale} \
      proximityScale=${params.proximityscale} \
      map=${finalmap} chromosome=$chrom calculateIntervals=order_${chrom}_${attempt}.int \
      numThreads=${params.numthreads} >order_${chrom}_${attempt}.txt
    """
  } else {
    """
    zcat $callsfile | java ${params.jvmoptions} -Xmx${params.javaheapsize} -cp ${params.lepmapdir} OrderMarkers2 data=- \
      useMorgan=1 \
      scale=${params.scale} \
      proximityScale=${params.proximityscale} \
      recombination2=0 \
      map=${finalmap} chromosome=$chrom calculateIntervals=order_${chrom}_${attempt}.int \
      numThreads=${params.numthreads} >order_${chrom}_${attempt}.txt
    """
  }
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





