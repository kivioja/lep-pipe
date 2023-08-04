#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
  do the variant calling and use variants to check the parents 
*/


params.numthreads = 4

// other defaults
params.allparentchildpairs = false
params.sampleibdfraction = 0.1


/*
    Calls variants from a subset of positions
 */
process callVarsSubset {
 
    input:
    path bam2ind

    output:
    path 'post.gz'

    script:
    """
    cut -f 1 -d ',' $bam2ind | \
        perl -ple 's{^}{${params.bamdir}/}' >bamlist.txt
        cut -f 2 -d ',' $bam2ind  > mapping.txt
    samtools mpileup -q 10 -Q 10 -s --bam-list bamlist.txt | \
    perl -ne 'print if (rand() < ${params.sampleibdfraction})' | \
    java -cp ${params.lepmapdir} Pileup2Likelihoods \
        mappingFile=mapping.txt | gzip >post.gz 
    """    
}



/* 
 * IBD calculation, should take sample (10%) if large number of markers (> 100,000) 
 */
process calcRelatedness {
  
  publishDir params.mapdir, mode: 'copy', overwrite: true

  input:
  file 'post.gz' 

  output:
  file 'ibd.txt'

  script:
  """
  zcat post.gz | java -Xmx${params.javaheapsize} -cp ${params.lepmapdir} IBD posteriorFile=- \
      numThreads=${params.numthreads} > ibd.txt
  """

}


/* 
 * Calculate Mendel errors for all parent/child pairs
 */
process calcParentMendelErrors {
  
    publishDir params.mapdir, mode: 'copy', overwrite: true

    input:
    path 'post.gz'
    path 'parents.txt'

    output:
    file 'parent_mendel_errors.txt'

    script:
    if (params.allparentchildpairs==true)
      """
      zcat post.gz | java -Xmx${params.javaheapsize} -cp ${params.lepmapdir} IBD numThreads=${params.numthreads} \
      parents=file:parents.txt allParentChildPairs=1 posteriorFile=- > parent_mendel_errors.txt
      """
    else
      """
      zcat post.gz | java -Xmx${params.javaheapsize} -cp ${params.lepmapdir} IBD numThreads=${params.numthreads} \
      parents=file:parents.txt allParentChildPairs=1 posteriorFile=- > parent_mendel_errors.txt
      """
}









