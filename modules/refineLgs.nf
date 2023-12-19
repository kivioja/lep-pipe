#!/usr/bin/env nextflow

/*
  run LepMAP after deciding which raw map to use as starting (from different LOD score limits)
  example
*/



/*
  take the given map and split the given chromosomes
  - remember that linkage group numbering changes unless specified not to 
*/
process splitChromosomes {

  publishDir params.mapdir, mode: 'copy', overwrite: true

  // these need to be processed in so that the previous output is the next input
  input:
  path callsfile 
  path 'map_split_in.txt'

  output:
  path 'map_refined.txt'

  // note that sc here is a bash variable, not groovy, and thus needs to be escaped 
  script:
  if (params.splitchromosomes != '()')
    """
    chroms=${params.splitchromosomes}
    lims=${params.splitlimits}
    numtosplit=\${#chroms[@]}
    for i in \$(seq 1 \$numtosplit)
    do
      sc=\${chroms[i-1]}
      lim=\${lims[i-1]}
      echo "splitting chromosome \$sc using lod limit \$lim"
      zcat $callsfile | java -Xmx${params.javaheapsize} -cp ${params.lepmapdir} SeparateChromosomes2 data=- \
        map=map_split_in.txt \
        distortionLod=1 \
        lodLimit=\$lim \
        lg=\$sc \
        renameLGs=0 \
        samplePairs=${params.samplepairsfraction} \
        numThreads=${params.numthreads} \
        >map_split_out.txt
      cp map_split_out.txt map_split_in.txt
    done
    cp map_split_out.txt map_refined.txt
    """
  else 
    """
    cp map_split_in.txt map_refined.txt
    """
}

 

