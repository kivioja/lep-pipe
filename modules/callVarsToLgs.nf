#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
  do the variant calling and separate to chromosomes
*/


/*
  Generate contig chunks as bed files for samtools to run variant calling in parallel
*/
process generateRefContigChunks {

  output:
  path 'contigchunk_*'

  // need to be a bed files 0..length 
  script:
  """
  cat $params.refgenome | awk -f ${params.lepanchordir}/contigLength.awk | shuf >contigs.shuf.length
  perl -ple 's/\\t/\\t0\\t/' contigs.shuf.length >contigs.bed
  split -l $params.contigchunksize contigs.bed contigchunk_
  """
}


/*
  Do both variant and parent calling
  TODO: output and collect also the variant calls
*/  
process callVarsAndParents {
 
  input:
  each path(contigchunk)
  path bam2ind
  path ped

  output:
  path 'data.chunk.call'

  script:
  if (params.xlimit == "inf" && params.zlimit == "inf" ) { 
    """
    cut -f 1 -d ',' $bam2ind | \
      perl -ple 's{^}{${params.bamdir}/}' >bamlist.txt
      cut -f 2 -d ',' $bam2ind  > mapping.txt
    samtools mpileup -l $contigchunk -q 10 -Q 10 -s --bam-list bamlist.txt | \
      java ${params.jvmoptions} -cp ${params.lepmapdir} Pileup2Likelihoods \
      numLowerCoverage=${params.LepMap_numLowerCoverage} \
      minAlleleFreq=${params.LepMap_minAlleleFreq} \
      mappingFile=mapping.txt | gzip  > chunkpost
    zcat chunkpost | java ${params.jvmoptions} -cp ${params.lepmapdir} ParentCall2 \
	    data=${ped} \
	    posteriorFile=- \
      halfSibs=${params.halfsibs} \
      removeNonInformative=${params.removeNonInformative} | gzip >data.chunk.call
    #perl /scratch/project_2005425/TEST/check_errors.pl < .command.err
    """
  } else if (params.xlimit != "inf") {
    """
    cut -f 1 -d ',' $bam2ind | \
      perl -ple 's{^}{${params.bamdir}/}' >bamlist.txt
      cut -f 2 -d ',' $bam2ind  > mapping.txt
    samtools mpileup -l $contigchunk -q 10 -Q 10 -s --bam-list bamlist.txt | \
      java ${params.jvmoptions} -cp ${params.lepmapdir} Pileup2Likelihoods \
      numLowerCoverage=${params.LepMap_numLowerCoverage} \
      minAlleleFreq=${params.LepMap_minAlleleFreq} \
      mappingFile=mapping.txt | gzip  > chunkpost
    zcat chunkpost | java ${params.jvmoptions} -cp ${params.lepmapdir} ParentCall2 \
	    data=${ped} \
	    posteriorFile=- \
      XLimit=${params.xlimit} \
      halfSibs=${params.halfsibs} \
      removeNonInformative=${params.removeNonInformative} | gzip >data.chunk.call
    """
  } else if (params.zlimit != "inf") {
    """
    cut -f 1 -d ',' $bam2ind | \
      perl -ple 's{^}{${params.bamdir}/}' >bamlist.txt
      cut -f 2 -d ',' $bam2ind  > mapping.txt
    samtools mpileup -l $contigchunk -q 10 -Q 10 -s --bam-list bamlist.txt | \
      java ${params.jvmoptions} -cp ${params.lepmapdir} Pileup2Likelihoods \
      numLowerCoverage=${params.LepMap_numLowerCoverage} \
      minAlleleFreq=${params.LepMap_minAlleleFreq} \
      mappingFile=mapping.txt | gzip  > chunkpost
    zcat chunkpost | java ${params.jvmoptions} -cp ${params.lepmapdir} ParentCall2 \
	    data=${ped} \
	    posteriorFile=- \
      ZLimit=${params.zlimit} \
      halfSibs=${params.halfsibs} \
      removeNonInformative=${params.removeNonInformative} | gzip >data.chunk.call
    """
  } else { 
      error "Could not figure out X/Z limits"
  }    
}


/*
 * collect parent call files into one file compressed file
 */
process collectCalls {

  publishDir params.mapdir, mode: 'copy', overwrite: true

  input:
  file 'data.first.chunk'
  file 'data.chunk*' 

  output:
  file 'data.call.gz'

  // need to take the first 7 lines from the very first line and all other lines from all files
  // - relies on the gzip property that allows to directly concatenate compressed files
  script:
  """
  zcat data.first.chunk | head -n 7 | gzip >data.call.gz
  for chunk in data.chunk*; do
    zcat \$chunk | tail -n +8 | gzip >>data.call.gz
  done
  """
}


/*
  Lep-MAP filtering + removing markers in consecutive contig positions (leaving only the first one of the run)
*/
process	filterCalls {

  publishDir params.mapdir, mode: 'copy', overwrite: true

  input:
  file 'data.call.gz'

  output:
  file 'data.call.filt.gz'

  script:
  """
  zcat data.call.gz | java ${params.jvmoptions} -cp ${params.lepmapdir} Filtering2 dataTolerance=${params.datatolerance} data=- | \
    filter_consecutive_pos_markers.pl | gzip >data.call.filt.gz
  """

}


/*
  separation into linkage groups / (pseudo)chromosomes
  TODO add informative mask options
*/
process	separateChromosomes {

  publishDir params.mapdir, mode: 'copy', overwrite: true

  input:
  file callfile
  each lodlim

  output:
  file "map_${lodlim}.txt"

  script:
  """
  zcat $callfile | java ${params.jvmoptions} -Xmx${params.javaheapsize} -cp ${params.lepmapdir} SeparateChromosomes2 data=- \
    samplePairs=${params.samplepairsfraction} \
    informativeMask=${params.informativemask} \
    lod3Mode=${params.LepMap_lod3Mode} \
    distortionLod=1 lodLimit=${lodlim} numThreads=${params.numthreads} >map_${lodlim}.txt
  """

}




/* 
  generate two auxillarily files needed for downstream analysis
*/
process genSeqNameAndSNPfiles {

  publishDir params.mapdir, mode: 'copy', overwrite: true

  input:
  path callsfile

  output:
  file 'snps.txt'
  file 'genome_seq_names.txt'
     
  shell:
  """
  zcat data.call.filt.gz | tail -n +7 | cut -f 1,2 >snps.txt
  seqkit seq -n ${params.refgenome} >genome_seq_names.txt
  """
}










