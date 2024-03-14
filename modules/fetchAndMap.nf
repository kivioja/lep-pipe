#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
  get fastq and map to reference genome, for large WGS sets
  - does not store fastqs to save space (if fetched here)
  - uses bwa-mem2 to save time
*/



params.numthreads = 4

// bam stats
statdir = params.bamdir + "/" + "stats"
params.samplebigbams = false


process fetchAndMapToRefGenomeSingle {

  publishDir params.bamdir, mode: 'move', overwrite: true
 
  input:
  val id

  output:
  path "${id}.bam"
  path "${id}.bam.bai"

  script:
  """
  prefetch $id
  fasterq-dump $id
  bwa-mem2 mem -t 4 ${params.refgenome} ${id}.fastq | \
    samtools sort -@ 4 -o ${id}.bam 
  samtools index ${id}.bam
  rm ${id}.fastq
  rm -r $id
  """

}


process fetchAndMapToRefGenomePaired {

  // if very big run can be better not stop but handle failed cases afterwards
  //errorStrategy 'ignore'

  publishDir params.bamdir, mode: 'move', overwrite: true
 
  input:
  val id

  output:
  path "${id}.bam" 
  path"${id}.bam.bai"

  script:
  """
  prefetch $id
  fasterq-dump --split-files $id
  bwa-mem2 mem -t ${params.numthreads}  ${params.refgenome} \
    ${id}_1.fastq \
    ${id}_2.fastq | \
    samtools sort -@ ${params.numthreads} -o ${id}.bam 
  samtools index ${id}.bam
  rm ${id}_*.fastq
  rm -r $id
  """

}


/*
 *  For locally stored single fastq files 
 */
process mapFastqsToRefGenomeSingle {

  publishDir params.bamdir, mode: 'move', overwrite: true

  input:
  tuple val(name), file(fastqfile)
  
  output:
  file "${name}.bam"

  script:
  """
  echo mapping ${name}
  bwa-mem2 mem -t ${params.numthreads} ${params.refgenome} $fastqfile | \
    samtools sort -@ ${params.numthreads} -o ${name}.bam
  samtools index ${name}.bam
  """

}


/*
 * For locally stored paired fastq files
 */
process mapFastqsToRefGenomePaired {

  publishDir params.bamdir, mode: 'move', overwrite: true
 
  input:
  tuple val(name), file(fastqpair)

  output:
  file "${name}.bam"
  
  script:
  """
  bwa-mem2 mem -t ${params.numthreads} ${params.refgenome} ${fastqpair} | \
    samtools sort -@ ${params.numthreads} -o ${name}.bam
  samtools index ${name}.bam 
  """

}


/*
  run stats and collect some of the key figures
  - use sampling at least if bams > 2GB
  TODO: add automated plotting, needs gnuplot
*/
process makeBamStatsAndCigars {

  publishDir statdir, mode: 'copy', overwrite: true

  input:
  tuple val(name), path(bamfile)

  output:
  path "${name}.bc"
  path "${name}.popular-cigars.txt"
  //path "${name}-coverage.png"

  script:
  if (params.samplebigbams == false) {
    """
    samtools stats -@ 3 $bamfile >${name}.bc
    #plot-bamstats -p ${name} ${name}.bc
    samtools flagstat $bamfile | grep 'primary mapped' > ${name}.popular-cigars.txt
    samtools view $bamfile | cut -f 6 | sort | uniq -c | sort -k1,2nr  | head >> ${name}.popular-cigars.txt
    """
  } else if (params.samplebigbams == true) {
    """
    samtools stats -@ 3 $bamfile >${name}.bc
    #plot-bamstats -p ${name} ${name}.bc
    samtools flagstat $bamfile | grep 'primary mapped' > ${name}.popular-cigars.txt
    samtools view --subsample 0.1 $bamfile | cut -f 6 | sort | uniq -c | sort -k1,2nr  | head >> ${name}.popular-cigars.txt
    cat tmp1.txt tmp2.txt > ${name}.popular-cigars.txt
    """
  } else {
    error "Could not figure bam sampling"
  }
}



/*
  Collect coverage stats from all samples to one file 
*/
process collectCoverageStats {

  publishDir statdir, mode: 'copy', overwrite: true  

  input:
  path statfilelist
  path cigarfilelist

  output:
  path "COV_all.tsv"
  path "CIGARS.txt"

  script:
  """
  grep '^COV' $statfilelist >COV_all.tsv
  head -n 4 $cigarfilelist >CIGARS.txt 
  """
}
