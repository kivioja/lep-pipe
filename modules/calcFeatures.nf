#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
  do the variant calling and use variants to check the parents 
  TODO: reduce disk usage
*/

params.numthreads = 4

/* 
  parse parent calls and lift over to new genome coordinates  
*/

process calcNucDiversity {

  input:
  each parentcall
  path bam2ind
  path agpfile

  output:
  path 'result.bed'

  script:
  """
  zcat $parentcall | head -n 7 > call_header.ped 
  zcat $parentcall | tail -n +8 | cut -f 1,2 > post.pos.txt
  awk -f ${params.lepanchordir}/liftover.awk $agpfile post.pos.txt > post.liftover.pos.txt
  zcat $parentcall | perl ${params.scriptdir}/parse_parent_calls.pl \
    call_header.ped post.pos.txt post.liftover.pos.txt >result.bed
  """

}



/*
  collect the files that contain heterozygotes and count in windows
  - merge and sort the bed files
  - count for individual in each window
    the number heterozygous and homozygous sites 
*/
process collectNucDiversities {

  publishDir params.outdir, mode: 'copy', overwrite: true

  input:
  path 'parent_calls*.bed'
  path 'intervals.bed'

  output:
  path 'parent_calls_interval_counts.txt'

  script:
  """
  sort -k1,1 -k2,2n intervals.bed >intervals_sorted.bed
  sort -k1,1 -k2,2n -T ${params.bigtempdir} parent_calls*.bed | \
    bedtools intersect -sorted -wo -a intervals_sorted.bed -b - | cut -f 1-3,7-8 | \
    sort -k1,1 -k2,2n -k4,4 -k5,5 -T ${params.bigtempdir}| uniq -c >parent_calls_interval_counts.txt  
  """
}








