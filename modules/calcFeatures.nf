#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * do the variant calling and use variants to check the parents 
 */

params.outdir = params.mapdir + "/" + "features"

/* 
 * Make the windows in which the features are counted and count basic features
 */
process makeGenomeWindows {
  publishDir params.outdir, mode: 'copy', overwrite: true

  input:
  val infiletemplate
  val bedoutfile
  val recoutfile
  val windowsize

  output:
  path "$bedoutfile"
  path "$recoutfile"

  script:
  """
  Rscript ${params.scriptdir}/collectMareyWindowsBed.R \
    -n $params.chromnum \
    -w $windowsize \
    -i $infiletemplate \
    -b $bedoutfile \
    -r $recoutfile   
  """
}



/*
 * Count the basic sequence features in the genomic windows 
 */
process countWindowSeqFeatures {
  publishDir params.outdir, mode: 'copy', overwrite: true
  
  input:
  path windowbed
  path genomefasta
  val seqfeatcountfile

  output:
  path "$seqfeatcountfile"
 
  script:
  """
  seqkit subseq --bed $windowbed $genomefasta > tmpseqfile.fa
  seqtk comp tmpseqfile.fa >testout.txt
  seqtk comp tmpseqfile.fa | perl -ple 's/^(LG\\d+)_(\\d+)-(\\d+):\\./\$1\\t\$2\\t\$3/' >$seqfeatcountfile
  """  
}



/* 
 * Call and parse to bed all variants, not only the ones useful for linkage mapping
 * and lift to chromosome coordinates
 * - long pipe to avoid huge temporary files
 * TODO: add more parent calling options?
 */ 
process callToBedAllParentVars {
 
  input:
  each path(contigchunk)
  path bam2ind
  path ped
  path agpfile

  output:
  path 'result.bed.gz'

  script:
  """
  cut -f 1 -d ',' $bam2ind | \
  perl -ple 's{^}{${params.bamdir}/}' >bamlist.txt
  cut -f 2 -d ',' $bam2ind  > mapping.txt
  samtools mpileup -l $contigchunk -q 10 -Q 10 -s --bam-list bamlist.txt | \
    java ${params.jvmoptions} -cp ${params.lepmapdir} Pileup2Likelihoods \
      numLowerCoverage=${params.LepMap_numLowerCoverage} \
      minAlleleFreq=0.0 \
      mappingFile=mapping.txt | \
    java ${params.jvmoptions} -cp ${params.lepmapdir} ParentCall2 \
	    data=${ped} \
	    posteriorFile=- \
      halfSibs=${params.halfsibs} \
      removeNonInformative=0 | \
  parse_parent_calls.pl $ped | gzip >result.contig.txt.gz
  awk -f ${params.lepanchordir}/liftover.awk $agpfile <(zcat result.contig.txt.gz) | \
    awk -v OFS='\t' '{print \$1, \$2-1, \$2, \$3, \$4}' | gzip > result.bed.gz  
  rm result.contig.txt.gz
  """ 
}



/*
  collect the files that contain heterozygotes and count in windows
  - merge and sort the bed files
  - count for individual in each window
    the number heterozygous and homozygous sites 
*/
process collectParentVariantWindowCounts {

  publishDir params.outdir, mode: 'copy', overwrite: true

  input:
  path 'parent_calls*.bed.gz'
  path 'intervals.bed'
  val countoutfile

  output:
  path "$countoutfile"

  script:
  """
  sort -k1,1 -k2,2n intervals.bed >intervals_sorted.bed
  zcat parent_calls*.bed.gz | sort -k1,1 -k2,2n -T ${params.bigtempdir} | \
    bedtools intersect -sorted -wo -a intervals_sorted.bed -b - | cut -f 1-3,7-8 | \
    sort -k1,1 -k2,2n -k4,4 -k5,5 -T ${params.bigtempdir}| uniq -c > $countoutfile  
  """
}








