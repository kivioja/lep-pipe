#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
  get fastq and map to reference genome, for large WGS sets
  - does not store fastqs to save space
  - uses bwa-mem2 to save time
*/

/*
  UNDER CONSTRUCTION
  ena-file-downloader sometimes fails for unknown reason
  - replace with enaDataGet?
*/

// defaults
params.fastqtype = "single"

process fetchAndMapToRefGenomeSingle {

  publishDir params.bamdir, mode: 'move', overwrite: true
 
  input:
  val id

  output:
  path "${id}.bam"
  path "${id}.bam.bai"

  script:
  """
  java -jar ${params.enaloaderpath}/ena-file-downloader.jar \
        --accessions=${id} \
        --format=READS_FASTQ \
        --location=\$PWD \
        --protocol=FTP --asperaLocation=null --email=NONE
  bwa-mem2 mem -t 4 ${params.refgenome} \
    reads_fastq/${id}/${id}.fastq.gz | \
    samtools sort -@ 4 -o ${id}.bam 
  samtools index ${id}.bam
  rm reads_fastq/${id}/${id}.fastq.gz
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
  java -jar ${params.enaloaderpath}/ena-file-downloader.jar \
        --accessions=${id} \
        --format=READS_FASTQ \
        --location=\$PWD \
        --protocol=FTP --asperaLocation=null --email=NONE
  bwa-mem2 mem -t 4 ${params.refgenome} \
    reads_fastq/${id}/${id}_1.fastq.gz \
    reads_fastq/${id}/${id}_2.fastq.gz | \
    samtools sort -@ 4 -o ${id}.bam 
  samtools index ${id}.bam
  rm reads_fastq/${id}/${id}_1.fastq.gz
  rm reads_fastq/${id}/${id}_2.fastq.gz
  """

}


