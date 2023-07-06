#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * Workflow for sequence preprocessing before linkage mapping
 * EXPERIMENTAL: the download from ENA fails occasionally  
 */

/*
  include processes 
*/
include { fetchAndMapToRefGenomeSingle; fetchAndMapToRefGenomePaired } from './modules/fetchAndMap.nf'

/*
 Fetch data from ENA and map to reference genome
 EXPERIMENTAL, download using ENA tools fails occasionally 
 */
workflow {

    // first column should contain the run names
    run_ch = Channel.fromPath(params.run2indfile)
                    .splitCsv(skip: 1)
                    .map{ row -> row[0] }

    // map to reference, single or paired-end reads 
    if (params.fastqtype == 'single') {
        fetchAndMapToRefGenomeSingle(run_ch)
    } else if (params.fastqtype == 'paired') {
        fetchAndMapToRefGenomePaired(run_ch)
    } else {
        error("Unknown fastq type")
    }
}
