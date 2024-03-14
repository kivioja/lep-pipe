#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * Workflow for sequence preprocessing before linkage mapping  
 */

/*
  include processes 
*/
include { fetchAndMapToRefGenomeSingle; fetchAndMapToRefGenomePaired; mapFastqsToRefGenomePaired; mapFastqsToRefGenomeSingle} from './modules/fetchAndMap.nf'


// by default assumes that the files are fetched from ENA
params.fetchfiles = true
params.fastqtype = "single"


workflow {

    // map fastqs to reference genome
    if (params.fetchfiles == false) {
        if (params.fastqtype == 'paired') {
            fastq_ch = Channel.fromFilePairs(params.fastqtemplate)
            mapFastqsToRefGenomePaired(fastq_ch)
        } else {
            // note, get filename without suffix
            fastq_ch = Channel.fromPath(params.fastqtemplate)
                .map { file -> tuple(file.simpleName, file) }
            mapFastqsToRefGenomeSingle(fastq_ch)
        }   
    } else {
        // fetch from ENA and map
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

}
