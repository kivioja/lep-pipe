# lep-pipe

Nextflow workflows for creating linkage maps and associated tasks using [Lep-MAP3](https://sourceforge.net/projects/lep-map3/) and [Lep-Anchor](https://sourceforge.net/projects/lep-anchor/). 

The workflows use as input sequence data (for example RAD-seq or WGS) from one or several families. 

Contains the following workflows

* main: workflow for generating linkage maps, assumes that read mapping has been done and bam the files exist 
* seq_fw: simple work workflow for mapping sequences to a reference genome
* rel_wf: various quality control processes to check for example the relatedness of individuals
* post_wf: some posprocessing of results




