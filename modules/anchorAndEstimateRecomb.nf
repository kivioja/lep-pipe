#!/usr/bin/env nextflow

/*
  run Lep-Anchor  
*/


/*
 * Lep-Anchor through the wrapper
 * - wrapper assumes that the other awk scripts are in the same directory
 *   and classes in ./bin
 */
process runAnchorWrapper {   
    
    publishDir params.anchordir, mode: 'copy', overwrite: true


    // Note that REF_LA.agp has both tabs and spaces as separator
    def final_agp_name = params.species + "_LA.agp"

    input: 
    path 'bin'
    path 'map_for_la.txt'

    output:
    path "marey*.png"
    path "chr*.agp"
    path  "$final_agp_name"
    path 'map.bed'
    path 'marey.data.gz'

    // if spaces to tabs, escaping within the perl regular expression does give the command
    // perl -ple 's/\s/\t/g' REF_LA.agp >Spotted_gar_LA.agp
    script:
    """
    cp ${params.lepanchordir}/lepanchor_wrapper2.sh .
    cp ${params.lepanchordir}/*.awk .
    cp ${params.lepanchordir}/plot_marey.R .
    cat ${params.refgenome} | awk -f ${params.lepanchordir}/contigLength.awk >contigs.length
    bash lepanchor_wrapper2.sh -n ${params.chromnum} -m map_for_la.txt -t ${params.numthreads}
    perl -ple 's/\\s/\\t/g' REF_LA.agp >${final_agp_name}
    """

}


/*
 * let Lep-Anchor to fix issues within contigs 
 */
process fixContigs {

    publishDir params.anchordir, mode: 'copy', overwrite: true

    def chroms = 1..params.chromnum

    input:
    path 'map_for_la.txt'
    path 'map.bed'
    each chrom 

    output:
    path "fixed_chr${chrom}.agp"
    path "fixed_marey_chr${chrom}.png"
    path "marey_chr${chrom}_fixed.data.gz"

    // note: need to separate groovy and awk variables
    shell:
    """
    java !{params.jvmoptions} -Xmx!{params.javaheapsize} -cp !{params.lepanchordir}/bin PlaceAndOrientContigs \
        map=map_for_la.txt bed=map.bed chromosome=!{chrom} \
        findContigErrors=1 minImprovement=!{params.LepAnchor_minImprovement} \
        >map_chr_fix.la
    awk -vlg=!{chrom} -vprefix=LG -f !{params.lepanchordir}/makeagp_full2.awk \
        map_chr_fix.la \
        >fixed_chr!{chrom}.agp 
    awk '(\$3==!{chrom})' map_for_la.txt | awk -f !{params.lepanchordir}/liftover.awk fixed_chr!{chrom}.agp - | \
        awk '(/LG/ && NF>=4) { \
            for (i=4;i<NF;i+=2) print \$1"\t"\$2"\t"\$3"\t"1"\t"\$(i)"\t"\$(i+1)}' | \
        gzip >marey_chr${chrom}_fixed.data.gz
    Rscript !{params.scriptdir}/plotMareyWithArgs.R \
        --data marey_chr${chrom}_fixed.data.gz \
        --agp fixed_chr \
        --out fixed_marey_chr
    """
}



/*
 * Lift-over the maps to the linkage group coordinates 
 * - use the fixed contigs, need to collect to one agp file first
 * - use the agp file to generate a new genome fasta file too
 */
process liftOverMaps {

    publishDir params.recdensitydir, mode: 'copy', overwrite: true

    def fixed_agp_name = params.species + "_LAfixed.agp"
    def fixed_genome_name = params.species + "_LAfixed.fa"

    input:
    file 'fixed_chr*.agp'
    file 'map_for_la.txt'
    file 'map_for_step.txt'

    output:
    file "$fixed_agp_name"
    file "map_for_la.liftover.txt"
    file "map_for_step.liftover.txt"
    file "$fixed_genome_name"

    // lift-over of input files to linkage group coordinates for fitting
    // awk [-vinverse=1] -f liftover.awk ref.agp chr_pos_file >chr_pos_file.liftover
    script:
    """    
    cat fixed_chr*.agp | sort -k1,1 -k2,2n > $fixed_agp_name
    awk -f ${params.lepanchordir}/liftover.awk $fixed_agp_name map_for_la.txt >map_for_la.liftover.txt
    awk -f ${params.lepanchordir}/liftover.awk $fixed_agp_name map_for_step.txt >map_for_step.liftover.txt
    seqkit seq -i $params.refgenome > tmpgenome.fa
    agptools assemble tmpgenome.fa $fixed_agp_name > $fixed_genome_name
    rm tmpgenome.fa
    """
}



process fitStepAndEstimateRecomb {

    publishDir params.recdensitydir, mode: 'copy', overwrite: true

    input:
    file 'map_for_la.liftover.txt'
    file 'map_for_step.liftover.txt'
    each chrom

    output:
    file "fit${chrom}.txt"
    file "estimates${chrom}.txt"
    file "density${chrom}.txt"
    

    script:
    if (params.bandwidth == 0) {
        """
        java ${params.jvmoptions} -cp ${params.lepanchordir}/bin FitStepFunction \
            intervals=map_for_la.liftover.txt \
            map=map_for_step.liftover.txt \
            chromosome=${chrom} >fit${chrom}.txt
        java ${params.jvmoptions} -cp ${params.lepanchordir}/bin EstimateRecombination \
            gridSize=${params.gridsize} \
            numSamples=${params.numsamples} \
            minSupport=${params.minsupport} \
            endSupport=${params.endsupport} \
            minLength=${params.minlength} \
            crossValidate=1 \
            map=fit${chrom}.txt >estimates${chrom}.txt
        grep 'density' estimates${chrom}.txt >density${chrom}.txt
        """
    } else {
        """
        java ${params.jvmoptions} -cp ${params.lepanchordir}/bin FitStepFunction \
            intervals=map_for_la.liftover.txt \
            map=map_for_step.liftover.txt \
            chromosome=${chrom} >fit${chrom}.txt
        java ${params.jvmoptions} -cp ${params.lepanchordir}/bin EstimateRecombination \
            gridSize=${params.gridsize} \
            numSamples=${params.numsamples} \
            minSupport=${params.minsupport} \
            endSupport=${params.endsupport} \
            minLength=${params.minlength} \
            bandwidth=${params.bandwidth} \
            map=fit${chrom}.txt >estimates${chrom}.txt
        grep 'density' estimates${chrom}.txt >density${chrom}.txt
        """
    } 
}
