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
    file "crossovers${chrom}.txt"
    

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
        awk -f ${params.scriptdir}/numRec3.awk fit${chrom}.txt >crossovers${chrom}.txt
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
        awk -f ${params.scriptdir}/numRec3.awk fit${chrom}.txt >crossovers${chrom}.txt
        """
    } 
}



process fitStepReferenceAndEstimateRecomb {

    publishDir params.recdensitydir, mode: 'copy', overwrite: true

    input:
    tuple val(chrom), val(refchrom), path(mapforla), path(mapforstep)

    output:
    file "fit${refchrom}.txt"
    file "estimates${refchrom}.txt"
    file "density${refchrom}.txt"
    file "crossovers${refchrom}.txt"
    file "command_fit_${refchrom}.err"

    script:
    if (params.bandwidth == 0) {
        """
        java ${params.jvmoptions} -cp ${params.lepanchordir}/bin FitStepFunction \
            intervals=${mapforla} \
            map=${mapforstep} \
            chromosome=${chrom} | perl -ple 's/^([\\w\\.]+)/CHR$refchrom/' >fit${refchrom}.txt
        java ${params.jvmoptions} -cp ${params.lepanchordir}/bin EstimateRecombination \
            gridSize=${params.gridsize} \
            numSamples=${params.numsamples} \
            minSupport=${params.minsupport} \
            endSupport=${params.endsupport} \
            minLength=${params.minlength} \
            crossValidate=1 \
            map=fit${refchrom}.txt >estimates${refchrom}.txt
        grep 'density' estimates${refchrom}.txt >density${refchrom}.txt
        awk -f ${params.scriptdir}/numRec3.awk fit${refchrom}.txt >crossovers${refchrom}.txt
        cp .command.err command_fit_${refchrom}.err
        """
    } else {
        """
        java ${params.jvmoptions} -cp ${params.lepanchordir}/bin FitStepFunction \
            intervals=${mapforla} \
            map=${mapforstep} \
            chromosome=${chrom} | perl -ple 's/^([\\w\\.]+)/CHR$refchrom/' >fit${refchrom}.txt
        java ${params.jvmoptions} -cp ${params.lepanchordir}/bin EstimateRecombination \
            gridSize=${params.gridsize} \
            numSamples=${params.numsamples} \
            minSupport=${params.minsupport} \
            endSupport=${params.endsupport} \
            minLength=${params.minlength} \
            bandwidth=${params.bandwidth} \
            map=fit${refchrom}.txt >estimates${refchrom}.txt
        grep 'density' estimates${chrom}.txt >density${refchrom}.txt
        awk -f ${params.scriptdir}/numRec3.awk fit${refchrom}.txt >crossovers${refchrom}.txt
        cp .command.err command_fit_${refchrom}.err
        """
    } 
}




/*
 * evaluate markers in the physical order, do not change order, 
 * and estimate recombination
 */
process evaluateOrderAndEstimateRecomb {
  publishDir params.recdensitydir, mode: 'copy', overwrite: true 
 
  input:
  path 'snps.txt'
  path callsfile
  path finalmap 
  each chrom 

  output:
  path "order${chrom}_phys_eval_coords.txt"
  path "estimates${chrom}.txt"
  path "density${chrom}.txt"
  path "crossovers${chrom}.txt"

  // The input file that has markers in physical order is done like shown in Lep-Anchor documentation 
  // except that the marker number is stored to .input file in the first step so that 
  // it carries to .phys file that is then be used after the evaluation to 
  // get the linkage group coordinates back  
  script:
  """
  awk -vn=${chrom} '(NR==FNR){map[NR-1]=\$0}(NR!=FNR){\$1=map[\$1] "\t" \$1 "\t";print}' snps.txt ${params.mapdir}/orders/order${chrom}.int >order${chrom}.input
  awk -f ${params.lepanchordir}/liftover ${params.mapdir}/anchoring/fixed_chr${chrom}.agp order${chrom}.input | sort -V | grep LG >order${chrom}.liftover
  awk -vinverse=1 -f ${params.lepanchordir}/liftover ${params.mapdir}/anchoring/fixed_chr${chrom}.agp order${chrom}.liftover | \
    awk '(NR==FNR){m[\$1"\t"(\$2+0)]=NR-1}(NR!=FNR){print m[\$1"\t"(\$2+0)]}' snps.txt - >order${chrom}.phys
  zcat $callsfile | java ${params.jvmoptions} -Xmx${params.javaheapsize} -cp ${params.lepmapdir} OrderMarkers2 data=- \
    useMorgan=1 \
    scale=${params.scale} \
    proximityScale=${params.proximityscale} \
    evaluateOrder=order${chrom}.phys improveOrder=0 outputPhasedData=2 phasingIterations=3 \
    map=${finalmap} chromosome=$chrom calculateIntervals=order_${chrom}_phys_eval.int \
    >order${chrom}_phys_eval.txt
  marker_nums_to_coords.pl order${chrom}.liftover < order${chrom}_phys_eval.txt >order${chrom}_phys_eval_coords.txt
  java ${params.jvmoptions} -cp ${params.lepanchordir}/bin EstimateRecombination \
    gridSize=${params.gridsize} \
    numSamples=${params.numsamples} \
    minSupport=${params.minsupport} \
    endSupport=${params.endsupport} \
    minLength=${params.minlength} \
    crossValidate=1 \
    map=order${chrom}_phys_eval_coords.txt >estimates${chrom}.txt
  grep 'density' estimates${chrom}.txt >density${chrom}.txt
  awk -f ${params.scriptdir}/numRec3.awk order${chrom}_phys_eval_coords.txt >crossovers${chrom}.txt
  """
}


/*
 * evaluate markers in the reference genome order, do not change order, 
 * and estimate recombination
 * sorts agp according to LG and uses that as the chromosome ordering
 * and after that moves to reference chromosome numbering 
 */
process evaluateReferenceOrderAndEstimateRecomb {
  publishDir params.recdensitydir, mode: 'copy', overwrite: true 
 
  input:
  tuple val(chrom), val(refchrom), path(snpfile), path(callsfile), path(finalmap)

  output:
  path "order${refchrom}_phys_eval_coords.txt"
  path "estimates${refchrom}.txt"
  path "crossovers${refchrom}.txt"
  path "density${refchrom}.txt"

  // The input file that has markers in physical order is done like shown in Lep-Anchor documentation 
  // except that the marker number is stored to .input file in the first step so that 
  // it carries to .phys file that is then be used after the evaluation to 
  // get the linkage group coordinates back  
  script:
  if (params.nofemalerecomb == false) {
    """
    sort -V $params.refagp | head -n $chrom | tail -n 1 >chrom.agp
    awk -vn=${chrom} '(NR==FNR){map[NR-1]=\$0}(NR!=FNR){\$1=map[\$1] "\t" \$1 "\t";print}' $snpfile ${params.mapdir}/orders/order${chrom}.int >order${chrom}.input
    awk -f ${params.lepanchordir}/liftover chrom.agp order${chrom}.input | sort -V | grep LG >order${chrom}.liftover
    awk -vinverse=1 -f ${params.lepanchordir}/liftover chrom.agp order${chrom}.liftover | \
      awk '(NR==FNR){m[\$1"\t"(\$2+0)]=NR-1}(NR!=FNR){print m[\$1"\t"(\$2+0)]}' $snpfile - >order${chrom}.phys
    zcat $callsfile | java ${params.jvmoptions} -Xmx${params.javaheapsize} -cp ${params.lepmapdir} OrderMarkers2 data=- \
      useMorgan=1 \
      scale=${params.scale} \
      proximityScale=${params.proximityscale} \
      evaluateOrder=order${chrom}.phys improveOrder=0 outputPhasedData=2 phasingIterations=3 \
      map=${finalmap} chromosome=$chrom calculateIntervals=order_${chrom}_phys_eval.int \
      >order${chrom}_phys_eval.txt
    marker_nums_to_coords.pl order${chrom}.liftover < order${chrom}_phys_eval.txt | \
      perl -ple 's/^LG$chrom/CHR$refchrom/'  > order${refchrom}_phys_eval_coords.txt
    java ${params.jvmoptions} -cp ${params.lepanchordir}/bin EstimateRecombination \
      gridSize=${params.gridsize} \
      numSamples=${params.numsamples} \
      minSupport=${params.minsupport} \
      endSupport=${params.endsupport} \
      minLength=${params.minlength} \
      crossValidate=1 \
      map=order${refchrom}_phys_eval_coords.txt >estimates${refchrom}.txt
    grep 'density' estimates${refchrom}.txt >density${refchrom}.txt
    awk -f ${params.scriptdir}/numRec3.awk order${refchrom}_phys_eval_coords.txt >crossovers${refchrom}.txt
    """
  } else {
    """
    sort -V $params.refagp | head -n $chrom | tail -n 1 >chrom.agp
    awk -vn=${chrom} '(NR==FNR){map[NR-1]=\$0}(NR!=FNR){\$1=map[\$1] "\t" \$1 "\t";print}' $snpfile ${params.mapdir}/orders/order${chrom}.int >order${chrom}.input
    awk -f ${params.lepanchordir}/liftover chrom.agp order${chrom}.input | sort -V | grep LG >order${chrom}.liftover
    awk -vinverse=1 -f ${params.lepanchordir}/liftover chrom.agp order${chrom}.liftover | \
      awk '(NR==FNR){m[\$1"\t"(\$2+0)]=NR-1}(NR!=FNR){print m[\$1"\t"(\$2+0)]}' $snpfile - >order${chrom}.phys
    zcat $callsfile | java ${params.jvmoptions} -Xmx${params.javaheapsize} -cp ${params.lepmapdir} OrderMarkers2 data=- \
      useMorgan=1 \
      scale=${params.scale} \
      proximityScale=${params.proximityscale} \
      recombination2=0 \
      evaluateOrder=order${chrom}.phys improveOrder=0 outputPhasedData=2 phasingIterations=3 \
      map=${finalmap} chromosome=$chrom calculateIntervals=order_${chrom}_phys_eval.int \
      >order${chrom}_phys_eval.txt
    marker_nums_to_coords.pl order${chrom}.liftover < order${chrom}_phys_eval.txt | \
      perl -ple 's/^LG$chrom/CHR$refchrom/'  > order${refchrom}_phys_eval_coords.txt
    java ${params.jvmoptions} -cp ${params.lepanchordir}/bin EstimateRecombination \
      gridSize=${params.gridsize} \
      numSamples=${params.numsamples} \
      minSupport=${params.minsupport} \
      endSupport=${params.endsupport} \
      minLength=${params.minlength} \
      crossValidate=1 \
      map=order${refchrom}_phys_eval_coords.txt >estimates${refchrom}.txt
    grep 'density' estimates${refchrom}.txt >density${refchrom}.txt
    awk -f ${params.scriptdir}/numRec3.awk order${refchrom}_phys_eval_coords.txt >crossovers${refchrom}.txt
    """
  }
}







process estimateGammaModel {

    publishDir params.gammamodeldir, mode: 'copy', overwrite: true

    input:
    each chrom
    path "*"

    output:
    file "estimates${chrom}.txt"
    
    script:
    """
    java ${params.jvmoptions} -cp ${params.lepanchordir}/bin EstimateRecombination \
        gammaModel=${params.gammamodel} \
        gridSize=${params.gridsize} \
        numSamples=${params.numsamples} \
        minSupport=${params.minsupport} \
        endSupport=${params.endsupport} \
        minLength=${params.minlength} \
        crossValidate=1 \
        map=fit${chrom}.txt >estimates${chrom}.txt
    """

}
