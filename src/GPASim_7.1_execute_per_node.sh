#!/bin/bash

### INPUT PARAMETERS
DIR=${1}                # ( 1) output directory
SRC_DIR=${2}            # ( 2) directory where quantine2 folder, plink, emmax, gemma, gcta, and genomic_prediction folder are located
rep=${3}                # ( 3) replicate number
nIndividuals=${4}       # ( 4) number of individuals to simulate
nLoci=${5}              # ( 5) number of loci to simulate (neutral + QTL)
nQTL=${6}               # ( 6) number of QTL among the all the loci simulated (should range from 1 to nLoci-1)
nBGS=${7}               # ( 7) number of background selection loci (should range from 1 to nLoci-1; conditional on nQTL where nQTL + nBGS <= nLoci)
nAlleles=${8}           # ( 8) number of alleles per loci, e.g. 5 for A,T,C,G, and DEL (NOTE: We're now assuming multi-allelic SNP loci are just noise so we're simulating only biallelic loci!)
allele_eff_model=${9}   # ( 9) probability distribution to sample the allele effects from (i.e. CHISQ, NORM, & UNIF)
nGen=${10}              # (10) number of generations to simulate (approximating the number of generations since the introduction of Lolium in Australia in the 1880s minus the time when extensive herbicide application started)
nPop=${11}              # (11) number of populations/subpopulations to simulate (NOTE: must have a natural number square-root)
migration=${12}         # (12) migration rate across the populations (surrently using the 1D stepping stone model see line 117)
selection=${13}         # (13) proportion of individual selected based on trait of interest (nQTL) defined as proportion of the most fit individuals selected, where fitness is defined using the generalized logistic curve (Richards, 1959): ranges from 0.00 to 1.00; lower s means more intense selection
bg_selection=${14}      # (14) proportion of individual selected based on background selection trait (nBGS) defined as proportion of the most fit individuals selected, where fitness is defined using the generalized logistic curve (Richards, 1959): ranges from 0.00 to 1.00; lower s means more intense selection
GRADIENT=${15}          # (15) uniformly distributed non-wildtype alleles
NTHREADS=${16}          # (16) number of computing cores or threads to run in parallel for some of the modules
NPOOLS=${17}            # (17) number of pools to genotype per population (NOTE: for across population sample 1 population:1 pool)
NLIB=${18}              # (18) maximum number of individually barcoded libraries for sequencing (Indi-seq and Pool-seq) (NOTE: May need to have a minimum limit at 1,000 to keep GCTA from failing due to small population size leading to REML convergence failure)
NSEQ_MAX=${19}          # (19) maximum number of multiplexed sequencing libraries (limiting factor in sampling strategy optimization)

# ###############
# ### TESTING ###
# DIR=/data/Lolium/Quantitative_Genetics/LOLSIM2020
# SRC_DIR=/data/Lolium/Softwares
# rep=1
# nIndividuals=1000
# nLoci=1000
# nQTL=10
# nBGS=10
# nAlleles=2
# allele_eff_model=CHISQ
# nGen=200
# nPop=16
# migration=0.001
# selection=0.05
# bg_selection=0.10
# GRADIENT=1
# NTHREADS=12; SPARTAN_CORES=12
# NPOOLS=5
# NLIB=384
# NSEQ_MAX=10
# ###############

echo "##################################"
echo "### PRE-COMPILE JULIA PACKAGES ###"
echo "##################################"
echo 'using GWAlpha
      using CSV
      using Distributions
      using DataFrames
      using DelimitedFiles
      using LinearAlgebra
      using Optim
      using UnicodePlots
      using ColorBrewer
      using ProgressMeter
      using Statistics
      using StatsBase
      using Distributed
      using SharedArrays
      using GeneticVariation
      using RCall
' > precompile-${rep}REP_${nQTL}NQTL_${nBGS}NBGS_${migration}MIG_${selection}SEL_${bg_selection}BGSEL_${GRADIENT}GRAD.jl
### run thrice just to be sure
julia precompile-${rep}REP_${nQTL}NQTL_${nBGS}NBGS_${migration}MIG_${selection}SEL_${bg_selection}BGSEL_${GRADIENT}GRAD.jl
julia precompile-${rep}REP_${nQTL}NQTL_${nBGS}NBGS_${migration}MIG_${selection}SEL_${bg_selection}BGSEL_${GRADIENT}GRAD.jl
julia precompile-${rep}REP_${nQTL}NQTL_${nBGS}NBGS_${migration}MIG_${selection}SEL_${bg_selection}BGSEL_${GRADIENT}GRAD.jl
rm precompile-${rep}REP_${nQTL}NQTL_${nBGS}NBGS_${migration}MIG_${selection}SEL_${bg_selection}BGSEL_${GRADIENT}GRAD.jl

echo "################################"
echo "### PREPARE OUTPUT DIRECTORY ###"
echo "################################"
### prefix for the quantiNemo2 initiation (*.ini) file and the output folder
OUTPREFIX=GPASim_${rep}REP_${nQTL}NQTL_${nBGS}NBGS_${migration}MIG_${selection}SEL_${bg_selection}BGSEL_${GRADIENT}GRAD
### output directory to dump all the output files and folder
OUTDIR=${DIR}/${OUTPREFIX}
### create the output directory
mkdir ${OUTDIR}
### navigate into the working directory (not critical we just want to save the nohup files there)
cd ${OUTDIR}

echo "################"
echo "### SIMULATE ###"
echo "################"
time \
nohup \
${SRC_DIR}/genomic_prediction/src/GPASim_1.0_simulate.sh \
      -s $SRC_DIR \
      -d $OUTDIR \
      -p $OUTPREFIX \
      -n $nIndividuals \
      -l $nLoci \
      -q $nQTL \
      -b $nBGS \
      -a $nAlleles \
      -r $allele_eff_model \
      -t $nGen \
      -P $nPop \
      -m $migration \
      -sq $selection \
      -sb $bg_selection \
      -g $GRADIENT
mv nohup.out GPASim_1.0_simulate.nohup

echo "#############"
echo "### PARSE ###"
echo "#############"
time \
nohup \
${SRC_DIR}/genomic_prediction/src/GPASim_2.0_parse.sh \
      -d $OUTDIR \
      -s $SRC_DIR \
      -t $NTHREADS \
      -p $NPOOLS
mv nohup.out GPASim_2.0_parse.nohup

echo "################"
echo "### SAMPLING ###"
echo "################"
time \
nohup \
${SRC_DIR}/genomic_prediction/src/GPASim_3.0_across_population_sampling.sh \
      -d $OUTDIR \
      -s $SRC_DIR \
      -t $NTHREADS \
      -l $NLIB
mv nohup.out GPASim_3.0_across_population_sampling.nohup

echo "#################################################"
echo "### FILTER GENOTYPE DATA AND BUILD COVARIATES ###"
echo "#################################################"
time \
nohup \
${SRC_DIR}/genomic_prediction/src/GPASim_4.0_filter_build_covariates.sh \
      -d $OUTDIR \
      -s $SRC_DIR \
      -t $NTHREADS
mv nohup.out GPASim_4.0_filter_build_covariates.nohup

echo "############"
echo "### GPAS ###"
echo "############"
time \
nohup \
${SRC_DIR}/genomic_prediction/src/GPASim_5.0_GPAS.sh \
      -d $OUTDIR \
      -s $SRC_DIR \
      -t $NTHREADS
mv nohup.out GPASim_5.0_GPAS.nohup

echo "################################################"
echo "### MERGE AND SAMPLING STRATEGY OPTIMIZATION ###"
echo "################################################"
time \
nohup \
${SRC_DIR}/genomic_prediction/src/GPASim_6.0_merging_and_abc_optim.sh \
    --dir ${OUTDIR} \
    --src-dir ${SRC_DIR} \
    --replicate-id ${rep} \
    --n-individuals ${nIndividuals} \
    --n-loci ${nLoci} \
    --n-QTL ${nQTL} \
    --n-BGS ${nBGS} \
    --n-alleles ${nAlleles} \
    --alleles-dist ${allele_eff_model} \
    --n-generations ${nGen} \
    --n-populations ${nPop} \
    --migration-rate ${migration} \
    --selection-intensity-QTL ${selection} \
    --selection-intensity-BGS ${bg_selection} \
    --diffusion-gradient ${GRADIENT} \
    --n-pools ${NPOOLS} \
    --n-libraries ${NLIB} \
    --max-sequencing ${NSEQ_MAX}
mv nohup.out GPASim_6.0_merging_and_abc_optim.nohup

echo "#################"
echo "### DEBUGGING ###"
echo "#################"
grep -i -f <(echo -e "err\nfail") GPASim_*.nohup > ERR.nohup
wc -l ERR.nohup

echo "################"
echo "### CLEAN-UP ###"
echo "################"
cd ${OUTDIR}
mkdir OUTPUT
mkdir RESOURCES
mkdir MISC
### remove gcta output directory
rm -R output/
### move main output csv files, summary figures, and nohup log files into OUTPUT/
if [ -f GPAS_OUTPUT.csv ]
then
      mv GPAS_OUTPUT.csv OUTPUT/
      mv GPAS_OUTPUT_STREAMLINED.csv OUTPUT/
      mv ABC_OPTIM_INPUT_SUMMSTATS.csv OUTPUT/
      mv ABC_OPTIM_OUTPUT-with_weights.csv OUTPUT/
      mv *.png OUTPUT/
      mv QTL_SPEC.csv OUTPUT/
      mv LANDSCAPE.stat OUTPUT/
      mv *.nohup OUTPUT/
      ### move the other files into MISC/
      ls | grep "[.]" | xargs mv -t MISC/
      ### move the *.bed, *.bim, and *.fam files into the RESOURCES/ directory
      mv MISC/*.bed RESOURCES/
      mv MISC/*.bim RESOURCES/
      mv MISC/*.fam RESOURCES/
      ### remove MISC/
      rm -R MISC/
      ### END
      echo "FINALLY DONE!"
else
      exit 1
fi