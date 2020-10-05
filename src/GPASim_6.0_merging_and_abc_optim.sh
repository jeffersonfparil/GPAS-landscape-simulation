#!/bin/bash

##########################################
###     GPAS PERFORMANCE MERGING       ###
### AND SAMPLING STRATEGY OPTIMIZATION ###
##########################################

### Input and output specifications:
while [ ! -z "$1" ]; do
  case "$1" in
    --dir|-d) shift
      if [ $(ls $1 | wc -l) == 0 ]; then
        echo "Directory is empty: $1"
      else
        DIR=$1
      fi;;
    --src-dir|-s) shift
      if [ $(ls $1 | wc -l) == 0 ]; then
        echo "Directory is empty: $1"
      else
        SRC_DIR=$1
      fi;;
    --replicate-id|-x) shift
      rep=$1;;
    --n-individuals|-n) shift
      nIndividuals=$1;;
    --n-loci|-l) shift
      nLoci=$1;;
    --n-QTL|-q) shift
      nQTL=$1;;
    --n-BGS|-b) shift
      nBGS=$1;;
    --n-alleles|-a) shift
      nAlleles=$1;;
    --alleles-dist|-r) shift
      allele_eff_model=$1;;
    --n-generations|-t) shift
      nGen=$1;;
    --n-populations|-P) shift
      nPop=$1;;
    --migration-rate|-m) shift
      migration=$1;;
    --selection-intensity-QTL|-sq) shift
      selection=$1;;
    --selection-intensity-BGS|-sb) shift
      bg_selection=$1;;
    --diffusion-gradient|-d) shift
      GRADIENT=$1;;
    --n-pools|-np) shift
      NPOOLS=$1;;
    --max-sequencing|-ms) shift
      NSEQ_MAX=$1;;
    --n-libraries|-nl) shift
      NLIB=$1;;
    --help|-h) shift
      echo "###################################################################"
      echo "### GPAS PERFORMANCE MERGING AND SAMPLING STRATEGY OPTIMIZATION ###"
      echo "###################################################################"
      echo ""
      echo ""
      echo '##############'
      echo '### INPUTS ###'
      echo '##############'
      echo "--dir|-d                        directory containing all the ouputs of GPASim_4.0_filter_build_covariates.sh"
      echo "--src-dir|-s                    directory where plink, gemma, gcta and genomic_prediction/ are located"
      echo "--n-individuals|-n              number of individuals to simulate"
      echo "--n-loci|-l                     number of loci to simulate (neutral + QTL)"
      echo "--n-QTL|-q                      number of QTL among the all the loci simulated"
      echo "--n-BGS|-b                      number of background selection loci (should be highly polygenic >= 100)"
      echo "--n-alleles|-a                  number of alleles per loci, e.g. 5 for A,T,C,G, and DEL"
      echo "--alleles-dist|-r               probability distribution to sample the allele effects from (i.e. CHISQ, NORM, & UNIF)"
      echo "--n-generations|-t              number of generations to simulate"
      echo "--n-populations|-P              number of populations/subpopulations to simulate (NOTE: must have a natural number square-root)"
      echo "--migration-rate|-m             migration rate across the populations (surrently using the 2D stepping stone model see line 117) [@Wang2017]"
      echo "--selection-intensity-QTL|-sq   selection intensity (nQTL) defined as the slope of the directional selection logistic curve (Richards, 1959): ranges from -Inf (select for low phen) to +Inf (select for high phen)"
      echo "--selection-intensity-BGS|-sb   background selection intensity (nBGS) defined as the slope of the directional selection logistic curve (Richards, 1959): ranges from -Inf (select for low phen) to +Inf (select for high phen)"
      echo "--diffusion-gradient|-d         gradient of QTL (foreground and background selection) across the landscape: 0 for uniform, 1 for first row of populations, and 2 for first and last rows of populations are the sources of the the non-zero (non-wild type) alleles"
      echo "--n-pools|-p                    number of pools per population"
      echo "--max-sequencing|-ms            maximum number of multiplexed sequencing libraries to use (NOTE: The limiting factor!)"
      echo "--n-libraries|-l                maximum number of individually barcoded libraries for Indi-seq or Pool-seq"
      echo "--help|-h                       help documentation"
      echo ""
      echo ""
      echo '###############'
      echo '### OUTPUTS ###'
      echo '###############'
      echo "(1) GPAS_OUTPUT.csv (comma-separated; header: |col01:REPLICATE|col02:EFFECTIVE_POP_SIZE|col03:NPOP|col04:NGEN|col05:NLOCI|col06:NALLELES|col07:NQTL|col08:NBGS|col09:QTL_EFF_DIST|col10:SELECTION_INTENSITY|\ "
      echo "                                              |col11:BGS_INTENSITY|col12:MIGRATION_RATE|col13:QTL_GRADIENT|col14:PHENO_MEAN_BETADIST_SHAPES|col15:PHENO_VAR_BETADIST_SHAPES|col16:GENO_FST_BETADIST_SHAPES|col17:NPOOLS_PER_POP|col18:NLIB_MAXIMUM|col19:SAMPLING_SCHEME|col20:GENOTYPING_SCHEME|\ "
      echo "                                              |col21:TRAIN_POP|col22:TEST_POP|col23:NTRAIN|col24:NTEST|col25:MODEL|col26:COVARIATE|col27:p0|col28:p1|col29:r2adj|col30:AUC|\ "
      echo "                                              |col31:FRACTION_QTL_FIXED|col32:AUC_CORRECTED|col33:BONFERRONI_5PERCENT_TRUE_POSITIVE|col34:BONFERRONI_5PERCENT_FALSE_POSITIVE|col35:BONFERRONI_5PERCENT_QTL_ID|col36:RMSE|col37:CORRELATION| "
      echo "(2) GPAS_OUTPUT_STREAMLINED.csv (auto-cross-validation and data-points with only 2 pools in the training population are filtered out;"
      echo "                                 comma-separated; header: |col01:SAMPLING_SCHEME|col02:GENOTYPING_SCHEME|col03:TRAIN_POP|col04:TEST_POP|col05:NTRAIN|\ "
      echo "                                                          |col06:NTEST|col07:MODEL|col08:COVARIATE|col09:AUC|col10:AUC_CORRECTED|\ "
      echo "                                                          |col11:BONFERRONI_5PERCENT_TRUE_POSITIVE|col12:BONFERRONI_5PERCENT_FALSE_POSITIVE|col13:BONFERRONI_5PERCENT_QTL_ID|col14:RMSE|col15:CORRELATION|) "
      echo "(3) ABC_OPTIM_INPUT_SUMMSTATS.csv (sampling strategies across random combinatorial sensible samples and the corresponsing GPAS performance summaries;"
      echo "                                  comma-separated; header: |col01:ACROSS_INDI|col02:ACROSS_POOL|col03:WITHIN_INDI|col04:WITHIN_POOL|col05:AUC_MEAN|col06:AUC_CORRECTED_MEAN|col07:RMSE_MEAN|col08:CORRELATION_MEAN|col09:TRUE_POSITIVE_RATE|col10:FALSE_DISCOVERY_RATE|)"
      echo "(4) ABC_OPTIM_OUTPUT-with_weights.csv (ABC optimized sampling scheme partitioning (Across_Indi vs Across_Pool vs Within_Indi vs Within_Pool) given a range of weights (0.00 --> 1.00) to the 2 metrics: AUC_CORRECTED and RMSE;"
      echo "                                       for determining the saturation point of the metrics (AUC and RMSE) as their corresponding weights change ;"
      echo "                                       this implies a step-wise optimization of sampling strategy (how many experiments devoted to ACROSS_INDI vs ACROSS_POOL vs WITHIN_INDI vs WITHIN_POOL?):"
      echo "                                                (1) ABC optim across 0.0 to 0.5 weights range for AUC and RMSE (weight_auc + weight_rmse = 1)"
      echo "                                                (2) Determine the weight combination at which AUC and RMSE saturate (maximal AUC and minimal RMSE)"
      echo "                                                (3) Extract the resulting sampling strategy (experiment partitioning) at this optimal weight combination;"
      echo "                                       comma-separated; header: |col01:WEIGHT_AUC|col02:WEIGHT_RMSE|\ "
      echo "                                                                |col03:ACROSS_INDI|col04:ACROSS_POOL|col05:WITHIN_INDI|col06:WITHIN_POOL|\ "
      echo "                                                                |col07:AUC_MEAN|col08:AUC_CORRECTED_MEAN|col09:RMSE_MEAN|col10:CORRELATION_MEAN|col11:TRUE_POSITIVE_RATE|col12:FALSE_DISCOVERY_RATE|\ "
      echo "                                                                |col13:BZ_AUC (normalized then [0,1]-bounded)|col13:BZ_MINUS_ONE (normalized, then [0,1]-bounded, then subtracted from 1)|)"
      echo ""
      echo ""
      echo '###############'
      echo '### EXAMPLE ###'
      echo '###############'
      echo 'rep=1'
      echo 'nIndividuals=1000'
      echo 'nLoci=1000'
      echo 'nQTL=10'
      echo 'nBGS=10'
      echo 'nAlleles=2'
      echo 'allele_eff_model=CHISQ'
      echo 'nGen=150'
      echo 'nPop=16'
      echo 'migration=0.001'
      echo 'selection=0.10'
      echo 'bg_selection=0.00'
      echo 'GRADIENT=1'
      echo 'NPOOLS=5'
      echo 'NLIB=1000'
      echo 'NSEQ_MAX=10'
      echo 'DIR=/data/Lolium/Quantitative_Genetics/LOLSIM2020/LOLSIM_${rep}rep_${nQTL}qtl_${migration}mr_${selection}fgs_${bg_selection}bgs_${GRADIENT}grad'
      echo 'SRC_DIR=/data/Lolium/Softwares'
      echo 'time \'
      echo '${SRC_DIR}/genomic_prediction/src/GPASim_6.0_merging_and_abc_optim.sh \'
      echo '        --dir ${DIR} \'
      echo '        --src-dir ${SRC_DIR} \'
      echo '        --replicate-id ${rep} \'
      echo '        --n-individuals ${nIndividuals} \'
      echo '        --n-loci ${nLoci} \'
      echo '        --n-QTL ${nQTL} \'
      echo '        --n-BGS ${nBGS} \'
      echo '        --n-alleles ${nAlleles} \'
      echo '        --alleles-dist ${allele_eff_model} \'
      echo '        --n-generations ${nGen} \'
      echo '        --n-populations ${nPop} \'
      echo '        --migration-rate ${migration} \'
      echo '        --selection-intensity-QTL ${selection} \'
      echo '        --selection-intensity-BGS ${bg_selection} \'
      echo '        --diffusion-gradient ${GRADIENT} \'
      echo '        --n-pools ${NPOOLS} \'
      echo '        --max-sequencing ${NSEQ_MAX} \'
      echo '        --n-libraries ${NLIB}'
      echo ""
      echo ""
      exit 0
      ;;
    *)
      echo "What is this?! $1"
      exit 1
      ;;
  esac
shift
done
### Input parameter check
echo -e "--dir\n--src-dir\n--replicate-id\n--n-individuals\n--n-loci\n--n-QTL\n--n-BGS\n--n-alleles\n--alleles-dist\n--n-generations\n--n-populations\n--migration-rate\n--selection-intensity-QTL\n--selection-intensity-BGS\n--diffusion-gradient\n--n-pools\n--n-libraries\n--n-sequencing" > input_parameter_names.temp
echo -e "DIR\nSRC_DIR\nrep\nnIndividuals\nnLoci\nnQTL\nnBGS\nnAlleles\nallele_eff_model\nnGen\nnPop\nmigration\nselection\nbg_selection\nGRADIENT\nNPOOLS\nNLIB\nNSEQ" > input_parameter_values.temp
for i in $(seq 1 $(cat input_parameter_names.temp | wc -l))
do
  input_name=$(head -n${i} input_parameter_names.temp | tail -n1)
  input_value=$(head -n${i} input_parameter_values.temp | tail -n1)
  if [ $(echo $input_value | wc -w) == 0 ]
  then
    echo "Please provide an input for: $input_name"
    exit 1
  fi
done
rm input_parameter_names.temp input_parameter_values.temp

### TEST:
# rep=1                                                         # replicate number
# nIndividuals=1000                                             # number of individuals to simulate
# nLoci=1000                                                    # number of loci to simulate (neutral + QTL)
# nQTL=10                                                       # number of QTL among the all the loci simulated (should range from 1 to nLoci-1)
# nBGS=10                                                       # number of background selection loci (should range from 1 to nLoci-1; conditional on nQTL where nQTL + nBGS <= nLoci)
# nAlleles=2                                                    # number of alleles per loci, e.g. 5 for A,T,C,G, and DEL (NOTE: We're now assuming multi-allelic SNP loci are just noise so we're simulating only biallelic loci!)
# allele_eff_model=CHISQ                                        # probability distribution to sample the allele effects from (i.e. CHISQ, NORM, & UNIF)
# nGen=150                                                      # number of generations to simulate (approximating the number of generations since the introduction of Lolium in Australia in the 1880s minus the time when extensive herbicide application started)
# nPop=16                                                       # number of populations/subpopulations to simulate (NOTE: must have a natural number square-root)
# migration=0.001                                               # migration rate across the populations (surrently using the 1D stepping stone model see line 117)
# selection=0.10                                                # selection intensity (trait of interest) defined as the slope of the directional selection logistic curve (Richards, 1959): ranges from -Inf (select for low phen) to +Inf (select for high phen)
# bg_selection=0.10                                             # background selection intensity defined as the slope of the directional selection logistic curve (Richards, 1959): ranges from -Inf (select for low phen) to +Inf (select for high phen)
# GRADIENT=1                                                    # uniformly distributed non-wildtype alleles
# NPOOLS=5                                                      # number of pools to genotype per population (NOTE: for across population sample 1 population:1 pool)
# NSEQ_MAX=10                                                   # maximum number of multiplexed sequencing libraries to use (NOTE: The limiting factor!)"
# NLIB=1000                                                     # maximum number of individually barcoded libraries for sequencing (Indi-seq and Pool-seq) (NOTE: May need to have a minimum limit at 1,000 to keep GCTA from failing due to small population size leading to REML convergence failure)
# DIR=/data/Lolium/Quantitative_Genetics/LOLSIM2020/GPASim_${rep}REP_${nQTL}NQTL_${nBGS}NBGS_${migration}MIG_${selection}SEL_${bg_selection}BGSEL_${GRADIENT}GRAD

### Set working directories
cd $DIR

### Prepare the comma-separated file of landscape simulation parameters, i.e. "population_parameters.csv":
echo "rep,nIndividuals,nLoci,nQTL,nBGS,nAlleles,allele_eff_model,nGen,nPop,migration,selection,bg_selection,GRADIENT,NPOOLS,NLIB" \
      > population_parameters.csv
echo ${rep},${nIndividuals},${nLoci},${nQTL},${nBGS},${nAlleles},${allele_eff_model},${nGen},${nPop},${migration},${selection},${bg_selection},${GRADIENT},${NPOOLS},${NLIB} \
      >> population_parameters.csv

echo "####################################################################"
echo "MERGE AND STREAMLINE LANDSCAPE SIMULATION PARAMETERS AND GPAS OUTPUT"
echo "####################################################################"
Rscript ${SRC_DIR}/genomic_prediction/src/GPASim_6.1_merge.r \
                        population_parameters.csv \
                        LANDSCAPE.stat \
                        LANDSCAPE.fst \
                        OUTPUT_MERGED_AUC.csv \
                        OUTPUT_MERGED_RMSE.csv
### Inputs:
### (1) population_parameters.csv (landscape simulation parameters)
### (2) LANDSCAPE.stat (Phenotypic summary statistics per population; comma-separated; header: |col1:POP|col2:SIZE|col3:MEAN|col4:VAR|)
### (3) LANDSCAPE.fst (Pairwise fixation indices; comma-separated; headerless)
### (4) OUTPUT_MERGED_AUC.csv (GWAS performances; comma-separated; |col01:SAMPLE_ID|col02:SAMPLING|col03:POPULATION|col04:MODEL|col05:COVARAITE|col06:AUC given QTL after filtering|col07:FRACTION_QTL_FIXED|col08:BONFERRONI_5PERCENT_TRUE_POSITIVE|col09:BONFERRONI_5PERCENT_FALSE_POSITIVE|col10:BONFERRONI_5PERCENT_QTL_ID|)"
### (5) OUTPUT_MERGED_RMSE.csv (GP performances; comma-separated; header: |col1:TRAIN_NAME|col02:TEST_NAME|col03:TRAIN_MODEL|col04:TRAIN_COVAR|col05:N_TRAIN|col06:n_test|col07:p0|col08:P1|col09:R2ADJ|col10:RMSE|col11:CORRELATION|)"
### Outputs
### (1) GPAS_OUTPUT.csv (comma-separated; header: |col01:REPLICATE|col02:EFFECTIVE_POP_SIZE|col03:NPOP|col04:NGEN|col05:NLOCI|col06:NALLELES|col07:NQTL|col08:NBGS|col09:QTL_EFF_DIST|col10:SELECTION_INTENSITY|\
###                                               |col11:BGS_INTENSITY|col12:MIGRATION_RATE|col13:QTL_GRADIENT|col14:PHENO_MEAN_BETADIST_SHAPES|col15:PHENO_VAR_BETADIST_SHAPES|col16:GENO_FST_BETADIST_SHAPES|col17:NPOOLS_PER_POP|col18:NLIB_MAXIMUM|col19:SAMPLING_SCHEME|col20:GENOTYPING_SCHEME|\
###                                               |col21:TRAIN_POP|col22:TEST_POP|col23:NTRAIN|col24:NTEST|col25:MODEL|col26:COVARIATE|col27:p0|col28:p1|col29:r2adj|col30:AUC|\
###                                               |col31:FRACTION_QTL_FIXED|col32:AUC_CORRECTED|col33:BONFERRONI_5PERCENT_TRUE_POSITIVE|col34:BONFERRONI_5PERCENT_FALSE_POSITIVE|col35:BONFERRONI_5PERCENT_QTL_ID|col36:RMSE|col37:CORRELATION|
### (2) GPAS_OUTPUT_STREAMLINED.csv (comma-separated; header: |col01:SAMPLING_SCHEME|col02:GENOTYPING_SCHEME|col03:TRAIN_POP|col04:TEST_POP|col05:NTRAIN|\
###                                                           |col06:NTEST|col07:MODEL|col08:COVARIATE|col09:AUC|col10:AUC_CORRECTED|\
###             maximum number of multiplexed sequencing libraries to use (NOTE: The limiting factor!)                                              |col11:BONFERRONI_5PERCENT_TRUE_POSITIVE|col12:BONFERRONI_5PERCENT_FALSE_POSITIVE|col13:BONFERRONI_5PERCENT_QTL_ID|col14:RMSE|col15:CORRELATION|)

echo "################################"
echo "BUILD THE ABC OPTIMIZATION INPUT"
echo "################################"
###   specifically the summary statistics of the GPAS performance metrics
###   across different combinations of the 4 sampling strategies:
###       (1) Across population sampling with individual genotyping (Indi-seq)
###       (2) Across population sampling with pool genotyping (Pool-seq)
###       (3) Within population sampling with Indi-seq
###       (4) Within population sampling with Pool-seq
### NOTE: Each instance of a sampling strategy combination has a fixed number of sampling-x-genotyping samples,
###       i.e. the number of sequencing experiments (the limiting factor or resource; NSEQ =  1 --> total number of populations)

### combinatorics, sampling and summary statistics (loopy section in parallel)
echo "using Distributed
      Distributed.addprocs(length(Sys.cpu_info()))
      @everywhere include(\"${SRC_DIR}/genomic_prediction/src/GPASim_6.3_abc_input.jl\")
      random_combinatorics_and_summary_statistics(ARGS[1], parse(Int64,ARGS[2]), parse(Int64,ARGS[3]), \"${SRC_DIR}/genomic_prediction/src/\", false)
" > build_abc_optimization_input.jl
time \
julia build_abc_optimization_input.jl \
                  GPAS_OUTPUT_STREAMLINED.csv \
                  ${NSEQ_MAX} \
                  ${nQTL}
### Inputs:
### (1) GPAS_OUTPUT_STREAMLINED.csv (comma-separated; header: |col01:SAMPLING_SCHEME|col02:GENOTYPING_SCHEME|col03:TRAIN_POP|col04:TEST_POP|col05:NTRAIN|\
###                                                           |col06:NTEST|col07:MODEL|col08:COVARIATE|col09:AUC|col10:AUC_CORRECTED|\
###                                                           |col11:BONFERRONI_5PERCENT_TRUE_POSITIVE|col12:BONFERRONI_5PERCENT_FALSE_POSITIVE|col13:BONFERRONI_5PERCENT_QTL_ID|col14:RMSE|col15:CORRELATION|)
### (2) maximum number of multiplexed sequencing libraries to use (NOTE: The limiting factor!)
### (3) number of simulated QTL
### Outputs:
### (1) ABC_OPTIM_INPUT_SUMMSTATS.csv (sampling strategies across random combinatorial sensible samples and the corresponsing GPAS performance summaries;
###                                  comma-separated; header: |col01:ACROSS_INDI|col02:ACROSS_POOL|col03:WITHIN_INDI|col04:WITHIN_POOL|col05:AUC_MEAN|col06:AUC_CORRECTED_MEAN|col07:RMSE_MEAN|col08:CORRELATION_MEAN|col09:TRUE_POSITIVE_RATE|col10:FALSE_DISCOVERY_RATE|)
### (2) POPULATION_COMBINATIONS_*.csv (output of the adaptive function in "population_combinatorics.jl"; comma-separated; headerless)
### (3) SUMMSTATS_*.csv (GPAS performance summary across random combinatorial sensible samples; comma-separated; headerless)

### clean-up (too many to remove for a single rm * argument therefore using xars with input redirecton from a file ls-grepped)
ls | grep "POPULATION_COMBINATIONS_" > remove_us.temp
ls | grep "SUMMSTATS_" >> remove_us.temp
cat remove_us.temp | xargs rm


echo "################"
echo "ABC OPTIMIZATION"
echo "################"
Rscript ${SRC_DIR}/genomic_prediction/src/GPASim_6.4_abc_with_weight.r \
                                                    ABC_OPTIM_INPUT_SUMMSTATS.csv
### Input:
### (1) ABC_OPTIM_INPUT_SUMMSTATS.csv (sampling strategies across random combinatorial sensible samples and the corresponsing GPAS performance summaries;
###                                    comma-separated; header: |col01:ACROSS_INDI|col02:ACROSS_POOL|col03:WITHIN_INDI|col04:WITHIN_POOL|col05:AUC_MEAN|col06:AUC_CORRECTED_MEAN|col07:RMSE_MEAN|col08:CORRELATION_MEAN|col09:TRUE_POSITIVE_RATE|col10:FALSE_DISCOVERY_RATE|)
### Outputs:
### (1) ABC_OPTIM_OUTPUT-with_weights.csv (ABC optimized sampling scheme partitioning (Across_Indi vs Across_Pool vs Within_Indi vs Within_Pool) given a range of weights (0.00 --> 1.00) to the 2 metrics: AUC_CORRECTED and RMSE;
###                                        for determining the saturation point of the metrics (AUC and RMSE) as their corresponding weights change ;
###                                        this implies a step-wise optimization of sampling strategy (how many experiments devoted to ACROSS_INDI vs ACROSS_POOL vs WITHIN_INDI vs WITHIN_POOL?):
###                                                 (1) ABC optim across 0.0 to 0.5 weights range for AUC and RMSE (weight_auc + weight_rmse = 1)
###                                                 (2) Determine the weight combination at which AUC and RMSE saturate (maximal AUC and minimal RMSE)
###                                                 (3) Extract the resulting sampling strategy (experiment partitioning) at this optimal weight combination;
###                                        comma-separated; header: |col01:WEIGHT_AUC|col02:WEIGHT_RMSE|\
###                                                                 |col03:ACROSS_INDI|col04:ACROSS_POOL|col05:WITHIN_INDI|col06:WITHIN_POOL|\
###                                                                 |col07:AUC_MEAN|col08:AUC_CORRECTED_MEAN|col09:RMSE_MEAN|col10:CORRELATION_MEAN|col11:TRUE_POSITIVE_RATE|col12:FALSE_DISCOVERY_RATE|\
###                                                                 |col13:BZ_AUC (normalized then [0,1]-bounded)|col13:BZ_MINUS_ONE (normalized, then [0,1]-bounded, then subtracted from 1)|)
### (2) ABC_OPTIM_OUTPUT-with_weights-*_vs_*.png (weights vs metrics and params vs metrics scatter plots)
### (3) ABC_OPTIM_PRELIM-*.png (preliminary sanity checking plots)
