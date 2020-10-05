#!/bin/bash

###################################################################
### GENOMIC PREDICTION AND GENOME-WIDE ASSOCIATION STUDY (GPAS) ###
###################################################################

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
    --n-threads|-t) shift
      NTHREADS=$1;;
    --help|-h) shift
      echo "###################################################################"
      echo "### GENOMIC PREDICTION AND GENOME-WIDE ASSOCIATION STUDY (GPAS) ###"
      echo "###################################################################"
      echo ""
      echo ""
      echo '##############'
      echo '### INPUTS ###'
      echo '##############'
      echo "--dir|-d          directory containing all the ouputs of GPASim_4.0_filter_build_covariates.sh"
      echo "--src-dir|-s      directory where plink, gemma, gcta and genomic_prediction/ are located"
      echo "--n-threads|-t    number of threads to use in parallel"
      echo "--help|-h         help documentation"
      echo ""
      echo ""
      echo '###############'
      echo '### OUTPUTS ###'
      echo '###############'
      echo "(1) OUTPUT_MERGED_AUC.csv (GWAS performances; comma-separated; |col01:SAMPLE_ID|col02:SAMPLING|col03:POPULATION|col04:MODEL|col05:COVARAITE|col06:AUC given QTL after filtering|col07:FRACTION_QTL_FIXED|col08:BONFERRONI_5PERCENT_TRUE_POSITIVE|col09:BONFERRONI_5PERCENT_FALSE_POSITIVE|col10:BONFERRONI_5PERCENT_QTL_ID|)"
      echo "(2) OUTPUT_MERGED_RMSE.csv (GP performances; comma-separated; header: |col1:TRAIN_NAME|col02:TEST_NAME|col03:TRAIN_MODEL|col04:TRAIN_COVAR|col05:N_TRAIN|col06:n_test|col07:p0|col08:P1|col09:R2ADJ|col10:RMSE|col11:CORRELATION|)"
      echo "(3) Intermediate outputs per sample (for GWAS) and per sample-pair (for GP)"
      echo "    - *.gwas (GWAS ouput; tab-delimited for Indi-GWAS and comma-separated (for Pool-GWAS)"
      echo "    - *.ranef (Pool-GWAS random effects; comma-separated"
      echo "    - *.gp (GP output:"
      echo "            + Indi-GP - tab-delimited; headerless: |col1:FAMID|col2:INDID|col3:PHENO|col4:NREF|col5:NALT|col6:SCORE|"
      echo "            + Pool-GP - tab-delimited; headerless: |col1:FAMID|col2:INDID|col3:PHENO|col4:NREF|col5:NALT|col6:SCORE|)"
      echo "    - *.p0p1 (linear link from polygenic score into phenotypic values:"
      echo "            + tab-delimited"
      echo "            + filename of the predictors of the y_pred~polygenic linear model; header: |col1:p0|col2:p1|col3:r2adj|"
      echo "            + if and only if: ${TRAIN_prefix}==${TEST_prefix} that is the second step of the 2-step prediction model building: y_pred = p0 + p1*plink_polygenic_score; where plink_polygenic_score was derived using GWAS output"
      echo "    - *.auc (GWAS performance statistics; tab-delimited; header: |col1:AUC|col2:FRACTION_OF_FIXED_QTL|)"
      echo "    - *.rmse (GP performance statistics; tab-delimited; header: |col1:train_name|col02:test_name|col03:train_model|col04:train_covar|col05:n_train|col06:n_test|col07:p0|col08:p1|col09:r2adj|col10:RMSE|col11:correlation|)"
      echo ""
      echo ""
      echo '###############'
      echo '### EXAMPLE ###'
      echo '###############'
      echo 'rep=1'
      echo 'nQTL=10'
      echo 'migration=0.00'
      echo 'selection=0.25'
      echo 'bg_selection=0.00'
      echo 'GRADIENT=1'
      echo 'DIR=/data/Lolium/Quantitative_Genetics/LOLSIM2020/LOLSIM_${rep}rep_${nQTL}qtl_${migration}mr_${selection}fgs_${bg_selection}bgs_${GRADIENT}grad'
      echo 'SRC_DIR=/data/Lolium/Softwares'
      echo 'NTHREADS=12'
      echo 'time \'
      echo '${SRC_DIR}/genomic_prediction/src/GPASim_5.0_GPAS.sh \'
      echo '      -d $DIR \'
      echo '      -s $SRC_DIR \'
      echo '      -t $NTHREADS'
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
echo -e "--dir\n--src-dir\n--n-threads" > input_parameter_names.temp
echo -e "$DIR\n$SRC_DIR\n$NTHREADS" > input_parameter_values.temp
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
### Set the amount of memory to allocate for each instance of plink execution
PLINK_MEM=$(echo $(head /proc/meminfo | grep "MemAvailable" | cut -d: -f2 | sed "s/[ ]//g" | sed "s/kB//g") / ${NTHREADS}000 | bc)
### Set working directories
cd $DIR
echo "######################"
echo "### Model building ###"
echo "######################"
### NOTE: - build iterative and non-iterative linear models
###       - perform significance tests to identify QTL
###       - performance metric: ROC plot AUC (area under the false positive rate vs true positive rate curve)
echo "#@@@@@@@@@@"
echo "# Indi-GWAS"
echo "#@@@@@@@@@@"
ls POP_*.bed MULTIPOP_*_*.bed  | sort | uniq | sed 's/.bed//g' > INDISEQ_PREFIXES.txt
time \
parallel -j ${NTHREADS} ${SRC_DIR}/genomic_prediction/src/GPASim_5.1_Indi-GWAS.sh \
                    {1} \
                    {1}_{2}.kin \
                    1 \
                    ${PLINK_MEM} \
                    ${SRC_DIR} \
                    ::: $(cat INDISEQ_PREFIXES.txt) \
                    ::: GRM STANDARDIZED
### Inputs:
### (1) prefix of the bed-bim-fam plink file of genotype and phenotype data (POP_* or MULTIPOP_*_*)
### (2) kinship file to use (headerless, tab-delimited, square-symmetric matrix) (GRM or STANDARDIZED)
###     NOTE: for GCTA the kinship matrix will be set as fixed effect covariate on top of ${prefix}_GRM_sparse*
### (3) number of compute cores or threads to use for parallel computation
### (4) path to the softwares directory containing plink, emmax, gcta, and genomic_prediction.git
### Outputs:
### (1) ${prefix}_EMMAX_${kinship_id}.gwas ### tab-delimited; headerless; |col1:SNP|col2:EFF|col3:SE|col4:PVAL|
### (2) ${prefix}_GEMMA_${kinship_id}.gwas ### tab-delimited; HEADER: |col01:chr|col02:rs|col03:n_miss|col04:allele1|col05:allele0|col06:af|col07:beta|col08:se|col09:logl_H1|col10:l_reml|col11:l_mle|col12:p_wald|col13:p_lrt|col14:p_score|
### (3) ${prefix}_GCTA_${kinship_id}.gwas ### for GCTA's LR without GCTA's GRM binary input ### tab-delimited; HEADER: |col01:CHR|col02:SNP|col03:POS|col04:A1|col05:A2|col06:N|col07:AF1|col08:BETA|col09:SE|col10:P|
### where:
### xiferp=$(echo $prefix | rev); pihsnik=$(echo ${kinship%.*} | rev)
### kinship_id=$(echo ${pihsnik%_$(echo $xiferp)*} | rev)
echo "#@@@@@@@@@@"
echo "# Pool-GWAS"
echo "#@@@@@@@@@@"
### NOTE: cv.glmnet fails for n_pools=2
### CONSEQUENCE: All Pool-GPAS models that are non-iterative failed! i.e. REML_[LS, RR, GLMNET, LASSO]_[HIVERT, WEIRCOCK]
ls POP_*.sync MULTIPOP_*_*.sync  | sort | uniq | sed 's/.sync//g' > POOLSEQ_PREFIXES.txt
echo '
  #!/bin/bash
  prefix=${1}
  SRC_DIR=${2}
  julia ${SRC_DIR}/genomic_prediction/src/GPASim_5.2_Pool-GWAS.jl ${prefix}
' > PoolGWAS.sh
chmod +x PoolGWAS.sh
time \
parallel -j ${NTHREADS} ./PoolGWAS.sh \
                      {} \
                      ${SRC_DIR} \
                      ::: $(cat POOLSEQ_PREFIXES.txt)
### Input:
### (1) filename of the training population's output of genome-wide association study using pool sequencing data (Pool-GWAS)
###     - string(prefix, "_GWAlpha_SNPwise.gwas")
###     - string(prefix, "_GWAlpha_MIXED_REML_", fst_id, ".gwas")
###     - string(prefix, "_GWAlpha_MIXED_REML_", fst_id, ".ranef")
###     - string(prefix, "_GWAlpha_RIDGE.gwas")
###     - string(prefix, "_GWAlpha_LASSO.gwas")
###     - string(prefix, "_GWAlpha_GLMNET.gwas")
### where:
# xiferp=$(echo $prefix | rev); tsf=$(echo ${fst%.*} | rev)
# fst_id=$(echo ${tsf%_$(echo $xiferp)*} | rev) ### (HIVERT or WEIRCOCK)
# files are comma-separated
# *.gwas: HEADER: |col1:CHROM|col2:POS|col3:ALLELES|col4:FREQ|col5:ALPHA|col6:PVALUES|col7:LOD|
# *.ranef: HEADERLESS: |col1:random effects|

echo "############################"
echo "### Phenotype Prediction ###"
echo "############################"
### NOTE: - use iterative and non-iterative linear models to predict phenotypes
###       - cross-validate on a test population
###       - performance metric: RMSE (root mean square error: sqrt[mean((y_pred - y_true)^2)])
echo "#@@@@@@"
echo "Indi-GP"
echo "#@@@@@@"
### Two-step genomic prediction using Indi-GWAS output
### Step 1: Calculate the polygenic "risk" score using test population's genotype data (.bed,.bim,&.fam) and training populations's SNP effects (.gwas)
### Step 2: Predict the test population's phenotypes using the training population's linear model connecting the actual phenotypes with plink-derived polygenic "risk" scores, i.e.:
###             - y_train ~ po + p1*polygenic_score_train
###             - y_test_predicted = p0 + p1*polygenic_score_test
### NOTE: Execute training population == test populations initially to compute p0 and p1
### auto-cross-validation
time \
parallel -j ${NTHREADS} ${SRC_DIR}/genomic_prediction/src/GPASim_5.3_Indi-GP.sh \
                      {1}_{2}_{3}.gwas \
                      {1} \
                      ${PLINK_MEM} \
                      ${SRC_DIR} \
                      ::: $(cat INDISEQ_PREFIXES.txt) \
                      ::: EMMAX GEMMA GCTA \
                      ::: GRM STANDARDIZED
### no-auto-cross-validation
### i.e. cross-validation such that training != test or validation population
### Iterating across testing/validataion populations to avoid conflicts of two different instances writing on the same file i.e. train==test and test==train
### ID mismatches in the --score file is due to the headers of the GWAS output file
echo '#!/bin/bash
  if [ ${1} == ${4} ]
  then
    exit
  fi
  ${6}/genomic_prediction/src/GPASim_5.3_Indi-GP.sh \
                        ${1}_${2}_${3}.gwas \
                        ${4} \
                        ${5} \
                        ${6}
' > IndiGP_noauto.sh
chmod +x IndiGP_noauto.sh
time \
parallel -j ${NTHREADS} ./IndiGP_noauto.sh \
                      {1} \
                      {2} \
                      {3} \
                      {4} \
                      ${PLINK_MEM} \
                      ${SRC_DIR} \
                      ::: $(cat INDISEQ_PREFIXES.txt) \
                      ::: EMMAX GEMMA GCTA \
                      ::: GRM STANDARDIZED \
                      ::: $(cat INDISEQ_PREFIXES.txt)
### Inputs:
### (1) filename of the training population's output of genome-wide association study using individual sequencing data (Indi-GWAS)
###     - ${TRAIN_prefix}_EMMAX_${TRAIN_kinship_id}.gwas ### tab-delimited; headerless; |col1:SNP|col2:EFF|col3:SE|col4:PVAL|
###     - ${TRAIN_prefix}_GEMMA_${TRAIN_kinship_id}.gwas ### tab-delimited; HEADER: |col01:chr|col02:rs|col03:n_miss|col04:allele1|col05:allele0|col06:af|col07:beta|col08:se|col09:logl_H1|col10:l_reml|col11:l_mle|col12:p_wald|col13:p_lrt|col14:p_score|
###     - ${TRAIN_prefix}_GCTA_${TRAIN_kinship_id}.gwas ### for GCTA's LR without GCTA's GRM binary input ### tab-delimited; HEADER: |col01:CHR|col02:SNP|col03:POS|col04:A1|col05:A2|col06:N|col07:AF1|col08:BETA|col09:SE|col10:P|
### (2) prefix of the test or validation population
###     - NOTE: We will be using the FULLSNPSET_${TEST_prefix}.bed, FULLSNPSET_${TEST_prefix}.bim, and FULLSNPSET_${TEST_prefix}.fam dataset
### (3) path to the softwares directory containing plink, emmax, gcta, and genomic_prediction.git
### Output:
### (1) ${TRAIN_prefix}_${TRAIN_model}_${TRAIN_kinship_id}-${TEST_prefix}.gp (tab-delimited; headerless: |col1:FAMID|col2:INDID|col3:PHENO|col4:NREF|col5:NALT|col6:SCORE|)

echo "#@@@@@@@"
echo "Pool-GP"
echo "#@@@@@@@"
### NOTE: during model building above, cv.glmnet failed for n_pools=2
### CONSEQUENCE: All Pool-GPAS models that are non-iterative failed! i.e. REML_[LS, RR, GLMNET, LASSO]_[HIVERT, WEIRCOCK]
### PARSE SYNC INTO CSV OF ALLELE FREQUENCIES
echo "using DelimitedFiles
      using Distributed
      Distributed.addprocs(length(Sys.cpu_info()))
      @everywhere using GWAlpha
      poolseq_prefixes = DelimitedFiles.readdlm(\"POOLSEQ_PREFIXES.txt\")[:,1]
      @time @sync @distributed for prefix in poolseq_prefixes
        fname_sync = string(\"FULLSNPSET_\", prefix, \".sync\")
        GWAlpha.sync_processing_module.sync_parse(fname_sync)
      end
" > parse_sync_parallel.jl
julia parse_sync_parallel.jl
### PREDICT
time parallel Rscript ${SRC_DIR}/genomic_prediction/src/GPASim_5.4_Pool-GP.r {1} FULLSNPSET_{2}_ALLELEFREQ.csv \
              ::: $(ls | grep ".gwas" | grep "GWAlpha") \
              ::: $(cat POOLSEQ_PREFIXES.txt)
### Inputs:
### (1) filename of the training population's output of genome-wide association study using pool sequencing data (Pool-GWAS)
###     - string(prefix, "_GWAlpha_SNPwise.gwas")
###     - string(prefix, "_GWAlpha_MIXED_REML_", fst_id, ".gwas")
###     - string(prefix, "_GWAlpha_MIXED_REML_", fst_id, ".ranef")
###     - string(prefix, "_GWAlpha_RIDGE.gwas")
###     - string(prefix, "_GWAlpha_LASSO.gwas")
###     - string(prefix, "_GWAlpha_GLMNET.gwas")
###     - *.gwas: HEADER: |col1:CHROM|col2:POS|col3:ALLELES|col4:FREQ|col5:ALPHA|col6:PVALUES|col7:LOD|
###     - *.ranef: HEADERLESS: |col1:random effects|
### (2) filename of the test population's full SNP data in sync format (FULLSNPSET_*.sync)
### Output:
### (1) Genomic prection cross-validation statistics (comma-separated; |col01:TRAIN_name|col02:TRAIN_model|col03:TRAIN_rancovar|col04:TEST_name|col05:n_train|col06:n_test|col07:p0|col08:p1|col09:r2adj|co10l:y_true|col11:y_pred|)
###     - For GWAlpha_SNPwise: string(TRAIN_name, "_", TRAIN_model, "-", TEST_name, ".gp")
###     - For GWAlpha_REML_[LS, RR, GLMNET, LASSO]_[HIVERT, WEIRCOCK]: string(TRAIN_name, "_", TRAIN_model, "_", TRAIN_rancovar, "-", TEST_name, ".gp")


echo "###############################"
echo "### GPAS PERFORMACE METRICS ###"
echo "###############################"
echo "#@@@@@@@@@@@@@@@"
echo "GWAS Metric: AUC"
echo "#@@@@@@@@@@@@@@@"
### (1) AUC (Area under the FPR vs TPR (receivier operating characteristic) curve) - conditional on the number of QTL captured that is not including the QTL which were filtered out because they were ~fixed.
### (2) AUC * (1.00 - Fraction of Fixed or Filtered-out QTL)
echo '#!/bin/bash
  f=$1
  SRC_DIR=$2
  Rscript ${SRC_DIR}/genomic_prediction/src/GPASim_5.5_AUC.r \
                          $f \
                          QTL_SPEC.csv \
                          GENOME_SPEC.csv \
                          1 \
                          FALSE
  ### Inputs:
  ### (1) filename of the GWAS output
  ### (2) filename of the QTL file specification
  ### (3) filenames of the simulated genome specification
  ### (4) linkage block size assumed constant across the genome (in kilobases)
  ### (5) save manhattan and ROC plots as png files (i.e. TRUE or FALSE)
  ### Output:
  ### (1) ${prefix}_${model}_${kinship_id}.auc ### tab-delimited; HEADER: |col1:AUC|col2:FRACTION_OF_FIXED_QTL|col3:BONFERRONI_5PERCENT_TRUE_POSITIVE|col4:BONFERRONI_5PERCENT_FALSE_POSITIVE|col5:BONFERRONI_5PERCENT_QTL_ID|
  ### where:
  ### model=[EMMAX, GCTA, GEMMA]
  ### xiferp=$(echo $prefix | rev); pihsnik=$(echo ${kinship%.*} | rev)
  ### kinship_id=$(echo ${pihsnik%_$(echo $xiferp)*} | rev)
' > AUC.sh
chmod +x AUC.sh
ls | grep ".gwas" > gwas_list.temp
time \
parallel ./AUC.sh {} ${SRC_DIR} ::: $(cat gwas_list.temp)
### concatenate GWAS metrics
ls | grep ".auc" > AUC_OUTPUT_FNAMES.txt
ls | grep ".temp" | xargs rm
touch SAMPLE_ID.temp; touch SAMPLING.temp; touch POPULATION.temp
touch MODEL.temp; touch COVARIATE.temp
touch AUC.temp
touch PFIXED_QTL.temp
touch BONFERRONI_5PERCENT_TRUE_POSITIVE.temp
touch BONFERRONI_5PERCENT_FALSE_POSITIVE.temp
touch BONFERRONI_5PERCENT_QTL_ID.temp
for i in $(seq 1 $(cat AUC_OUTPUT_FNAMES.txt | wc -l))
do
  # echo $i
  COL1=$(cut -d_ -f1 AUC_OUTPUT_FNAMES.txt | head -n${i} | tail -n1)
  if [ $COL1 == "MULTIPOP" ]
  then
    echo "ACROSS" >> SAMPLING.temp
    cut -d_ -f1-3 AUC_OUTPUT_FNAMES.txt | head -n${i} | tail -n1 >> SAMPLE_ID.temp
    cat $(cut -d_ -f1-3 AUC_OUTPUT_FNAMES.txt | head -n${i} | tail -n1).idx | cut -f1 | sed ':a;N;$!ba;s/\n/:/g' >> POPULATION.temp
    MODEL_COVAR=$(cut -d_ -f4- AUC_OUTPUT_FNAMES.txt | head -n${i} | tail -n1 | sed 's/.auc//g')
  else
    echo "WITHIN" >> SAMPLING.temp
    cut -d_ -f1-2 AUC_OUTPUT_FNAMES.txt | head -n${i} | tail -n1 >> SAMPLE_ID.temp
    cut -d_ -f1-2 AUC_OUTPUT_FNAMES.txt | head -n${i} | tail -n1 >> POPULATION.temp
    MODEL_COVAR=$(cut -d_ -f3- AUC_OUTPUT_FNAMES.txt | head -n${i} | tail -n1 | sed 's/.auc//g')
  fi
  if [ $(echo $MODEL_COVAR | grep "GWAlpha" | wc -l) -eq 1 ]
  then
    echo $MODEL_COVAR | cut -d_ -f1-3 >> MODEL.temp
    echo $MODEL_COVAR | cut -d_ -f4 >> COVARIATE.temp
  else
    echo $MODEL_COVAR | cut -d_ -f1 >> MODEL.temp
    echo $MODEL_COVAR | cut -d_ -f2 >> COVARIATE.temp
  fi
  tail -n+2 $(head -n${i} AUC_OUTPUT_FNAMES.txt | tail -n1) | cut -f1 >> AUC.temp
  tail -n+2 $(head -n${i} AUC_OUTPUT_FNAMES.txt | tail -n1) | cut -f2 >> PFIXED_QTL.temp
  tail -n+2 $(head -n${i} AUC_OUTPUT_FNAMES.txt | tail -n1) | cut -f3 >> BONFERRONI_5PERCENT_TRUE_POSITIVE.temp
  tail -n+2 $(head -n${i} AUC_OUTPUT_FNAMES.txt | tail -n1) | cut -f4 >> BONFERRONI_5PERCENT_FALSE_POSITIVE.temp
  tail -n+2 $(head -n${i} AUC_OUTPUT_FNAMES.txt | tail -n1) | cut -f5 >> BONFERRONI_5PERCENT_QTL_ID.temp
done
echo -e "SAMPLE_ID,SAMPLING,POPULATION,MODEL,COVARIATE,AUC,FRACTION_QTL_FIXED,BONFERRONI_5PERCENT_TRUE_POSITIVE,BONFERRONI_5PERCENT_FALSE_POSITIVE,BONFERRONI_5PERCENT_QTL_ID" \
          > OUTPUT_MERGED_AUC.csv
paste -d, SAMPLE_ID.temp \
          SAMPLING.temp \
          POPULATION.temp \
          MODEL.temp \
          COVARIATE.temp \
          AUC.temp \
          PFIXED_QTL.temp \
          BONFERRONI_5PERCENT_TRUE_POSITIVE.temp \
          BONFERRONI_5PERCENT_FALSE_POSITIVE.temp \
          BONFERRONI_5PERCENT_QTL_ID.temp \
          >> OUTPUT_MERGED_AUC.csv
ls | grep ".temp" | xargs rm
### Input:
### (1) GWAS output data (*.gwas; tab/comma-separated)
### Output:
### (1) OUTPUT_MERGED_AUC.txt (comma-separated; |col01:SAMPLE_ID|col02:SAMPLING|col03:POPULATION|col04:MODEL|col05:COVARAITE|col06:AUC given QTL after filtering|col07:FRACTION_QTL_FIXED|col08:BONFERRONI_5PERCENT_TRUE_POSITIVE|col09:BONFERRONI_5PERCENT_FALSE_POSITIVE|col10:BONFERRONI_5PERCENT_QTL_ID|)

echo "#@@@@@@@@@@@@@@@"
echo "GP Metric: RMSE"
echo "#@@@@@@@@@@@@@@@"
echo '#!/bin/bash
  fname_gp=$1
  autocross=$2
  SRC_DIR=$3
  # fname_gp=POP_01_EMMAX_GRM-POP_01.gp
  # fname_gp=POP_01_EMMAX_GRM-POP_02.gp
  # autocross=TRUE
  # autocross=FALSE
  SAMPLING_ID=$(cut -d_ -f1 <<<$(echo $fname_gp))
  if [ $SAMPLING_ID == "POP" ]
  then
    TRAIN_prefix=$(cut -d_ -f1-2 <<<$(echo $fname_gp))
  elif [ $SAMPLING_ID == "MULTIPOP" ]
  then
    TRAIN_prefix=$(cut -d_ -f1-3 <<<$(echo $fname_gp))
  else
    echo "INVALID SAMPLING ID: " ${SAMPLING_ID}
  fi
  TEST_prefix=$(cut -d- -f2 <<<$(echo $fname_gp) | sed s/.gp//g)
  if [ $autocross == "TRUE" ]
  then
    if [ $TRAIN_prefix != $TEST_prefix ]
    then
      exit
    fi
  else
    if [ $TRAIN_prefix == $TEST_prefix ]
    then
      exit
    fi
  fi
  Rscript ${SRC_DIR}/genomic_prediction/src/GPASim_5.6_RMSE.r ${fname_gp}
' > rmse.sh
chmod +x rmse.sh

### Indi-GP
time xargs parallel ./rmse.sh {} TRUE ${SRC_DIR} ::: <<<$(ls | grep ".gp" | grep -v "GWAlpha") ### auto-cross-validation
time xargs parallel ./rmse.sh {} FALSE ${SRC_DIR} ::: <<<$(ls | grep ".gp" | grep -v "GWAlpha") ### no auto-cross-validation
### Pool-GP
time xargs parallel ./rmse.sh {} TRUE ${SRC_DIR} ::: <<<$(ls | grep ".gp" | grep "GWAlpha") ### no need to perform auto-cross-validation
time xargs parallel ./rmse.sh {} FALSE ${SRC_DIR} ::: <<<$(ls | grep ".gp" | grep "GWAlpha") ### but I'd like to re-use the bash code above for simplicity
### Merge RMSE output and other GP statistics
head -n1 $(ls | grep ".rmse" | head -n1) | sed 's/[[:blank:]]/,/g' > OUTPUT_MERGED_RMSE.csv
for rmse in $(ls | grep ".rmse")
do
  tail -n+2 ${rmse} | sed 's/[[:blank:]]/,/g' >> OUTPUT_MERGED_RMSE.csv
done
### Input:
### (1) ${UNIQUE_ID}.gp - filename of polygenic risk score summing-up output from plink --score
### Outputs:
### (1) ${UNIQUE_ID}.p0p1 - filename of the predictors of the y_pred~polygenic linear model; header: |col1:p0|col2:p1|col3:r2adj|
###     - if and only if: ${TRAIN_prefix}==${TEST_prefix} that is the second step of the 2-step prediction model building: y_pred = p0 + p1*plink_polygenic_score; where plink_polygenic_score was derived using GWAS output
### (2) ${UNIQUE_ID}.rmse - filename of the RSME output; header: |col1:TRAIN_NAME|col02:TEST_NAME|col03:TRAIN_MODEL|col04:TRAIN_COVAR|col05:N_TRAIN|col06:N_TEST|col07:P0|col08:P1|col09:R2ADJ|col10:RMSE|col11:CORRELATION|
### where:
### UNIQUE_ID=${TRAIN_prefix}_${TRAIN_model}_${TRAIN_kinship_id}-${TEST_prefix} for Indi-GP; and
### UNIQUE_ID=string(TRAIN_name, "_", TRAIN_model, "_", TRAIN_rancovar, "-", TEST_name) for Pool-GP
### NOTE: These were merged into a single file!
echo "#########################################################################################"
