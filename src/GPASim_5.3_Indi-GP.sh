#!/bin/bash
###################################################
###          GENOMIC PREDICTION (GP)            ###
### using Individual sequencing (Indi-seq) data ###
###################################################
### 2-step genomic prediction using Indi-GWAS output
### Step 1: Calculate the polygenic "risk" score using test population's genotype data (.bed,.bim,&.fam) and training populations's SNP effects (.gwas)
### Step 2: Predict the test population's phenotypes using the training population's linear model connecting the actual phenotypes with plink-derived polygenic "risk" scores, i.e.:
###             - y_train ~ po + p1*polygenic_score_train
###             - y_test_predicted = p0 + p1*polygenic_score_test
### NOTE: Execute this script for training population == test populations initially to compute p0 and p1

### INPUTS:
fname_train_gwas=$1 ### filename of the training population's output of genome-wide association study using individual sequencing data (Indi-GWAS)
                    ### - ${TRAIN_prefix}_EMMAX_${TRAIN_kinship_id}.gwas ### tab-delimited; headerless; |col1:SNP|col2:EFF|col3:SE|col4:PVAL|
                    ### - ${TRAIN_prefix}_GEMMA_${TRAIN_kinship_id}.gwas ### tab-delimited; HEADER: |col01:chr|col02:rs|col03:n_miss|col04:allele1|col05:allele0|col06:af|col07:beta|col08:se|col09:logl_H1|col10:l_reml|col11:l_mle|col12:p_wald|col13:p_lrt|col14:p_score|
                    ### - ${TRAIN_prefix}_GCTA_${TRAIN_kinship_id}.gwas ### for GCTA's LR without GCTA's GRM binary input ### tab-delimited; HEADER: |col01:CHR|col02:SNP|col03:POS|col04:A1|col05:A2|col06:N|col07:AF1|col08:BETA|col09:SE|col10:P|
TEST_prefix=$2      ### prefix of the test or validation population
                    ### NOTE: We will be using the FULLSNPSET_${TEST_prefix}.bed, FULLSNPSET_${TEST_prefix}.bim, and FULLSNPSET_${TEST_prefix}.fam dataset
PLINK_MEM=$3        ### amount of memory to allocate to each instance of plink execution
SRC_DIR=$4          ### path to the softwares directory containing plink, emmax, gcta, and genomic_prediction.git
# fname_train_gwas=POP_01_EMMAX_GRM.gwas
# fname_train_gwas=POP_01_GEMMA_GRM.gwas
# fname_train_gwas=MULTIPOP_16_5665778165567716942_EMMAX_GRM.gwas
# TEST_prefix=POP_15
# TEST_prefix=MULTIPOP_16_5665778165567716942
# SRC_DIR=/data/Lolium/Softwares

### OUTPUTS:
### (1) ${TRAIN_prefix}_${TRAIN_model}_${TRAIN_kinship_id}-${TEST_prefix}.gp (tab-delimited; headerless: |col1:FAMID|col2:INDID|col3:PHENO|col4:NREF|col5:NALT|col6:SCORE|)
### where:
SAMPLING_ID=$(cut -d_ -f1 <<<$(echo $fname_train_gwas))
if [ $SAMPLING_ID == "POP" ]
then
  TRAIN_prefix=$(cut -d_ -f1-2 <<<$(echo $fname_train_gwas))
  TRAIN_model=$(cut -d_ -f3 <<<$(echo $fname_train_gwas))
  TRAIN_kinship_id=$(cut -d_ -f4 <<<$(echo $fname_train_gwas))
  TRAIN_kinship_id=${TRAIN_kinship_id%.*}
elif [ $SAMPLING_ID == "MULTIPOP" ]
then
  TRAIN_prefix=$(cut -d_ -f1-3 <<<$(echo $fname_train_gwas))
  TRAIN_model=$(cut -d_ -f4 <<<$(echo $fname_train_gwas))
  TRAIN_kinship_id=$(cut -d_ -f5 <<<$(echo $fname_train_gwas))
  TRAIN_kinship_id=${TRAIN_kinship_id%.*}
else
  echo "INVALID SAMPLING ID: " ${SAMPLING_ID}
fi
echo $TRAIN_prefix
echo $TRAIN_model
echo $TRAIN_kinship_id
echo $TEST_prefix
UNIQUE_ID=${TRAIN_prefix}_${TRAIN_model}_${TRAIN_kinship_id}-${TEST_prefix}

# ### GENVAL-->PLINK-SCORE
  # ### y --> BLUP(y)
  # ${SRC_DIR}/gcta_1.93.0beta/gcta64  --reml --grm ${TRAIN_prefix}_${TRAIN_kinship_id} --pheno ${TRAIN_prefix}.pheno --reml-pred-rand --out ${TRAIN_prefix}
  # ### output: ${TRAIN_prefix}.indi.blp (tab-delimited; headerless: |col1:FAMID|col2:INDID|col3:INTERM_VAR1|col4:GENVAL|col5:INTERM_VAR2|col6:RESID|)
  # ### BLUP(y) = GRM + X_train * b + e
  # ${SRC_DIR}/gcta_1.93.0beta/gcta64  --bfile ${TRAIN_prefix} --blup-snp ${TRAIN_prefix}.indi.blp --out ${TRAIN_prefix}
  # ### output: ${TRAIN_prefix}.snp.blp (tab-delimited; headerless: |col1:SNP_ID|col2:REF|col3:EFF|col4:RESID|)
  # # y_pred = X_test * b_bhat (1: SNP_ID at col1, 2: ALLELE_ID at col2, and SNP_SCORES at col3)
  # ${SRC_DIR}/plink --bfile FULLSNPSET_${TEST_prefix} --score ${TRAIN_prefix}.snp.blp 1 2 3 --memory ${PLINK_MEM} --out ${TRAIN_prefix}-${TEST_prefix}
  # mv ${TRAIN_prefix}-${TEST_prefix}.profile ${TRAIN_prefix}-${TEST_prefix}.gp
  # rm ${TRAIN_prefix}-${TEST_prefix}.log ${TRAIN_prefix}-${TEST_prefix}.nopred
  # ### output: ${TRAIN_prefix}-${TEST_prefix}.gp (tab-delimited; headerless: |col1:FAMID|col2:INDID|col3:PHENO|col4:NREF|col5:NALT|col6:SCORE|)

### GWAS-->PLINK-SCORE
### polygenic_score = X_test * b_bhat (1: SNP_ID at col1, 2: ALLELE_ID at col2, and SNP_SCORES at col3)
if [ ${TRAIN_model} == "EMMAX" ]
then
  # add the allele1 column (effects in emmax are all for allele1)
  cp ${TRAIN_prefix}_${TRAIN_model}_${TRAIN_kinship_id}.gwas ${UNIQUE_ID}.gwas.temp
  printf '2\n%.0s' $(seq 1 $(cut -f1 ${UNIQUE_ID}.gwas.temp | wc -l)) > ${UNIQUE_ID}.allele1.temp
  paste ${UNIQUE_ID}.gwas.temp ${UNIQUE_ID}.allele1.temp > ${UNIQUE_ID}.gwas
  rm ${UNIQUE_ID}.gwas.temp ${UNIQUE_ID}.allele1.temp
  echo "1 5 2" > ${UNIQUE_ID}.gwas.colIDX
elif [ ${TRAIN_model} == "GEMMA" ]
then
  echo "2 5 8" > ${UNIQUE_ID}.gwas.colIDX
  cp ${TRAIN_prefix}_${TRAIN_model}_${TRAIN_kinship_id}.gwas ${UNIQUE_ID}.gwas
elif [ ${TRAIN_model} == "GCTA" ]
then
  echo "2 4 8" > ${UNIQUE_ID}.gwas.colIDX
  cp ${TRAIN_prefix}_${TRAIN_model}_${TRAIN_kinship_id}.gwas ${UNIQUE_ID}.gwas
else
  echo "Invalid model:" ${model}
  exit
fi
${SRC_DIR}/plink \
              --bfile FULLSNPSET_${TEST_prefix} \
              --memory ${PLINK_MEM} \
              --score ${UNIQUE_ID}.gwas $(cat ${UNIQUE_ID}.gwas.colIDX) sum \
              --out ${UNIQUE_ID}
mv ${UNIQUE_ID}.profile ${UNIQUE_ID}.gp
rm ${UNIQUE_ID}.gwas
rm ${UNIQUE_ID}.log
rm ${UNIQUE_ID}.gwas.colIDX

# ### PLINK-SCORE-->PREDICTED-PHENOTYPE
# ### y_predicted = p0 + p1 * polygenic_score
# Rscript ${SRC_DIR}/genomic_prediction/src/GPASim_5.6_RMSE.r ${UNIQUE_ID}.gp
# ### Input:
# ### (1) ${UNIQUE_ID}.gp - filename of polygenic risk score summing-up output from plink --score
# ### Outputs:
# ### (1) ${UNIQUE_ID}.p0p1 - filename of the predictors of the y_pred~polygenic linear model; header: |col1:p0|col2:p1|col3:r2adj|
# ###     - if and only if: ${TRAIN_prefix}==${TEST_prefix} that is the second step of the 2-step prediction model building: y_pred = p0 + p1*plink_polygenic_score; where plink_polygenic_score was derived using GWAS output
# ### (2) ${UNIQUE_ID}.rmse - filename of the RSME output; header: |col01:train_pop|col02:test_pop|col03:train_model|col04:n_train|col05:n_test|col06:p0|col07:p1|col08:r2adj|col09:RMSE|col10:CORR|
