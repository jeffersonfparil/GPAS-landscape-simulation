#!/bin/bash
###################################################
###   GENOME-WIDE ASSOCIATION STUDIES (GWAS)    ###
### using Individual sequencing (Indi-seq) data ###
###################################################
### INPUTS:
prefix=$1     ### prefix of the bed-bim-fam plink file of genotype and phenotype data (POP_* or MULTIPOP_*_*)
kinship=$2    ### kinship file to use (headerless, tab-delimited, square-symmetric matrix) (GRM or STANDARDIZED)
              ### NOTE: for GCTA the kinship matrix will AUTOMATICALLY be the GCTA-specific GRM sparse files instead of the headerless tab-delimited kinship file
NTHREADS=$3   ### number of compute cores or threads to use for parallel computation
PLINK_MEM=$4  ### amount of memory to allocate to each instance of plink execution
SRC_DIR=$5    ### path to the softwares directory containing plink, emmax, gcta, and genomic_prediction.git
# prefix=POP_01
# kinship=${prefix}_GRM.kin
# kinship=${prefix}_STANDARDIZED.kin
# NTHREADS=1
# SRC_DIR=/data/Lolium/Softwares
### OUTPUTS:
### (1) ${prefix}_EMMAX_${kinship_id}.gwas ### tab-delimited; headerless; |col1:SNP|col2:EFF|col3:SE|col4:PVAL|
### (2) ${prefix}_GEMMA_${kinship_id}.gwas ### tab-delimited; HEADER: |col01:chr|col02:rs|col03:n_miss|col04:allele1|col05:allele0|col06:af|col07:beta|col08:se|col09:logl_H1|col10:l_reml|col11:l_mle|col12:p_wald|col13:p_lrt|col14:p_score|
### (3) ${prefix}_GCTA_${kinship_id}.gwas ### for GCTA's LR without GCTA's GRM binary input ### tab-delimited; HEADER: |col01:CHR|col02:SNP|col03:POS|col04:A1|col05:A2|col06:N|col07:AF1|col08:BETA|col09:SE|col10:P|
### (4) ${prefix}_GCTA_GRM-MLM.gwas ### for GCTA's MLM with GCTA's GRM binary input
### where:
xiferp=$(echo $prefix | rev); pihsnik=$(echo ${kinship%.*} | rev)
kinship_id=$(echo ${pihsnik%_$(echo $xiferp)*} | rev)
echo "#############"
echo "### EMMAX ###"
echo "#############"
# transpose genotype data
${SRC_DIR}/plink  --bfile ${prefix} \
                  --recode12 \
                  --transpose \
                  --threads ${NTHREADS} \
                  --memory ${PLINK_MEM} \
                  --out ${prefix}_${kinship_id}_transposed
# extract phenotype data
cut -d' ' -f1,2,6 ${prefix}_${kinship_id}_transposed.tfam > ${prefix}_${kinship_id}_transposed.pheno
# run emmax association
${SRC_DIR}/emmax-intel64 -d ${NTHREADS} \
                         -t ${prefix}_${kinship_id}_transposed \
                         -p ${prefix}_${kinship_id}_transposed.pheno \
                         -k ${kinship} \
                         -o ${prefix}_EMMAX_${kinship_id}
# output: tab-delimited; headerless; |col1:SNP|col2:EFF|col3:SE|col4:PVAL|
mv ${prefix}_EMMAX_${kinship_id}.ps ${prefix}_EMMAX_${kinship_id}.gwas
# clean-up
rm ${prefix}_${kinship_id}_transposed* ${prefix}_EMMAX_${kinship_id}.log ${prefix}_EMMAX_${kinship_id}.reml
echo "#############"
echo "### GEMMA ###"
echo "#############"
# run gemma with -lmm=4 where Wald's test, likelihood ratio test, & score test are all performed
${SRC_DIR}/gemma-0.98.1-linux-static  -bfile ${prefix} \
                                      -k ${kinship} \
                                      -lmm 4 \
                                      -o ${prefix}_GEMMA_${kinship_id}
# output: tab-delimited; HEADER: |col01:chr|col02:rs|col03:n_miss|col04:allele1|col05:allele0|col06:af|col07:beta|col08:se|col09:logl_H1|col10:l_reml|col11:l_mle|col12:p_wald|col13:p_lrt|col14:p_score|
mv output/${prefix}_GEMMA_${kinship_id}.assoc.txt ${prefix}_GEMMA_${kinship_id}.gwas
# clean-up
rm output/${prefix}_GEMMA_${kinship_id}.*
echo "############"
echo "### GCTA ###"
echo "############" ### NOTE: Plink-derived Balding-Nichols kinship matrix keeps explaining 100% of the phenotypic variance! That's why it's excluded and we're only including emmax's standardized kinship and GCTA's GRM matrices
# extract phenotype data
cut -d' ' -f1,2,6 ${prefix}.fam > ${prefix}_${kinship_id}.pheno
# add identifiers to the kinship matrix
cut -d' ' -f1 ${prefix}.fam > ${prefix}_${kinship}_col1_temp
cut -d' ' -f2 ${prefix}.fam > ${prefix}_${kinship}_col2_temp
paste ${prefix}_${kinship}_col1_temp ${prefix}_${kinship}_col2_temp ${kinship} > ${kinship}_WITH_ID.temp
# run gcta using ${kinship} file as fixed covariates on top of the gcta-derived GRM (kinship matrix)
${SRC_DIR}/gcta_1.93.0beta/gcta64 --bfile ${prefix} \
                                  --grm-sparse ${prefix}_${kinship_id}_sparse \
                                  --fastGWA-mlm \
                                  --pheno ${prefix}_${kinship_id}.pheno \
                                  --threads ${NTHREADS} \
                                  --out ${prefix}_GCTA_${kinship_id} || \
echo "ERROR ${prefix}_GCTA_${kinship_id}! Sample size is probably too small - bump it up to at least 1000 individuals. Or genotype data too small resulting to not enough null SNPs."
# output: tab-delimited; HEADER: |col01:CHR|col02:SNP|col03:POS|col04:A1|col05:A2|col06:N|col07:AF1|col08:BETA|col09:SE|col10:P|
mv ${prefix}_GCTA_${kinship_id}.fastGWA ${prefix}_GCTA_${kinship_id}.gwas || \
echo "ERROR ${prefix}_GCTA_${kinship_id}! Sample size is probably too small - bump it up to at least 1000 individuals. Or genotype data too small resulting to not enough null SNPs."
# clean-up
rm ${prefix}_${kinship_id}.pheno ${kinship}_WITH_ID.temp ${prefix}_GCTA_${kinship_id}.log ${prefix}_${kinship}_col*_temp
