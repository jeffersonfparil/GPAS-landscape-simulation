#!/bin/bash
###################################################
### ESTIMATE KINSHIP MATRICES FOR INDI-SEQ DATA ###
################################################### NOTE: Also filters by MAF=1/pop_size

### INPUTS:
prefix=$1
# prefix=POP_01
SRC_DIR=$2 ### directory where plink and genomic_prediction/ are located
PLINK_MEM=$3 ### set the amount of memory to allocate for each instance of plink execution
NTHREADS=1 ### set to 1 because we will be parallelizing across Indi-seq data

### OUTPUTS:
### MAF-filtered bed, bim, & fam files (old files were renamed with the prefix: FULLSNPSET_*):
### (1) ${prefix}.bed
### (2) ${prefix}.bim
### (3) ${prefix}.fam
### Headerless square symmetric kinship matrices and GCTA GRM files:
### (4) ${prefix}_STANDARDIZED.kin (GCTA-specific GRM files: ${prefix}_STANDARDIZED.grm.bin,${prefix}_STANDARDIZED.grm.id, & ${prefix}_STANDARDIZED.grm.N.bin)
### (5) ${prefix}_GRM.kin (GCTA-specific GRM files: ${prefix}_GRM.grm.bin,${prefix}_GRM.grm.id, & ${prefix}_GRM.grm.N.bin)

###################################
### Filter by MAF of 1/pop_size ###
###################################
MAF=$(echo "scale=10; 1 / $(cat ${prefix}.fam | wc -l)" | bc)
${SRC_DIR}/plink  --bfile ${prefix} \
                  --maf ${MAF}\
                  --snps-only \
                  --threads ${NTHREADS} \
                  --memory ${PLINK_MEM} \
                  --make-bed \
                  --out ${prefix}_MAF_filtered
### rename original dataset (to be used in genomic prediction cross-validation)
mv ${prefix}.bed FULLSNPSET_${prefix}.bed
mv ${prefix}.bim FULLSNPSET_${prefix}.bim
mv ${prefix}.fam FULLSNPSET_${prefix}.fam
### rename filtered dataset
mv ${prefix}_MAF_filtered.bed ${prefix}.bed
mv ${prefix}_MAF_filtered.bim ${prefix}.bim
mv ${prefix}_MAF_filtered.fam ${prefix}.fam
rm ${prefix}_MAF_filtered.log

################################################
### gemma: standardized kinship matrix XX'/p ###
################################################
${SRC_DIR}/gemma-0.98.1-linux-static -bfile ${prefix} -gk 2 -o ${prefix}_STANDARDIZED
mv output/${prefix}_STANDARDIZED.sXX.txt ${prefix}_STANDARDIZED.kin
# txt STANDARDIZED to bin STANDARDIZED
Rscript txt2bin.r ${prefix}
# generate sparse GRM as gcta --fastGWA-lmm input (keep 95% of pairs)
${SRC_DIR}/gcta_1.93.0beta/gcta64 --grm ${prefix}_STANDARDIZED \
                                  --make-bK-sparse 0.05 \
                                  --thread-num ${NTHREADS} \
                                  --out ${prefix}_STANDARDIZED_sparse

##############################################################
### gcta: genetic relatedness matrix (see Yang et al 2010) ###
##############################################################
# dummy male to prevent: "Empty sample_subset is not currently permitted" error!
cp ${prefix}.fam ${prefix}.fam.bk
cut -d' ' -f1-4 ${prefix}.fam.bk > col1234.temp.${prefix}
cut -d' ' -f5 ${prefix}.fam.bk | tail -n+2 > col5.temp.${prefix}
echo -e "1" >> col5.temp.${prefix}
cut -d' ' -f6 ${prefix}.fam.bk > col6.temp.${prefix}
paste -d' ' col1234.temp.${prefix} \
            col5.temp.${prefix} \
            col6.temp.${prefix} > ${prefix}.fam
${SRC_DIR}/gcta_1.93.0beta/gcta64 --bfile ${prefix} \
                                  --make-grm \
                                  --autosome \
                                  --thread-num ${NTHREADS} \
                                  --out ${prefix}_GRM
# bin GRM to txt GRM
Rscript bin2txt.r ${prefix}
# generate sparse GRM as gcta --fastGWA-lmm input (keep 95% of pairs)
${SRC_DIR}/gcta_1.93.0beta/gcta64 --grm ${prefix}_GRM \
                                  --make-bK-sparse 0.05 \
                                  --thread-num ${NTHREADS} \
                                  --out ${prefix}_GRM_sparse
################
### clean-up ###
################
rm ${prefix}.fam.bk
rm ${prefix}_STANDARDIZED*.log ${prefix}_GRM*.log
rm *.temp.${prefix}
