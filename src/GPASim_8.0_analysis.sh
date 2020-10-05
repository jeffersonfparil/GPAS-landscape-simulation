#!/bin/bash

### Input parameter description
# ( 1) ${SRC_DIR}           directory where qunatinemo2 folder, plink, emmax, gemma, gcta, and genomic_prediction folder are located
# ( 2) ${INPUT_DIR}         location of the simulation output

SRC_DIR=/data/Lolium/Softwares/genomic_prediction/src/
PLINK_DIR=/data/Lolium/Softwares/plink-1.07-x86_64/
INPUT_DIR=/data/Lolium/Quantitative_Genetics/LOLSIM_2020

# ### (0) Fix test_pop train_pop mix-up during AUC and RMSE datasets merging
# echo '
#     ### Input:
#     ### (1) GPAS_OUTPUT.csv
#     ### Output:
#     ### (1) GPAS_OUTPUT.csv (remedied)
#     args = commandArgs(trailingOnly=TRUE)
#     # args = c("GPASim_1REP_10NQTL_100NBGS_0.0001MIG_0.50SEL_0.10BGSEL_1GRAD/OUTPUT/GPAS_OUTPUT.csv")
#     fname_input = args[1]
#     ### load data
#     dat = read.csv(fname_input)
#     auto_dat = droplevels(dat[dat$TRAIN_POP == dat$TEST_POP, ])
#     for (i in 1:nrow(auto_dat)){
#         # i = 1
#         idx =   (dat$TRAIN_POP==auto_dat$TRAIN_POP[i]) & 
#                 (dat$SAMPLING_SCHEME==auto_dat$SAMPLING_SCHEME[i]) &
#                 (dat$GENOTYPING_SCHEME==auto_dat$GENOTYPING_SCHEME[i]) &
#                 (dat$MODEL==auto_dat$MODEL[i]) & 
#                 (dat$COVARIATE==auto_dat$COVARIATE[i])
#         dat$AUC[idx]                                = auto_dat$AUC[i]
#         dat$FRACTION_QTL_FIXED[idx]                 = auto_dat$FRACTION_QTL_FIXED[i]
#         dat$AUC_CORRECTED[idx]                      = auto_dat$AUC_CORRECTED[i]
#         dat$BONFERRONI_5PERCENT_TRUE_POSITIVE[idx]  = auto_dat$BONFERRONI_5PERCENT_TRUE_POSITIVE[i]
#         dat$BONFERRONI_5PERCENT_FALSE_POSITIVE[idx] = auto_dat$BONFERRONI_5PERCENT_FALSE_POSITIVE[i]
#         dat$BONFERRONI_5PERCENT_QTL_ID[idx]         = auto_dat$BONFERRONI_5PERCENT_QTL_ID[i]
#     }
#     write.table(dat, file=fname_input, sep=",", quote=FALSE, row.names=FALSE, col.names=TRUE)
# ' > remedy_ouput.R
# ### bacxk-up first
# for f in $(find ${INPUT_DIR}/GPASim_*/OUTPUT/GPAS_OUTPUT.csv)
# do
#     echo $f
#     cp $f ${f}.bk
# done
# ### remedy
# time \
# parallel \
# Rscript remedy_ouput.R {} ::: $(find ${INPUT_DIR}/GPASim_*/OUTPUT/GPAS_OUTPUT.csv)
# rm remedy_ouput.R

# ### (1) Merge complete simulation output and parse into an RDS binary file while also restricting maximum RMSE to 1 for computational ease (since an RMSE of 1 simply says that the genomic prediction error maximized, i.e. predicting 1 when the truth is zero and vice versa)
# ### Outputs: "MERGED_GPAS_OUTPUT.csv" and "MERGED_GPAS_OUTPUT.rds"
# ### streamline the dataset including intra and inter-population data
# echo 'args = commandArgs(trailingOnly=TRUE) ### args = c("GPAS_OUTPUT.csv")
#       dat = read.csv(args[1])
#       ### restrict maximum RMSE to 1
#       dat$RMSE[dat$RMSE > 1] = 1
#       ### remove the fileds we will not be using
#       dat = droplevels(dat[, c(1:13, 22:29, 33, 39)])
#       gc()
#       ### remove auto cross-validataions
#       dat = droplevels(dat[dat$TRAIN_POP != dat$TEST_POP, ])
#       gc()
#       ### remove at least one common population in the training and validation sets
#       dat$TRAIN_POP = as.character(dat$TRAIN_POP)
#       dat$TEST_POP = as.character(dat$TEST_POP)
#       idx = c()
#       pb = txtProgressBar(min=0, max=nrow(dat), initial=0, style=3)
#       for (i in 1:nrow(dat)){
#           # i = 1
#           # print(i)
#           setTxtProgressBar(pb, i)
#           train_pops = unlist(strsplit(dat$TRAIN_POP[i], ":"))
#           test_pops = unlist(strsplit(dat$TEST_POP[i], ":"))
#           if ( ((sum(train_pops %in% test_pops)==0) & (sum(test_pops %in% train_pops)==0)) == TRUE ){
#               idx = c(idx, i)
#           }
#       }
#       close(pb)
#       dat = droplevels(dat[idx, ])
#       ### save csv
#       write.table(dat, file=paste0(dirname(args[1]), "/STREAMLINED_OUTPUT.csv"), sep=",", quote=FALSE, row.names=FALSE, col.names=TRUE)
# ' > streamline_parallel.r
# time \
# parallel \
# Rscript streamline_parallel.r {} ::: $(find ${INPUT_DIR}/GPASim_*/OUTPUT/GPAS_OUTPUT.csv)
# ### prepare the header: (1)REPLICATE,(2)EFFECTIVE_POP_SIZE,(3)NPOP,(4)NGEN,(5)NLOCI,(6)NALLELES,(7)NQTL,(8)NBGS,(9)QTL_EFF_DIST,(10)SELECTION_INTENSITY,(11)BGS_INTENSITY,(12)MIGRATION_RATE,(13)QTL_GRADIENT,(14)PHENO_MEAN_BETADIST_SHAPE1,(15)PHENO_MEAN_BETADIST_SHAPE2,(16)PHENO_VAR_BETADIST_SHAPE1,(17)PHENO_VAR_BETADIST_SHAPE2,(18)GENO_FST_BETADIST_SHAPE1,(19)GENO_FST_BETADIST_SHAPE2,(20)NPOOLS_PER_POP,(21)NLIB_MAXIMUM,(22)SAMPLING_SCHEME,(23)GENOTYPING_SCHEME,(24)TRAIN_POP,(25)TEST_POP,(26)NTRAIN,(27)NTEST,(28)MODEL,(29)COVARIATE,(30)p0,(31)p1,(32)r2adj,(33)AUC,(34)FRACTION_QTL_FIXED,(35)AUC_CORRECTED,(36)BONFERRONI_5PERCENT_TRUE_POSITIVE,(37)BONFERRONI_5PERCENT_FALSE_POSITIVE,(38)BONFERRONI_5PERCENT_QTL_ID,(39)RMSE,(40)CORRELATION
# head -n1 <$(find ${INPUT_DIR}/GPASim_*/OUTPUT/STREAMLINED_OUTPUT.csv | head -n1) > MERGED_GPAS_OUTPUT.csv
# ### concatenate while excluding the header from each file
# for f in $(find ${INPUT_DIR}/GPASim_*/OUTPUT/STREAMLINED_OUTPUT.csv)
# do
#     tail -n+2 $f >> MERGED_GPAS_OUTPUT.csv
# done
# ### parse into an R-readable binary file
# echo 'dat = read.csv("MERGED_GPAS_OUTPUT.csv")
#       saveRDS(dat, file="MERGED_GPAS_OUTPUT.rds")
# ' > parse_into_rds.r
# time \
# Rscript parse_into_rds.r
# rm parse_into_rds.r

### (2) Extract GPAS performance as a function of the number of merged independent GPAS experiments per population
### Output: "GPAS_OUTPUT_NPOP_PERFORMANCE.csv"
time \
parallel \
Rscript ${SRC_DIR}/GPASim_8.1_performance_npop.r {} ::: $(find ${INPUT_DIR}/GPASim_*/OUTPUT/GPAS_OUTPUT.csv)

### (3) Merge the extracted GPAS performance as function of the number of population samples (1 to 100)
### Outputs: "GPAS_OUTPUT_NPOP_PERFORMANCE.csv" and "GPAS_OUTPUT_NPOP_PERFORMANCE.rds"
### prepare the header: (1)REPLICATE,(2)EFFECTIVE_POP_SIZE,(3)NPOP,(4)NGEN,(5)NLOCI,(6)NALLELES,(7)NQTL,(8)NBGS,(9)QTL_EFF_DIST,(10)SELECTION_INTENSITY,(11)BGS_INTENSITY,(12)MIGRATION_RATE,(13)QTL_GRADIENT,(14)PHENO_MEAN_BETADIST_SHAPE1,(15)PHENO_MEAN_BETADIST_SHAPE2,(16)PHENO_VAR_BETADIST_SHAPE1,(17)PHENO_VAR_BETADIST_SHAPE2,(18)GENO_FST_BETADIST_SHAPE1,(19)GENO_FST_BETADIST_SHAPE2,(20)
head -n1 <$(find ${INPUT_DIR}/GPASim_*/OUTPUT/GPAS_OUTPUT_NPOP_PERFORMANCE.csv | head -n1) > MERGED_GPAS_OUTPUT_NPOP_PERFORMANCE.csv
### concatenate while excluding the header from each file
for f in $(find ${INPUT_DIR}/GPASim_*/OUTPUT/GPAS_OUTPUT_NPOP_PERFORMANCE.csv)
do
    tail -n+2 $f >> MERGED_GPAS_OUTPUT_NPOP_PERFORMANCE.csv
done
### parse into an R-readable binary file
echo 'dat = read.csv("MERGED_GPAS_OUTPUT_NPOP_PERFORMANCE.csv")
      saveRDS(dat, file="MERGED_GPAS_OUTPUT_NPOP_PERFORMANCE.rds")
' > parse_into_rds.r
time \
Rscript parse_into_rds.r
rm parse_into_rds.r

### (4) Extract genetic variability summary statistics
### Output: "MERGED_VARIATION_SUMMSTATS.csv"
###     i.e. Mean phenotypes, phenotype variance, number of QTL captured, expected heterozygosity, and observed heterozygosity
### Compute the heterozogosities
echo '#!/bin/bash
      FAM=$1
      PLINK_DIR=$2
      # FAM=/data/Lolium/Quantitative_Genetics/LOLSIM_2020/GPASim_5REP_50NQTL_100NBGS_0.01MIG_0.50SEL_0.10BGSEL_2GRAD/RESOURCES/POP_001.fam
      # PLINK_DIR=/data/Lolium/Softwares/plink-1.07-x86_64
      DIR=$(dirname $FAM)
      BASENAME=$(basename $FAM)
      ROOTNAME=${BASENAME%.fam*}
      ### Calculate heterozygosities
      ${PLINK_DIR}/plink --bfile ${DIR}/${ROOTNAME} --hardy --noweb --out ${DIR}/${ROOTNAME}
      sed "s/ \+/,/g" ${DIR}/${ROOTNAME}.hwe | sed "s/^,//" | sed "s/,$//g" > ${DIR}/${ROOTNAME}.hwe.csv
' > compute_heterozigosities.sh
chmod +x compute_heterozigosities.sh
time \
for dir in $(ls ${INPUT_DIR} | grep GPASim)
do
    parallel ./compute_heterozigosities.sh {} ${PLINK_DIR} ::: $(ls ${dir}/RESOURCES/POP_*.fam)
done
rm compute_heterozigosities.sh
### Compute QTL allele frequencies
echo '#!/bin/bash
      FAM=$1
      PLINK_DIR=$2
      # FAM=/data/Lolium/Quantitative_Genetics/LOLSIM_2020/GPASim_5REP_50NQTL_100NBGS_0.01MIG_0.50SEL_0.10BGSEL_2GRAD/RESOURCES/OUT.fam
      # PLINK_DIR=/data/Lolium/Softwares/plink-1.07-x86_64
      DIR=$(dirname $FAM)
      BASENAME=$(basename $FAM)
      ROOTNAME=${BASENAME%.fam*}
      ### Extract population clustering ID
      cut -d" " -f1,2 ${FAM} > ${DIR}/within_id.col12
      cut -d" " -f1   ${FAM} > ${DIR}/within_id.col3
      paste -d" " ${DIR}/within_id.col12 ${DIR}/within_id.col3 > ${DIR}/within_id.txt
      ### Calculate allele frequencies of the QTL
      ${PLINK_DIR}/plink --bfile ${DIR}/${ROOTNAME} --freq --within ${DIR}/within_id.txt --noweb --out ${DIR}/${ROOTNAME}
      sed "s/ \+/,/g" ${DIR}/${ROOTNAME}.frq.strat | sed "s/^,//" | sed "s/,$//g" > ${DIR}/${ROOTNAME}.frq.csv
      rm ${DIR}/within_id.col12 ${DIR}/within_id.col3 ${DIR}/within_id.txt
' > compute_QTL_freq.sh
chmod +x compute_QTL_freq.sh
time \
parallel ./compute_QTL_freq.sh {} ${PLINK_DIR} ::: $(ls ${INPUT_DIR}/*/RESOURCES/OUT.fam)
rm compute_QTL_freq.sh

### compute the phenotype means, variances, number of QTL captured, mean expected and observed heterozygosities per population, and QTL allele frequencies
time \
parallel \
Rscript ${SRC_DIR}/GPASim_8.2_variabilities_per_pop.r {} ::: $(ls ${INPUT_DIR}/ | grep GPASim_)
### merge across landscapes
head -n1 <$(find ${INPUT_DIR}/GPASim_*/OUTPUT/VARIATION_PER_POP.csv | head -n1) > MERGED_VARIATION_SUMMSTATS.csv
### concatenate while excluding the header from each file
for f in $(find ${INPUT_DIR}/GPASim_*/OUTPUT/VARIATION_PER_POP.csv)
do
    tail -n+2 $f >> MERGED_VARIATION_SUMMSTATS.csv
done


### (5) Analysis in R to answer the following questions:
###### (1) How many populations to sample?
###### (2) When to use Indi-seq or Pool-seq?
###### (3) Which populations to select?
### Outputs: main figures and tables (0.*_*.svg to 3.*_*.svg, and 0.*_*.csv to 3.*_*.csv)
time \
nohup \
Rscript ${SRC_DIR}/GPASim_8.3_analysis.r \
    MERGED_GPAS_OUTPUT_NPOP_PERFORMANCE.rds \
    MERGED_VARIATION_SUMMSTATS.csv &






###@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@###
### Miscellaneous show that for one-SNP-at-a-time polygenic score-based prediction is linear:
### y_train = s0 + s1*polygenic_score_train
###@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@###
################
### INDI-SEQ ###
################
DIR=/data/Lolium/Quantitative_Genetics/LOLSIM_2020/GPASim_1REP_50NQTL_100NBGS_0.001MIG_0.50SEL_0.10BGSEL_0GRAD/RESOURCES
cd ${DIR}
time \
for prefix in $(ls POP_*.fam | cut -d'.' -f1)
do
    # prefix=POP_001
    echo $prefix
    ### estimate STD
    /data/Lolium/Softwares/gemma-0.98.1-linux-static -bfile ${prefix} -gk 2 -o ${prefix}_STANDARDIZED
    mv output/${prefix}_STANDARDIZED.sXX.txt ${prefix}_STANDARDIZED.kin
    rm -R output/
    ### GEMMA
    kinship=${prefix}_STANDARDIZED.kin
    xiferp=$(echo $prefix | rev); pihsnik=$(echo ${kinship%.*} | rev)
    kinship_id=$(echo ${pihsnik%_$(echo $xiferp)*} | rev)
    /data/Lolium/Softwares/gemma-0.98.1-linux-static  -bfile ${prefix} \
                                        -k ${kinship} \
                                        -lmm 4 \
                                        -o ${prefix}_GEMMA_${kinship_id}
    mv output/${prefix}_GEMMA_${kinship_id}.assoc.txt ${prefix}_GEMMA_${kinship_id}.gwas
    rm -R output/
    ### predice polygenic score
    /data/Lolium/Softwares/plink \
                --bfile ${prefix} \
                --memory 200000 \
                --score ${prefix}_GEMMA_${kinship_id}.gwas 2 5 8 sum \
                --out ${prefix}_GEMMA_${kinship_id}
    mv ${prefix}_GEMMA_${kinship_id}.profile ${prefix}_GEMMA_${kinship_id}.gp
    ### plot observed phenotypic values vs polygenic scores
    echo '
        args = commandArgs(trailingOnly=TRUE)
        fname_scores = args[1]
        # fname_scores = "POP_001_GEMMA_STANDARDIZED.gp"
        dat = read.table(fname_scores, header=TRUE)
        mod = lm(PHENO ~ SCORESUM, data=dat)
        p0 = coef(mod)[1]
        p1 = coef(mod)[2]
        r2adj = summary(mod)$adj.r.squared
        n_train = nrow(dat)
        svg(sub(".gp", "_OBSxPOLY.svg", fname_scores), width=7, height=7)
        plot(x=dat$PHENO, y=dat$SCORESUM, xlab="Observed Phenotypic Value", ylab="Polygenic Score", pch=20, col=rgb(0.5,0.5,0.5, 0.5))
        score = seq(from=min(dat$SCORESUM), to=max(dat$SCORESUM), length=100)
        pheno = predict(mod, newdata=data.frame(SCORESUM=score))
        lines(x=pheno, y=score, lty=2, lwd=2)
        grid()
        legend("topleft", legend=paste0("n=", n_train))
        dev.off()
        write.table(r2adj, file=sub(".gp", "_OBSxPOLY.r2adj", fname_scores), col.names=FALSE, row.names=FALSE, sep=",")
    ' > plot_obsPhen_vs_polyScore.r
    Rscript plot_obsPhen_vs_polyScore.r ${prefix}_GEMMA_${kinship_id}.gp
done
### summarise R2_adjusted
cat POP_*GEMMA*.r2adj > merged_GEMMA.r2adj
rm POP_*GEMMA*.r2adj
echo '
dat = read.csv("merged_GEMMA.r2adj", header=FALSE)
print(mean(dat$V1))
print(sd(dat$V1)*sqrt(nrow(dat)-1)/nrow(dat))
' > mean_se_r2adj.r
Rscript mean_se_r2adj.r
### clean-up
rm plot_obsPhen_vs_polyScore.r mean_se_r2adj.r


################
### POOL-SEQ ###
################
### script wrapper to convert vcf into sync with 5 pools per population
echo '#!/bin/bash
  i=$1
  NPOOLS=$2
  if [ $NPOOLS == 1 ]
  then
    ### add prefix to single pool per population
    cp ${i}.vcf ONEPOOL_${i}.vcf
    julia /data/Lolium/Softwares/genomic_prediction/src/GPASim_2.1_vcf2sync.jl ONEPOOL_${i}.vcf ${i}.fam ${NPOOLS}
    rm ONEPOOL_${i}.vcf
  else
    julia /data/Lolium/Softwares/genomic_prediction/src/GPASim_2.1_vcf2sync.jl ${i}.vcf ${i}.fam ${NPOOLS}
  fi
' > vcf2sync.sh; chmod +x vcf2sync.sh
### main script wrapper
echo '
    using GWAlpha
    using DelimitedFiles
    using RCall
    # ARGS = ["POP_001"]
    prefix = ARGS[1]
    filename_sync = string(prefix, ".sync")
    filename_phen_py = string(prefix, ".py")
    filename_phen_csv = string(prefix, ".csv")
    filename_gwas = string(prefix, "_GWAlpha_SNPwise.gwas")
    filename_geno_csv = string(prefix, "_ALLELEFREQ.csv")
' > GWAlpha_polygenic_score_calc.jl
echo "
    ### Calculate MAF and set DEPTH=1 for simplicity
    MAF = 1.00/(2*sum(DelimitedFiles.readdlm(filename_phen_csv, ',')[:,1]))
    DEPTH = 1
" >> GWAlpha_polygenic_score_calc.jl
echo '
    @time GWAlpha.PoolGPAS(filename_sync=filename_sync, filename_phen=filename_phen_py, model="GWAlpha", maf=MAF, depth=DEPTH, fpr=0.01)
    # output: comma-delimited; HEADER: |col1:CHROM|col2:POS|col3:ALLELES|col4:FREQ|col5:ALPHA|col6:PVALUES|col7:LOD|
    mv(string(prefix, "-GWAlpha-OUTPUT.csv"), filename_gwas, force=true)
    GWAlpha.sync_processing_module.sync_parse(filename_sync)
' >> GWAlpha_polygenic_score_calc.jl
echo "
    GENO = DelimitedFiles.readdlm(filename_geno_csv, ',')

    ### load allelic effects from the training population
    gwas = DelimitedFiles.readdlm(filename_gwas, ',', header=true)[1]
    idx = [ sum(sum(reshape(GENO[i, 1:3],1,3) .== gwas[:,1:3], dims=2) .== 3) > 0 for i in 1:size(GENO,1) ]

    ### extract the testing genotype and the training allelic effects
    X = convert(Array{Float64,2}, GENO[idx,4:end])'
    b_hat = convert(Array{Float64,1}, gwas[:,5])
    polygenic_score = X * b_hat
    y_test_true = DelimitedFiles.readdlm(filename_phen_csv, ',')[:,2]

    @rput polygenic_score y_test_true filename_gwas
" >> GWAlpha_polygenic_score_calc.jl
echo '
    R"
    PHENO = y_test_true
    SCORESUM = polygenic_score
    mod = lm(PHENO ~ SCORESUM)
    r2adj = summary(mod)$adj.r.squared
    n_train = length(PHENO)
    svg(sub(\".gwas\", \"_OBSxPOLY.svg\", filename_gwas), width=7, height=7)
    plot(x=PHENO, y=SCORESUM, xlab=\"Observed Phenotypic Value\", ylab=\"Polygenic Score\", pch=19, , col=rgb(0.5,0.5,0.5, 0.5))
    score = seq(from=min(SCORESUM), to=max(SCORESUM), length=100)
    pheno = predict(mod, newdata=data.frame(SCORESUM=score))
    lines(x=pheno, y=score, lty=2, lwd=2)
    grid()
    legend(\"topleft\", legend=paste0(\"n=\", n_train))
    dev.off()
    write.table(r2adj, file=sub(\".gwas\", \"_OBSxPOLY.r2adj\", filename_gwas), col.names=FALSE, row.names=FALSE, sep=\",\")
    "
' >> GWAlpha_polygenic_score_calc.jl
### parallelizable script
echo '#!/bin/bash
    prefix=$1
    echo $prefix
    /data/Lolium/Softwares/plink  --bed ${prefix}.bed \
                    --bim ${prefix}.bim \
                    --fam ${prefix}.fam \
                    --recode vcf \
                    --threads 1 \
                    --memory 8750 \
                    --out ${prefix}
    ### convert vcf into sync with 5 pools per population
    ./vcf2sync.sh ${prefix} 5
    ### perform GWAlpha, plot and compute R2_adjusted
    julia GWAlpha_polygenic_score_calc.jl ${prefix}
' > GWAlpha_execute_in_parallel.sh; chmod +x GWAlpha_execute_in_parallel.sh
### execute
time \
parallel ./GWAlpha_execute_in_parallel.sh {} ::: $(ls POP_*.fam | cut -d'.' -f1)
### summary R2_adjusted
cat POP_*GWAlpha*.r2adj > merged_GWAlpha.r2adj
rm POP_*GWAlpha*.r2adj
echo '
dat = read.csv("merged_GWAlpha.r2adj", header=FALSE)
print(mean(dat$V1))
print(sd(dat$V1)*sqrt(nrow(dat)-1)/nrow(dat))
' > mean_se_r2adj.r
Rscript mean_se_r2adj.r
### clean-up
rm vcf2sync.sh GWAlpha_polygenic_score_calc.jl GWAlpha_execute_in_parallel.sh mean_se_r2adj.r
