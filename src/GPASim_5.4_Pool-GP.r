#############################################
###         GENOMIC PREDICTION (GP)       ###
### using Pool sequencing (Pool-seq) data ###
#############################################
### 2-step genomic prediction using "GWAlpha_SNPwise" output
### Step 1: Calculate the polygenic "risk" score using test population's genotype data (.bed,.bim,&.fam) and training populations's SNP effects (.gwas)
### Step 2: Predict the test population's phenotypes using the training population's linear model connecting the actual phenotypes with plink-derived polygenic "risk" scores, i.e.:
###             - y_train ~ po + p1*polygenic_score_train
###             - y_test_predicted = p0 + p1*polygenic_score_test
### Straight-forward predition for the non-itertive models (i.e. ["ML", "REML"]_["LS", "GLMNET"])
#############################################
### Inputs:
### (1) filename of the training population's output of genome-wide association study using pool sequencing data (Pool-GWAS)
###     - string(prefix, "_GWAlpha_SNPwise.gwas")
###     - string(prefix, "_GWAlpha_REML_LS_", fst_id, ".gwas")
###     - string(prefix, "_GWAlpha_REML_LS_", fst_id, ".ranef")
###     - string(prefix, "_GWAlpha_REML_RR_", fst_id, ".gwas")
###     - string(prefix, "_GWAlpha_REML_RR_", fst_id, ".ranef")
###     - string(prefix, "_GWAlpha_REML_GELMNET_", fst_id, ".gwas")
###     - string(prefix, "_GWAlpha_REML_GELMNET_", fst_id, ".ranef")
###     - string(prefix, "_GWAlpha_REML_LASSO_", fst_id, ".gwas")
###     - string(prefix, "_GWAlpha_REML_LASSO_", fst_id, ".ranef")
###     - *.gwas: HEADER: |col1:CHROM|col2:POS|col3:ALLELES|col4:FREQ|col5:ALPHA|col6:PVALUES|col7:LOD|
###     - *.ranef: HEADERLESS: |col1:random effects|
### (2) filename of the test population's full SNP data in sync format (FULLSNPSET_*.sync)
### Output:
### (1) Genomic prection cross-validation statistics (tab_delimited; |col01:TRAIN_name|col02:TRAIN_model|col03:TRAIN_rancovar|col04:TEST_name|col05:n_train|col06:n_test|col07:p0|col08:p1|col09:r2adj|co10l:y_true|col11:y_pred|)
###     - For GWAlpha_SNPwise: string(TRAIN_name, "_", TRAIN_model, "-", TEST_name, ".gp")
###     - For GWAlpha_REML_[LS, RR, GLMNET, LASSO]_[HIVERT, WEIRCOCK]: string(TRAIN_name, "_", TRAIN_model, "_", TRAIN_rancovar, "-", TEST_name, ".gp")

args = commandArgs(trailingOnly=TRUE)
# args = c("POP_01_GWAlpha_SNPwise.gwas", "FULLSNPSET_POP_01_ALLELEFREQ.csv")
# args = c("POP_01_GWAlpha_RIDGE.gwas", "FULLSNPSET_POP_01_ALLELEFREQ.csv")
# args = c("POP_01_GWAlpha_MIXED_REML_HIVERT.gwas", "FULLSNPSET_POP_01_ALLELEFREQ.csv")
# args = c("MULTIPOP_8_9961816979469537538_GWAlpha_SNPwise.gwas", "FULLSNPSET_POP_01_ALLELEFREQ.csv")
# args = c("MULTIPOP_8_9961816979469537538_GWAlpha_RIDGE.gwas", "FULLSNPSET_POP_01_ALLELEFREQ.csv")
# args = c("MULTIPOP_8_9961816979469537538_GWAlpha_MIXED_REML_HIVERT.gwas", "FULLSNPSET_POP_01_ALLELEFREQ.csv")
# args = c("POP_01_GWAlpha_SNPwise.gwas", "FULLSNPSET_MULTIPOP_8_9961816979469537538_ALLELEFREQ.csv")
# args = c("POP_01_GWAlpha_RIDGE.gwas", "FULLSNPSET_MULTIPOP_8_9961816979469537538_ALLELEFREQ.csv")
# args = c("POP_01_GWAlpha_MIXED_REML_HIVERT.gwas", "FULLSNPSET_MULTIPOP_8_9961816979469537538_ALLELEFREQ.csv")
TRAIN_fname_gwas = args[1]
TEST_fname_geno = args[2]
### population names
if (length(grep("MULTIPOP_", TRAIN_fname_gwas))==0){
    TRAIN_name = paste(unlist(strsplit(TRAIN_fname_gwas, "_"))[1:2], collapse="_")
} else {
    TRAIN_name = paste(unlist(strsplit(TRAIN_fname_gwas, "_"))[1:3], collapse="_")
}
if (length(grep("MULTIPOP_", TEST_fname_geno))==0){
    TEST_name = paste(unlist(strsplit(TEST_fname_geno, "_"))[2:3], collapse="_")
} else {
    TEST_name = paste(unlist(strsplit(TEST_fname_geno, "_"))[2:4], collapse="_")
}
### model names with or without covar
TRAIN_model_split = unlist(strsplit(unlist(strsplit(unlist(strsplit(TRAIN_fname_gwas, TRAIN_name))[2], ".gwas"))[1], "_"))
MODEL_COVAR_name = paste(TRAIN_model_split[2:length(TRAIN_model_split)], collapse="_")
### training and test data filenames
TRAIN_fname_geno = paste0("FULLSNPSET_", TRAIN_name, "_ALLELEFREQ.csv")
TRAIN_fname_pheno = paste0(TRAIN_name , ".csv")
TEST_fname_pheno = paste0(TEST_name , ".csv")
### merge gwas-inferred allele effects and genotype data
GWAS = read.csv(TRAIN_fname_gwas)
GENO = read.csv(TEST_fname_geno, header=FALSE)
n_pools = ncol(GENO) - 3
colnames(GENO) = c("CHROM", "POS", "ALLELE", paste0("POOL_", 1:n_pools))
MERGED = merge(GWAS, GENO, by=c("CHROM", "POS", "ALLELE"))
### number of observations in the training and test sets
n_train = nrow(read.csv(TRAIN_fname_pheno, header=FALSE))
n_test = n_pools
### observed test population phenotypes
y_true = read.csv(TEST_fname_pheno, header=FALSE)[,2]
### SNPwise vs non-iterative models
if (length(grep("SNPwise", TRAIN_fname_gwas))!=0){
  ### step-2 training
  TRAIN_GENO = read.csv(TRAIN_fname_geno, header=FALSE)
  n_pools = ncol(TRAIN_GENO) - 3
  colnames(TRAIN_GENO) = c("CHROM", "POS", "ALLELE", paste0("POOL_", 1:n_pools))
  TRAIN_MERGED = merge(GWAS, TRAIN_GENO, by=c("CHROM", "POS", "ALLELE"))
  TRAIN_polygenic_score = t(as.matrix(TRAIN_MERGED[,8:ncol(TRAIN_MERGED)])) %*% TRAIN_MERGED[,5]
  TRAIN_y = read.csv(TRAIN_fname_pheno, header=FALSE)[,2]
  mod = lm(TRAIN_y ~ TRAIN_polygenic_score)
  p0 = mod$coefficients[1]
  p1 = mod$coefficients[2]
  r2adj = summary(mod)$adj.r.sq
  polygenic_score = t(as.matrix(MERGED[,8:ncol(MERGED)])) %*% MERGED[,5]
  y_pred = p0 + (p1 * polygenic_score)
} else {
  ### for non-iterative models
  p0 = NA
  p1 = NA
  r2adj = NA
  intercept = GWAS$ALPHA[GWAS$CHROM=="Intercept"]
  y_pred = intercept + t(as.matrix(MERGED[,8:ncol(MERGED)])) %*% MERGED[,5]
}
if (length(unlist(strsplit(MODEL_COVAR_name, "_"))) > 2){
  TRAIN_model = paste(unlist(strsplit(MODEL_COVAR_name, "_"))[1:3], collapse="_")
  TRAIN_rancovar = unlist(strsplit(MODEL_COVAR_name, "_"))[4]
} else {
  TRAIN_model = MODEL_COVAR_name
  TRAIN_rancovar = NA
}
### ouput
OUT = data.frame(TRAIN_name = rep(TRAIN_name, length(y_pred)),
                TRAIN_model = rep(TRAIN_model, length(y_pred)),
                TRAIN_rancovar = rep(TRAIN_rancovar, length(y_pred)),
                TEST_name = rep(TEST_name, length(y_pred)),
                n_train = rep(n_train, length(y_pred)),
                n_test = rep(n_test, length(y_pred)),
                p0 = rep(p0, length(y_pred)),
                p1 = rep(p1, length(y_pred)),
                r2adj = rep(r2adj, length(y_pred)),
                y_true = y_true,
                y_pred = y_pred
                )
fname_output = paste0(TRAIN_name, "_", MODEL_COVAR_name, "-", TEST_name, ".gp")
write.table(OUT, file=fname_output, col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)
