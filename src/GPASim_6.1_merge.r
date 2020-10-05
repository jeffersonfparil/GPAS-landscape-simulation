####################################################################################################
### MERGE GENOMIC PREDICTION AND GENOME-WIDE ASSOCIATION OUTPUT (GP:RMSE AND GWAS:AUC_CORRECTED) ###
###                   TOGETHER WITH LANDSCAPE SIMULATION INPUT PARAMETERS,                       ###
###   AS WELL AS PHENOTYPIC MEANS, AND VARIATION AD GENETIC DIFFERENTION SUMMARY STATISTICS       ###
###      ALSO STREAMLINE: REMOVE AUTO-CROSS-VALIDATION AND LESS 2-POOL TEST POPULATIONS          ###
####################################################################################################

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
###                                                           |col11:BONFERRONI_5PERCENT_TRUE_POSITIVE|col12:BONFERRONI_5PERCENT_FALSE_POSITIVE|col13:BONFERRONI_5PERCENT_QTL_ID|col14:RMSE|col15:CORRELATION|)
args = commandArgs(trailingOnly=TRUE)
# args = c("population_parameters.csv", "LANDSCAPE.stat", "LANDSCAPE.fst", "OUTPUT_MERGED_AUC.csv", "OUTPUT_MERGED_RMSE.csv")

### Load datasets
PARAMETERS = read.csv(args[1])
LANDSCAPE_PHENOTYPES = read.csv(args[2])
FST = read.csv(args[3])
AUC = read.csv(args[4])
RMSE = read.csv(args[5])

### rename POP* into TEST_POP* for merging with the cross-validation results of RMSE
colnames(AUC) = c("TRAIN_SAMPLE_ID", "SAMPLING", "TRAIN_POP", "MODEL", "COVARIATE", "AUC", "FRACTION_QTL_FIXED", "BONFERRONI_5PERCENT_TRUE_POSITIVE", "BONFERRONI_5PERCENT_FALSE_POSITIVE", "BONFERRONI_5PERCENT_QTL_ID")

### compute AUC corrected by the fraction of QTL that are discoverable in the experiment (i.e. fraction of the QTL which were not filtered out because they were fixed or lost)
AUC$AUC_CORRECTED = AUC$AUC * (1-AUC$FRACTION_QTL_FIXED)

### convert Pool-GPAS - GWAlpha_SNPwise COVARIATE from "" into an explicit "NONE"
levels(AUC$COVARIATE)[levels(AUC$COVARIATE)==""] = "NONE"

### rename column for consistency prior to merging
colnames(RMSE) = c("TRAIN_SAMPLE_ID", "TEST_SAMPLE_ID", "MODEL", "COVARIATE", "N_TRAIN", "N_TEST", "p0", "p1", "r2adj", "RMSE", "CORRELATION")
### convert Pool-GPAS - GWAlpha_SNPwise COVARIATE from NA into an explicit "NONE"
RMSE$COVARIATE = as.character(RMSE$COVARIATE)
RMSE$COVARIATE[is.na(RMSE$COVARIATE)] = "NONE"
RMSE$COVARIATE = as.factor(RMSE$COVARIATE)

### list training and test population *_SAMPLE_ID and *_POP
SAMPLE_ID_AND_POP_NAMES = unique(data.frame(TRAIN_SAMPLE_ID=AUC$TRAIN_SAMPLE_ID, TRAIN_POP=AUC$TRAIN_POP, TEST_SAMPLE_ID=AUC$TRAIN_SAMPLE_ID, TEST_POP=AUC$TRAIN_POP))

#### list all training and test population conbinations as observed in RMSE
CV_POP_NAMES = unique(merge(merge(RMSE[,1:2], SAMPLE_ID_AND_POP_NAMES[,1:2], by="TRAIN_SAMPLE_ID"), SAMPLE_ID_AND_POP_NAMES[,3:4], by="TEST_SAMPLE_ID")) ### TEST_SAMPLE_ID, TRAIN_SAMPLE_ID, TEST_SAMPLE_ID, TRAIN_POP, TEST_POP

### add *_POP (population names) into RMSE
RMSE_new = merge(RMSE, CV_POP_NAMES, by=c("TRAIN_SAMPLE_ID", "TEST_SAMPLE_ID"))

### clean-up
rm(RMSE); gc()

### merge RMSE_new with AUC duplicating AUC entries per test population across all the training populations
RMSE_AUC = merge(RMSE_new, AUC, by=c("TRAIN_SAMPLE_ID", "TRAIN_POP", "MODEL", "COVARIATE"), all=TRUE)

### add genotyping column (Indi-seq and Pool-seq)
RMSE_AUC$GENOTYPING = rep("INDIVIDUAL", times=nrow(RMSE_AUC))
RMSE_AUC$GENOTYPING[grep("GWAlpha", RMSE_AUC$MODEL)] = "POOL"
RMSE_AUC$GENOTYPING = as.factor(RMSE_AUC$GENOTYPING)

### correct SAMPLING_SCHEME to correspond to the training population's sampling scheme (and not the test population's)
RMSE_AUC$TRAIN_SAMPLE_ID = as.character(RMSE_AUC$TRAIN_SAMPLE_ID)
RMSE_AUC$SAMPLING = as.character(RMSE_AUC$SAMPLING)
within_and_across_names = matrix(unlist(strsplit(RMSE_AUC$TRAIN_SAMPLE_ID, "POP_")), ncol=2, byrow=TRUE)[,1] ### withins: ""; across: "MULTI"
RMSE_AUC$SAMPLING[within_and_across_names == ""] = "WITHIN"
RMSE_AUC$SAMPLING[within_and_across_names == "MULTI"] = "ACROSS"


### clean-up
rm(list=c("RMSE_new", "AUC", "SAMPLE_ID_AND_POP_NAMES", "CV_POP_NAMES")); gc()

# ### prepare population parameter columns
# PARAMETERS_COLUMNS = as.data.frame(matrix(unlist(rep(PARAMETERS, times=nrow(RMSE_AUC))), byrow=TRUE, nrow=nrow(RMSE_AUC)))
# colnames(PARAMETERS_COLUMNS) = colnames(PARAMETERS)

### estimate summary statistics for landscape phenotype means, variances and Fst by modelling as beta distribution
landscape_summstats = function(x, epsilon=1e-20, PLOT=FALSE) {
  # add a very small value so that we avoid infinitesimals in the likelihood function of the beta distribution below
  x = x + epsilon
  # define the beta-distribution likelihood function that we would like to minimize
  beta_pdf = function(par, data){
    # par = c(10, 0.01)
    # data = LANDSCAPE_PHENOTYPES$MEAN
    -sum(dbeta(data, shape1=par[1], shape2=par[2], log=TRUE))
  }
  # parameter (shapes1 and 2) estimation via maximum likelihood (minimization of the -loglik(dbeta|x) via Nelder-Mead optimization)
  optim_out = suppressWarnings(optim(par=c(1,1), fn=beta_pdf, data=x))
  # extract the parameters
  beta_shapes = optim_out$par
  # plotting histograms
  if (PLOT==TRUE) {
    par(mfrow=c(2,1))
    hist(x); legend("topright", legend=c(paste0("MIN=",round(min(x),4)), paste0("MAX=",round(max(x),4))))
    hist(rbeta(1000, shape1=beta_shapes[1], shape2=beta_shapes[2]))
  }
  return(beta_shapes)
}
# phenotype mean distribution`
PHENO_MEAN_BETADIST_SHAPES = as.data.frame(matrix(rep(landscape_summstats(LANDSCAPE_PHENOTYPES$MEAN), times=nrow(RMSE_AUC)), byrow=TRUE, nrow=nrow(RMSE_AUC)))
colnames(PHENO_MEAN_BETADIST_SHAPES) = c("SHAPE1", "SHAPE2")
# phenotype variance distribution`
PHENO_VAR_BETADIST_SHAPES = as.data.frame(matrix(rep(landscape_summstats(LANDSCAPE_PHENOTYPES$VAR), times=nrow(RMSE_AUC)), byrow=TRUE, nrow=nrow(RMSE_AUC)))
colnames(PHENO_VAR_BETADIST_SHAPES) = c("SHAPE1", "SHAPE2")
# genetic differentiation distribution
GENO_FST_BETADIST_SHAPES = as.data.frame(matrix(rep(landscape_summstats(as.vector(unlist(FST))), times=nrow(RMSE_AUC)), byrow=TRUE, nrow=nrow(RMSE_AUC)))
colnames(GENO_FST_BETADIST_SHAPES) = c("SHAPE1", "SHAPE2")

### merge everything
MERGED = data.frame(REPLICATE =                                 rep(PARAMETERS$rep, nrow(RMSE_AUC)),
                    EFFECTIVE_POP_SIZE =                        rep(PARAMETERS$nIndividuals, nrow(RMSE_AUC)),
                    NPOP =                                      rep(PARAMETERS$nPop, nrow(RMSE_AUC)),
                    NGEN =                                      rep(PARAMETERS$nGen, nrow(RMSE_AUC)),
                    NLOCI =                                     rep(PARAMETERS$nLoci, nrow(RMSE_AUC)),
                    NALLELES =                                  rep(PARAMETERS$nAlleles, nrow(RMSE_AUC)),
                    NQTL =                                      rep(PARAMETERS$nQTL, nrow(RMSE_AUC)),
                    NBGS =                                      rep(PARAMETERS$nBGS, nrow(RMSE_AUC)),
                    QTL_EFF_DIST =                              rep(PARAMETERS$allele_eff_model, nrow(RMSE_AUC)),
                    SELECTION_INTENSITY =                       rep(PARAMETERS$selection, nrow(RMSE_AUC)),
                    BGS_INTENSITY =                             rep(PARAMETERS$bg_selection, nrow(RMSE_AUC)),
                    MIGRATION_RATE =                            rep(PARAMETERS$migration, nrow(RMSE_AUC)),
                    QTL_GRADIENT =                              rep(PARAMETERS$GRADIENT, nrow(RMSE_AUC)),
                    PHENO_MEAN_BETADIST_SHAPE1 =                PHENO_MEAN_BETADIST_SHAPES$SHAPE1,
                    PHENO_MEAN_BETADIST_SHAPE2 =                PHENO_MEAN_BETADIST_SHAPES$SHAPE2,
                    PHENO_VAR_BETADIST_SHAPE1 =                 PHENO_VAR_BETADIST_SHAPES$SHAPE1,
                    PHENO_VAR_BETADIST_SHAPE2 =                 PHENO_VAR_BETADIST_SHAPES$SHAPE2,
                    GENO_FST_BETADIST_SHAPE1 =                  GENO_FST_BETADIST_SHAPES$SHAPE1,
                    GENO_FST_BETADIST_SHAPE2 =                  GENO_FST_BETADIST_SHAPES$SHAPE2,
                    NPOOLS_PER_POP =                            rep(PARAMETERS$NPOOLS, nrow(RMSE_AUC)),
                    NLIB_MAXIMUM =                              rep(PARAMETERS$NLIB, nrow(RMSE_AUC)),
                    SAMPLING_SCHEME =                           RMSE_AUC$SAMPLING,
                    GENOTYPING_SCHEME =                         RMSE_AUC$GENOTYPING,
                    TRAIN_POP =                                 RMSE_AUC$TRAIN_POP,
                    TEST_POP =                                  RMSE_AUC$TEST_POP,
                    NTRAIN =                                    RMSE_AUC$N_TRAIN,
                    NTEST =                                     RMSE_AUC$N_TEST,
                    MODEL =                                     RMSE_AUC$MODEL,
                    COVARIATE =                                 RMSE_AUC$COVARIATE,
                    p0 =                                        RMSE_AUC$p0,
                    p1 =                                        RMSE_AUC$p1,
                    r2adj =                                     RMSE_AUC$r2adj,
                    AUC =                                       RMSE_AUC$AUC,
                    FRACTION_QTL_FIXED =                        RMSE_AUC$FRACTION_QTL_FIXED,
                    AUC_CORRECTED =                             RMSE_AUC$AUC_CORRECTED,
                    BONFERRONI_5PERCENT_TRUE_POSITIVE =         RMSE_AUC$BONFERRONI_5PERCENT_TRUE_POSITIVE,
                    BONFERRONI_5PERCENT_FALSE_POSITIVE =        RMSE_AUC$BONFERRONI_5PERCENT_FALSE_POSITIVE,
                    BONFERRONI_5PERCENT_QTL_ID =                RMSE_AUC$BONFERRONI_5PERCENT_QTL_ID,
                    RMSE =                                      RMSE_AUC$RMSE,
                    CORRELATION =                               RMSE_AUC$CORRELATION
                    )
write.table(MERGED, file="GPAS_OUTPUT.csv", sep=",", quote=FALSE, row.names=FALSE, col.names=TRUE)

# ### Quick plots
# library(violinplotter)
# plot_me_df = data.frame(SAMPLING_STRATEGY=paste0(MERGED$SAMPLING_SCHEME, "-", MERGED$GENOTYPING_SCHEME),
#                         MODEL_COVARIATE=paste0(MERGED$MODEL, "-", MERGED$COVARIATE),
#                         AUC_CORRECTED=MERGED$AUC_CORRECTED,
#                         RMSE=MERGED$RMSE)
# png("GPAS_OUTPUT_metrics_X_samplings.png", width=1000, height=1000)
# par(mfrow=c(2,1))
# violinplotter(AUC_CORRECTED ~ SAMPLING_STRATEGY, data=plot_me_df, 
#               TITLE="Area under the ROC Plot\n(Correction Factor (1-P(fixed QTL)))",
#               XLAB="Sampling Strategies",
#               YLAB="AUC (corrected)",
#               VIOLIN_COLOURS=c("#EDF8B1", "#C7E9B4", "#7FCDBB", "#1D91C0"))
# violinplotter(RMSE ~ SAMPLING_STRATEGY, data=plot_me_df, 
#               TITLE="Root Mean Square Error",
#               XLAB="Sampling Strategies",
#               YLAB="RMSE",
#               VIOLIN_COLOURS=c("#EDF8B1", "#C7E9B4", "#7FCDBB", "#1D91C0"))
# dev.off()
# png("GPAS_OUTPUT_metrics_X_models.png", width=2000, height=1000)
# par(mfrow=c(2,1))
# violinplotter(AUC_CORRECTED ~ MODEL_COVARIATE, data=plot_me_df, 
#               TITLE="Area under the ROC Plot\n(Correction Factor (1-P(fixed QTL)))",
#               XLAB="GPAS Models",
#               YLAB="AUC (corrected)")
# violinplotter(RMSE ~ MODEL_COVARIATE, data=plot_me_df, 
#               TITLE="Root Mean Square Error",
#               XLAB="GPAS Models",
#               YLAB="RMSE")
# dev.off()

### STREAMLINE MERGED DATA:
### remove aut-cross-validation datapoints
train_pop_list = strsplit(as.character(MERGED$TRAIN_POP), ":")
test_pop_list = strsplit(as.character(MERGED$TEST_POP), ":")
idx = c()
for (i in 1:nrow(MERGED)) {
  if ( sum(unlist(train_pop_list[i]) %in% unlist(test_pop_list[i]))==0 ) {
    idx = c(idx, i)
  }
}
MERGED2 = droplevels(MERGED[idx, ])
rm(MERGED); gc()

### remove 2-pool validation sets (and incidentally all test populations with less than 2 individuals which don't really exist here)
DATA = droplevels(MERGED2[MERGED2$NTEST != 2, ])
rm(MERGED2); gc()

# ### a quick look
# violinplotter(AUC_CORRECTED ~ SAMPLING + MODEL*COVARIATE, data=DATA)
# violinplotter(RMSE ~ SAMPLING + MODEL*COVARIATE, data=DATA)

### streamlined dataframe for sampling - include only the relevant parameters and GPAS metrics
write.table(data.frame(SAMPLING_SCHEME=DATA$SAMPLING_SCHEME,
                       GENOTYPING_SCHEME=DATA$GENOTYPING_SCHEME,
                       TRAIN_POP=DATA$TRAIN_POP,
                       TEST_POP=DATA$TEST_POP,
                       NTRAIN=DATA$NTRAIN,
                       NTEST=DATA$NTEST,
                       MODEL=DATA$MODEL,
                       COVARIATE=DATA$COVARIATE,
                       AUC=DATA$AUC,
                       AUC_CORRECTED=DATA$AUC_CORRECTED,
                       BONFERRONI_5PERCENT_TRUE_POSITIVE=DATA$BONFERRONI_5PERCENT_TRUE_POSITIVE,
                       BONFERRONI_5PERCENT_FALSE_POSITIVE=DATA$BONFERRONI_5PERCENT_FALSE_POSITIVE,
                       BONFERRONI_5PERCENT_QTL_ID=DATA$BONFERRONI_5PERCENT_QTL_ID,
                       RMSE=DATA$RMSE,
                       CORRELATION=DATA$CORRELATION),
            file="GPAS_OUTPUT_STREAMLINED.csv", sep=",", row.names=FALSE, col.names=TRUE, quote=FALSE)

### Outputs
### (1) GPAS_OUTPUT.csv (comma-separated; header: |col01:REPLICATE|col02:EFFECTIVE_POP_SIZE|col03:NPOP|col04:NGEN|col05:NLOCI|col06:NALLELES|col07:NQTL|col08:NBGS|col09:QTL_EFF_DIST|col10:SELECTION_INTENSITY|\
###                                               |col11:BGS_INTENSITY|col12:MIGRATION_RATE|col13:QTL_GRADIENT|col14:PHENO_MEAN_BETADIST_SHAPES|col15:PHENO_VAR_BETADIST_SHAPES|col16:GENO_FST_BETADIST_SHAPES|col17:NPOOLS_PER_POP|col18:NLIB_MAXIMUM|col19:SAMPLING_SCHEME|col20:GENOTYPING_SCHEME|\
###                                               |col21:TRAIN_POP|col22:TEST_POP|col23:NTRAIN|col24:NTEST|col25:MODEL|col26:COVARIATE|col27:p0|col28:p1|col29:r2adj|col30:AUC|\
###                                               |col31:FRACTION_QTL_FIXED|col32:AUC_CORRECTED|col33:BONFERRONI_5PERCENT_TRUE_POSITIVE|col34:BONFERRONI_5PERCENT_FALSE_POSITIVE|col35:BONFERRONI_5PERCENT_QTL_ID|col36:RMSE|col37:CORRELATION|
### (2) GPAS_OUTPUT_STREAMLINED.csv (comma-separated; header: |col01:SAMPLING_SCHEME|col02:GENOTYPING_SCHEME|col03:TRAIN_POP|col04:TEST_POP|col05:NTRAIN|\
###                                                           |col06:NTEST|col07:MODEL|col08:COVARIATE|col09:AUC|col10:AUC_CORRECTED|\
###                                                           |col11:BONFERRONI_5PERCENT_TRUE_POSITIVE|col12:BONFERRONI_5PERCENT_FALSE_POSITIVE|col13:BONFERRONI_5PERCENT_QTL_ID|col14:RMSE|col15:CORRELATION|)
