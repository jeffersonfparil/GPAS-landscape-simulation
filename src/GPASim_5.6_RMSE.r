######################################################
### RMSE AND GENOMIC PREDICTION SUMMARY STATISTICS ###
######################################################
### Predict phenotypes using the plink-derived polygenic score and
### the training population's linear predictors connecting the actual phenotypes with
### the plink-derived polygenic scores
### Input:
### (1) ${TRAIN_prefix}_${TRAIN_model}_${TRAIN_kinship_id}-${TEST_prefix}.gp - filename of polygenic risk score summing-up output from plink --score
### Outputs:
### (1) ${TRAIN_prefix}_${TRAIN_model}_${TRAIN_kinship_id}-${TEST_prefix}.p0p1 - filename of the predictors of the y_pred~polygenic linear model; header: |col1:p0|col2:p1|col3:r2adj|
###     - if and only if: ${TRAIN_prefix}==${TEST_prefix} that is the second step of the 2-step prediction model building: y_pred = p0 + p1*plink_polygenic_score; where plink_polygenic_score was derived using GWAS output
### (2) ${TRAIN_prefix}_${TRAIN_model}_${TRAIN_kinship_id}-${TEST_prefix}.rmse - filename of the RSME output; header: |col01:train_pop|col02:test_pop|col03:train_model|col04:n_train|col05:n_test|col06:p0|col07:p1|col08:r2adj|col09:RMSE|col10:CORR|

### input filename of polygenic risk score summing-up output from plink --score (.gp)
args = commandArgs(trailingOnly=TRUE)
# args = c("POP_01_EMMAX_GRM-POP_01.gp")
# args = c("MULTIPOP_16_5665778165567716942_GCTA_STANDARDIZED-MULTIPOP_16_5665778165567716942.gp")
# args = c("MULTIPOP_8_9961816979469537538_EMMAX_GRM-MULTIPOP_8_9961816979469537538.gp")
# args = c("POP_01_GWAlpha_SNPwise-POP_01.gp")
# args = c("POP_01_GWAlpha_SNPwise-POP_02.gp")
# args = c("MULTIPOP_8_9961816979469537538_GWAlpha_REML_GLMNET_WEIRCOCK-POP_02.gp")

if (sum(grepl("GWAlpha", unlist(strsplit(unlist(strsplit(sub(".gp", "", basename(args[1])), "-"))[1], "_")))) == 0) {
  print("@@@@@@@@@@@@@@@@@@@@@@@@@@@")
  print("Indi-GP cross-validation")
  print("@@@@@@@@@@@@@@@@@@@@@@@@@@@")
  ### training and testing population idenitifying names
  fname_scores = args[1]
  train_id = unlist(strsplit(unlist(strsplit(sub(".gp", "", basename(fname_scores)), "-"))[1], "_")) ### string-splitted pop+model names
  train_name = paste(train_id[1:(length(train_id)-2)], collapse="_") ### training population name
  train_model = train_id[length(train_id)-1] ### training model name
  train_covar = train_id[length(train_id)] ### training model name
  test_name = unlist(strsplit(sub(".gp", "", basename(fname_scores)), "-"))[2] ### validation or testing population name
  ### output filenames
  #(1) filename of the derived y_pred = p0 + p1*plink_polygenic_score for 2-step prediction using GWAS output header: |col1:p0|col2:p1|col3:r2adj|
  fname_p0p1 = paste0(dirname(fname_scores), "/", paste(train_id, collapse="_"), "-", train_name, ".p0p1")
  #(2) filename of the RSME output; header: |col01:train_pop|col02:test_pop|col03:train_model|col04:n_train|col05:n_test|col06:p0|col07:p1|col08:r2adj|col09:RMSE|col10:CORR|
  fname_rmse = paste0(dirname(fname_scores), "/", sub(".gp", "", basename(fname_scores)), ".rmse")
  ### load scores plink output
  dat = read.table(fname_scores, header=TRUE)
  ### solve for p0 and p1 for 2-step prediction using GWAS output for the training population only!
  if (train_name == test_name) {
    mod = lm(PHENO ~ SCORESUM, data=dat)
    p0 = coef(mod)[1]
    p1 = coef(mod)[2]
    r2adj = summary(mod)$adj.r.squared
    n_train = nrow(dat)
    write.table(data.frame(p0,p1,r2adj,n_train), file=fname_p0p1, col.names=TRUE, row.names=FALSE, quote=FALSE)
  }
  ### predict phenotypes
  if (!file.exists(fname_p0p1)){
    print(paste0("The file of coefficients of the PLINK-polygenic-score-to-phenotype model: '", fname_p0p1, "' does not exist! "))
    print((paste0("Please re-run with '", train_name, "' as both training and validation sets.")))
    quit(save="no")
  }
  p0 = read.table(fname_p0p1, header=TRUE)$p0
  p1 = read.table(fname_p0p1, header=TRUE)$p1
  r2adj = read.table(fname_p0p1, header=TRUE)$r2adj
  n_train = read.table(fname_p0p1, header=TRUE)$n_train
  n_test = nrow(dat)
  y_pred = p0 + (p1*dat$SCORESUM)
  ### extract summary statistics and write-out
  # plot(dat$PHENO, y_pred)
  dev = dat$PHENO - y_pred
  RMSE = sqrt((t(dev) %*% dev)/n_test)
  correlation = cor(dat$PHENO, y_pred)
} else {
  print("@@@@@@@@@@@@@@@@@@@@@@@@@@@")
  print("Pool-GP cross-validation")
  print("@@@@@@@@@@@@@@@@@@@@@@@@@@@")
  ### prepare names
  fname_gp = args[1]
  fname_rmse = paste0(dirname(fname_gp), "/", sub(".gp", "", basename(fname_gp)), ".rmse")
  ### load the Pool-GP ouput
  dat = read.table(fname_gp, header=TRUE)
  ### prepare output
  train_name = dat$TRAIN_name[1]
  test_name = dat$TEST_name[1]
  train_model = dat$TRAIN_model[1]
  train_covar = dat$TRAIN_rancovar[1]
  n_train = dat$n_train[1]
  n_test = dat$n_test[1]
  p0 = dat$p0[1]
  p1 = dat$p1[1]
  r2adj = dat$r2adj[1]
  dev = dat$y_true - dat$y_pred
  ### calculate RMSE and correlation
  RMSE = sqrt((t(dev) %*% dev)/n_test)
  correlation = cor(dat$y_true, dat$y_pred)
}
### Output
write.table(data.frame(TRAIN_NAME=train_name, TEST_NAME=test_name, TRAIN_MODEL=train_model, TRAIN_COVAR=train_covar, N_TRAIN=n_train, N_TEST=n_test, P0=p0, P1=p1, R2ADJ=r2adj, RMSE=RMSE, CORRELATION=correlation), file=fname_rmse, row.names=FALSE, col.names=TRUE, quote=FALSE)
