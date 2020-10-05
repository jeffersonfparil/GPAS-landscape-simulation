########################
### ABC OPTIMIZATION ACROSS 0.0 TO 0.50 WEIGHTS ###
########################
### Input:
### (1) ABC_OPTIM_INPUT_SUMMSTATS.csv (sampling strategies across random combinatorial sensible samples and the corresponsing GPAS performance summaries;
###                                    comma-separated; header: col01:NSEQ||col02:ACROSS_INDI|col03:ACROSS_POOL|col04:WITHIN_INDI|col05:WITHIN_POOL|col06:AUC_MEAN|col07:AUC_CORRECTED_MEAN|col08:RMSE_MEAN|col09:CORRELATION_MEAN|col10:TRUE_POSITIVE_RATE|col11:FALSE_DISCOVERY_RATE|)
### Outputs:
### (1) ABC_OPTIM_OUTPUT-with_weights.csv (ABC optimized sampling scheme partitioning (Across_Indi vs Across_Pool vs Within_Indi vs Within_Pool) given a range of weights (0.00 --> 1.00) to the 2 metrics: AUC_CORRECTED and RMSE;
###                                        for determining the saturation point of the metrics (AUC and RMSE) as their corresponding weights change ;
###                                        this implies a step-wise optimization of sampling strategy (how many experiments devoted to ACROSS_INDI vs ACROSS_POOL vs WITHIN_INDI vs WITHIN_POOL?):
###                                                 (1) ABC optim across 0.0 to 1.0 weights range for AUC and RMSE (weight_auc + weight_rmse = 1)
###                                                 (2) Determine the weight combination at which AUC and RMSE saturate (maximal AUC and minimal RMSE)
###                                                 (3) Extract the resulting sampling strategy (experiment partitioning) at this optimal weight combination;
###                                        comma-separated; header: |col01:WEIGHT_AUC|col02:WEIGHT_RMSE|\
###                                                                 |col03:ACROSS_INDI|col04:ACROSS_POOL|col05:WITHIN_INDI|col06:WITHIN_POOL|\
###                                                                 |col07:AUC_MEAN|col08:AUC_CORRECTED_MEAN|col09:RMSE_MEAN|col10:CORRELATION_MEAN|col11:TRUE_POSITIVE_RATE|col12:FALSE_DISCOVERY_RATE|\
###                                                                 |col13:BZ_AUC (normalized then [0,1]-bounded)|col13:BZ_MINUS_ONE (normalized, then [0,1]-bounded, then subtracted from 1)|)
### (2) ABC_OPTIM_OUTPUT-with_weights-*_vs_*.png (weights vs metrics and params vs metrics scatter plots)
### (3) ABC_OPTIM_PRELIM-*.png (preliminary sanity checking plots)

### load input
args = commandArgs(trailingOnly=TRUE)
# args = c("ABC_OPTIM_INPUT_SUMMSTATS.csv")
fname_input = args[1]
dat = read.csv(fname_input, header=TRUE)

### load libraries
library(abc)

### define parameters and metrics
parameters = colnames(dat)[2:5]
metrics = colnames(dat)[6:11]
names_parameters = c("Across population sampling\nIndividual sequencing",
                     "Across population sampling\nPool sequencing",
                     "Within population sampling\nIndividual sequencing",
                     "Within population sampling\nPool sequencing")
colors_parameters = c("#EDF8B1",
                      "#C7E9B4",
                      "#7FCDBB",
                      "#1D91C0")
 names_metrics = c("Area under the ROC plot (AUC)",
                   "Corrected AUC",
                   "Root Mean Squared Error (RMSE)",
                   "Correlation",
                   "True Positive Rate",
                   "False Discovery Rate")
colors_metrics = c("#FB8072",
                   "#FDAE61",
                   "#A6D96A",
                   "#DE77AE",
                   "#9970AB",
                   "#80CDC1")

### sanity checks
png("ABC_OPTIM_PRELIM-paramsHist.png", width=900, height=700)
par(mfrow=c(2,2)) ### parameters
for (i in 1:length(parameters)){
  p = parameters[i]
  name_p = names_parameters[i]
  color_p = colors_parameters[i]
  eval(parse(text=paste0("hist(dat$", p, ", main='", name_p, "', xlab='", name_p, "', col='", color_p, "', bord=FALSE)")))
}
dev.off()

png("ABC_OPTIM_PRELIM-metricsHist.png", width=1200, height=700)
par(mfrow=c(2,3)) ### metrics
for (i in 1:length(metrics)){
  p = metrics[i]
  name_p = names_metrics[i]
  color_p = colors_metrics[i]
  eval(parse(text=paste0("hist(dat$", p, ", main='", name_p, "', xlab='", name_p, "', col='", color_p, "', bord=FALSE)")))
}
dev.off()

png("ABC_OPTIM_PRELIM-metricsXmetrics.png", width=1000, height=1000)
par(mfrow=c(4,4)) ### metrics bi-plots
for (i in 1:(length(metrics)-1)) {
  for (j in (i+1):length(metrics)) {
    x = eval(parse(text=paste0("dat$", metrics[i])))
    y = eval(parse(text=paste0("dat$", metrics[j])))
    names_xy = c(names_metrics[i], names_metrics[j])
    colors_xy = c(colors_metrics[i], colors_metrics[j])
    mod = lm(y ~ x)
    b0 = round(coef(mod)[1], 4)
    b1 = round(coef(mod)[2], 4)
    r2adj = summary(mod)$adj.r.squared
    plot(x, y, xlab=names_xy[1], ylab=names_xy[2], main=paste(names_xy, collapse=c("\nX\n")), type="p", pch=20, col=colors_xy)
    grid(col="grey")
    lines(x=c(min(x), max(x)), y=b0+(c(min(x), max(x))*b1), lty=2, lwd=2, col="black")
    legend("topleft", legend=c(paste0("R2adj=", round(r2adj*100), "%"),
                               paste0("Intercept=", b0),
                               paste0("Slope=", b1)))
  }
}
dev.off()

png("ABC_OPTIM_PRELIM-paramsXmetrics.png", width=1500, height=1100)
par(mfrow=c(length(parameters), length(metrics)), mar=c(5,5,6,2)) ### parameters x metrics
for (i in 1:length(parameters)) {
  for (j in 1:length(metrics)) {
    x = eval(parse(text=paste0("dat$", parameters[i])))
    y = eval(parse(text=paste0("dat$", metrics[j])))
    names_xy = c(names_parameters[i], names_metrics[j])
    colors_xy = c(colors_parameters[i], colors_metrics[j])
    mod = lm(y ~ x)
    b0 = round(coef(mod)[1], 4)
    b1 = round(coef(mod)[2], 4)
    r2adj = summary(mod)$adj.r.squared
    plot(x, y, xlab=names_xy[1], ylab=names_xy[2], main=paste(names_xy, collapse=c("\nX\n")), type="p", pch=20, col=colors_xy)
    grid(col="grey")
    lines(x=c(min(x), max(x)), y=b0+(c(min(x), max(x))*b1), lty=2, lwd=2, col="black")
    legend("topleft", legend=c(paste0("R2adj=", round(r2adj*100), "%"),
                               paste0("Intercept=", b0),
                               paste0("Slope=", b1)))
  }
}
dev.off()

### ABC with only 2 metrics: AUC_CORRECTED_MEAN and RMSE_MEAN
idx = c(2,3)
metrics = metrics[idx]
names_metrics = names_metrics[idx]
colors_metrics = colors_metrics[idx]

### NOTE: put both metrics on the same range so they can be comparable
### i.e standardardize then map to 0->1 range
## NO LONGER PERFORMED: mapping to 0->1 range for an easy weighting scheme
# Z_AUC = scale(dat$AUC_CORRECTED_MEAN, center=TRUE, scale=TRUE)
# Z_RMSE = scale(dat$RMSE_MEAN, center=TRUE, scale=TRUE)
# dat$BZ_AUC = (Z_AUC - min(Z_AUC,na.rm=TRUE)) / (max(Z_AUC,na.rm=TRUE) - min(Z_AUC,na.rm=TRUE)) ### COLUMN 11
# dat$ONE_MINUS_BZ_RMSE = 1 - ((Z_RMSE - min(Z_RMSE,na.rm=TRUE)) / (max(Z_RMSE,na.rm=TRUE) - min(Z_RMSE,na.rm=TRUE))) ### COLUMN 12
dat$Z_AUC = scale(dat$AUC_CORRECTED_MEAN, center=TRUE, scale=TRUE)
dat$Z_RMSE = scale(1 - dat$RMSE_MEAN, center=TRUE, scale=TRUE)

### ABC optimization putting 0->1 weights on the 2 metrics
### metric weights repeated twice:
###     (1) in the first half target_auc is fixed to 1 and target_rmse inreases from 0 to 1 as weight_rmse increases from 0.0 to 0.5
###     (2) in the final half target_rmse is fixed to 1 and target_auc inreases from 0 to 1 as weight_auc increases from 0.0 to 0.5
### this ensures we iterate across 0.0 to 1.0 weights for both metrics: AUC and RMSE, ...
### alternating the fixed target value between AUC and RMSE as the weights increase from 0.0 to 0.5...
### where at weight_auc = weight_rmse = 0.5: target_auc = target_rmse = 1
### and we can't just range the weight from 0 to 1 as one of the target values are fixed which means the other target value will be greater than 0 which is beyond the range of possible values for our taget values since we bounded
WEIGHTS = rep(seq(from=0.0, to=0.5, by=1/max(dat$NSEQ)), times=2)
OUT = matrix(0, nrow=length(WEIGHTS), ncol=2+ncol(dat))
### heuristically determine an appropriate tolerance value via 100-fold cross-validation
### columns 2:5 are the parameters (experiment partitionining across the 4 sampling schemes)
### columns 12:13 are the bounded standardized AUC and 1-bounded standardized RMSE
cv = tryCatch(
  cv4abc(param=dat[,2:5], sumstat=dat[,12:13], nval=100, tols=c(0.0001, 0.001, 0.01, 0.1), statistic="mean", method="rejection"),
  error=function(e){
    cv4abc(param=dat[,2:5], sumstat=dat[,12:13], nval=100, tols=c(0.001, 0.01), statistic="mean", method="rejection")
  })
tol = as.numeric(names(which.min(rowMeans(scale(summary(cv))))))
### (STEP 1/2) increase the weight of 1-f(RMSE) from 0.0 to 0.5, where the target value of f(AUC)=1 is fixed (at weight_rmse=0.5 target_auc=target_rmse=1)
for (i in 1:(length(WEIGHTS)/2)) {
  ### weights (weight_rmse increasing from 0.0 to 0.5)
  weight_auc = 1-WEIGHTS[i]
  weight_rmse = WEIGHTS[i]
  ### target values (target_auc is fixed since weight_rmse is increasing)
  target_auc = max(dat$Z_AUC, na.rm=TRUE) * 1
  target_rmse = max(dat$Z_RMSE, na.rm=TRUE) * (weight_rmse / (1-weight_rmse)) ### from: weight_rmse = target_rmse / (target_rmse + target_auc), where target_auc=1
  target_vector = c(target_auc,
                    target_rmse)
  ### optimize
  out = abc(target=target_vector, param=dat[,2:5], sumstat=dat[,12:13], tol=tol, method="rejection")
  DATA_SELECTED = droplevels(dat[out$region,])
  OUT[i, ] = c(WEIGHT_AUC=weight_auc, WEIGHT_RMSE=weight_rmse, colMeans(DATA_SELECTED))
}
### (STEP 2/2) increase the weight of f(AUC)) from 0.0 to 0.5, where the target value of 1-f(RMSE)=1 is fixed (at weight_auc=0.5 target_auc=target_rmse=1)
for (i in ((length(WEIGHTS)/2)+1):length(WEIGHTS)) {
  ### weights (weight_auc increasing from 0.0 to 0.5)
  weight_auc = WEIGHTS[i]
  weight_rmse = 1-WEIGHTS[i]
  ### target values (target_auc is fixed since weight_rmse is increasing)
  target_auc = max(dat$Z_AUC, na.rm=TRUE) *  (weight_auc / (1-weight_auc)) ### from: weight_auc = target_auc / (target_auc + target_rmse), where target_rmse=1
  target_rmse = max(dat$Z_RMSE, na.rm=TRUE) * 1
  target_vector = c(target_auc,
                    target_rmse)
  ### optimize
  out = abc(target=target_vector, param=dat[,2:5], sumstat=dat[,12:13], tol=tol, method="rejection")
  DATA_SELECTED = droplevels(dat[out$region,])
  OUT[i, ] = c(WEIGHT_AUC=weight_auc, WEIGHT_RMSE=weight_rmse, colMeans(DATA_SELECTED))
}

### prepare output
colnames(OUT) = c("WEIGHT_AUC", "WEIGHT_RMSE", colnames(dat))
OUT = as.data.frame(OUT)
fnames_prefix_output = "ABC_OPTIM_OUTPUT-with_weights"

### plot outputs
png(paste0(fnames_prefix_output, "-WEIGHTS_vs_METRICS.png"), width=2000, height=1200)
par(mfrow=c(2,2), mar=c(5,5,6,2), cex=1.5)
for (i in 1:length(metrics)){
  y = cbind(OUT$AUC_CORRECTED_MEAN, OUT$RMSE_MEAN)[,i] ### metrics
  x_weight_auc = OUT$WEIGHT_AUC ### AUC weights
  x_weight_rmse = OUT$WEIGHT_RMSE ### RMSE weights
  ylab = names_metrics[i]
  color = colors_metrics[i]
  plot(x=x_weight_auc, y=y, xlab="AUC Weights", ylab=ylab, main=paste0("AUC weights", "\nX\n", ylab), col=color, pch=20)
  plot(x=x_weight_rmse, y=y, xlab="RMSE Weights", ylab=ylab, main=paste0("RMSE weights", "\nX\n", ylab), col=color, pch=20)
}
dev.off()
png(paste0(fnames_prefix_output, "-PARAMS_vs_METRICS.png"), width=2000, height=1200)
par(mfrow=c(2,4), mar=c(5,5,6,2), cex=1.5)
for (i in 1:length(metrics)){
  for (j in 1:length(parameters)){
    y = cbind(OUT$AUC_CORRECTED_MEAN, OUT$RMSE_MEAN)[,i] ### metrics
    x = cbind(OUT$ACROSS_INDI, OUT$ACROSS_POOL, OUT$WITHIN_INDI, OUT$WITHIN_POOL)[,j] ### parameters
    ylab = names_metrics[i]
    xlab = names_parameters[j]
    colors = c(colors_metrics[i], colors_parameters[j])
    plot(x=x, y=y, xlab=xlab, ylab=ylab, main=paste0(xlab, "\nX\n", ylab), col=colors, pch=20)
  }
}
dev.off()

### write-out
write.table(OUT, file=paste0(fnames_prefix_output, ".csv"), sep=",", quote=FALSE, row.names=FALSE, col.names=TRUE)

# ### Alternative optimization with L-BFGS-B (limited memory box-constrained Broyden, Fletcher, Goldfarb and Shanno quasi-Newton optimisation algorithm, Byrd et al, 1995)
# ### NOTE: PROBLEM IS THE VERY LIMITED SEARCH SPACE BECAUSE WE ARE USING OBERVED (SIMULATED) DATA!!! MAY JUST NEED TO STICK WITH ABC OR USE THIS TYPE OF OPTIMIZATION ON THE FULL DATASET ACROSS SIMULATED LANDSCAPES!
# cost_function = function(par, weights, data, epsilon=0.1){
#   # par=c(1, 1, 2, 2); weights=c(0.5, 0.5); data=dat; epsilon=0.01
#   idx_across_indi = ((par[1]/sum(par)) - data$ACROSS_INDI) < epsilon
#   idx_across_pool = ((par[2]/sum(par)) - data$ACROSS_POOL) < epsilon
#   idx_within_indi = ((par[3]/sum(par)) - data$WITHIN_INDI) < epsilon
#   idx_within_pool = ((par[4]/sum(par)) - data$WITHIN_POOL) < epsilon
#   sub = data[idx_across_indi & idx_across_pool & idx_within_indi & idx_within_pool, ]
#   sub = droplevels(sub[complete.cases(sub), ])
#   if (nrow(sub)>0){
#     loss = 1 - ((mean(sub$BZ_AUC) * weights[1]) + (mean(sub$ONE_MINUS_BZ_RMSE) * weights[2]))
#   } else {
#     loss = Inf
#   }
#   return(loss)
# }
# weights = c(0.5, 0.5)
# initial_par = c(1,1,1,1)
# optim(initial_par, fn=cost_function, weights=weights, data=dat, method="L-BFGS-B", lower=c(0,0,0,0), upper=c(1,1,1,1))
