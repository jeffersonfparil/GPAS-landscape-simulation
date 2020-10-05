####################
### ROC PLOTTING ###
####################
args = commandArgs(trailingOnly=TRUE)
# ### ONE POP TESTS
# args = c("POP_03_EMMAX_GRM.gwas", "QTL_SPEC.csv", "GENOME_SPEC.csv", "1.0", "TRUE")
# args = c("POP_03_GEMMA_GRM.gwas", "QTL_SPEC.csv", "GENOME_SPEC.csv", "1.0", "TRUE")
# args = c("POP_03_GCTA_GRM.gwas", "QTL_SPEC.csv", "GENOME_SPEC.csv", "1.0", "TRUE")
# args = c("POP_03_GCTA_STANDARDIZED.gwas", "QTL_SPEC.csv", "GENOME_SPEC.csv", "1.0", "TRUE")
# args = c("POP_03_GWAlpha_SNPwise.gwas", "QTL_SPEC.csv", "GENOME_SPEC.csv", "1.0", "TRUE")
# args = c("POP_03_GWAlpha_REML_LS_HIVERT.gwas", "QTL_SPEC.csv", "GENOME_SPEC.csv", "1.0", "TRUE")
# args = c("POP_03_GWAlpha_REML_LS_WEIRCOCK.gwas", "QTL_SPEC.csv", "GENOME_SPEC.csv", "1.0", "TRUE")
# args = c("POP_03_GWAlpha_REML_RR_WEIRCOCK.gwas", "QTL_SPEC.csv", "GENOME_SPEC.csv", "1.0", "TRUE")
# args = c("POP_03_GWAlpha_REML_GLMNET_WEIRCOCK.gwas", "QTL_SPEC.csv", "GENOME_SPEC.csv", "1.0", "TRUE")
# args = c("POP_03_GWAlpha_REML_LASSO_WEIRCOCK.gwas", "QTL_SPEC.csv", "GENOME_SPEC.csv", "1.0", "TRUE")
# ### MULTIPLOT TESTS
# args = c("MULTIPOP_8_9961816979469537538_EMMAX_GRM.gwas", "QTL_SPEC.csv", "GENOME_SPEC.csv", "1.0", "TRUE")
# args = c("MULTIPOP_8_9961816979469537538_GEMMA_GRM.gwas", "QTL_SPEC.csv", "GENOME_SPEC.csv", "1.0", "TRUE")
# args = c("MULTIPOP_8_9961816979469537538_GCTA_GRM.gwas", "QTL_SPEC.csv", "GENOME_SPEC.csv", "1.0", "TRUE")
# args = c("MULTIPOP_8_9961816979469537538_GCTA_STANDARDIZED.gwas", "QTL_SPEC.csv", "GENOME_SPEC.csv", "1.0", "TRUE")
# args = c("MULTIPOP_8_9961816979469537538_GWAlpha_SNPwise.gwas", "QTL_SPEC.csv", "GENOME_SPEC.csv", "1.0", "TRUE")
# args = c("MULTIPOP_8_9961816979469537538_GWAlpha_REML_LS_HIVERT.gwas", "QTL_SPEC.csv", "GENOME_SPEC.csv", "1.0", "TRUE")
# args = c("MULTIPOP_8_9961816979469537538_GWAlpha_REML_LS_WEIRCOCK.gwas", "QTL_SPEC.csv", "GENOME_SPEC.csv", "1.0", "TRUE")
# args = c("MULTIPOP_8_9961816979469537538_GWAlpha_REML_RR_WEIRCOCK.gwas", "QTL_SPEC.csv", "GENOME_SPEC.csv", "1.0", "TRUE")
# args = c("MULTIPOP_8_9961816979469537538_GWAlpha_REML_GLMNET_WEIRCOCK.gwas", "QTL_SPEC.csv", "GENOME_SPEC.csv", "1.0", "TRUE")
# args = c("MULTIPOP_8_9961816979469537538_GWAlpha_REML_LASSO_WEIRCOCK.gwas", "QTL_SPEC.csv", "GENOME_SPEC.csv", "1.0", "TRUE")
# for (i in 1:9) {
# args = c(paste0("POP_0", i,"_GCTA_STANDARDIZED.gwas"), "QTL_SPEC.csv", "GENOME_SPEC.csv", "1.0", "TRUE")

### Inputs
fname_input = args[1]       ### filename of the GWAS output
fname_qtl = args[2]         ### filename of the QTL file specification
fname_genome = args[3]      ### filenames of the simulated genome specification
LD_kb = as.numeric(args[4]) ### linkage block size assumed constant across the genome (in kilobases)
PLOT = as.logical(args[5])  ### save manhattan and ROC plots as png files (i.e. TRUE or FALSE)

### Output
### (1) ${prefix}_${model}_${kinship_id}.auc ### tab-delimited; HEADER: |col1:AUC|col2:FRACTION_OF_FIXED_QTL|col3:BONFERRONI_5PERCENT_TRUE_POSITIVE|col4:BONFERRONI_5PERCENT_FALSE_POSITIVE|col5:BONFERRONI_5PERCENT_QTL_ID|
### where:
### model=[EMMAX, GCTA, GEMMA]
### xiferp=$(echo $prefix | rev); pihsnik=$(echo ${kinship%.*} | rev)
### kinship_id=$(echo ${pihsnik%_$(echo $xiferp)*} | rev)

### NOTE: Since we've simulate each loci to be equidistant from each other within each chromosome then we have equal number of linked loci per QTL!
###       This means that our true positive rate is unbiased in terms of marker density variation across the genome!
###       And that our ROC AUC method described below should be applicable to empirical dataset when SNPs are LD trimmed or SNPs are equidistant per chromosome or scaffold.
### NOTE: The ROC plot null line (i.e. diagonal 1:1 line) assume that ALL THE QTL WERE PART OF THE ANALYSIS!
###       Therefore the curve will plunge below this null line if not all the QTL (including linked loci) were part of the analysis:
###       i.e. QTL (including linked loci) were filtered out because they were fixed to 0 or 1.

### Ouput prefix
split_vec = unlist(strsplit(fname_input, "[.]"))
prefix_output = paste(split_vec[1:(length(split_vec)-1)], collapse=".")

### Load GWAS output (Minimum info: CHROM, POS (or SNP=CHROM_POS), EFF, and PVAL)
if (grepl("EMMAX", fname_input)) {
  dat = read.delim(fname_input, sep="\t", header=FALSE)
  colnames(dat) = c("SNP", "EFF", "SE", "PVAL")
  dat$CHROM = as.factor(matrix(unlist(strsplit(as.character(dat$SNP), "_")), byrow=TRUE, ncol=2)[,1])
  dat$POS = as.numeric(matrix(unlist(strsplit(as.character(dat$SNP), "_")), byrow=TRUE, ncol=2)[,2])
} else if (grepl("GEMMA", fname_input)) {
  dat = read.delim(fname_input, sep="\t", header=TRUE)
  colnames(dat) = c("CHROM", "SNP", "POS", "NMISS", "ALT", "REF", "FREQ", "EFF", "SE", "LLH1", "LREML", "LMLE", "PWALD", "PLRT", "PVAL")
} else if (grepl("GCTA", fname_input)) {
  dat = read.delim(fname_input, sep="\t", header=TRUE)
  colnames(dat) = c("CHROM", "SNP", "POS", "REF", "ALT", "N", "FREQ", "EFF", "SE", "PVAL")
} else if (grepl("GWAlpha", fname_input)) {
  dat = read.delim(fname_input, sep=",", header=TRUE)
  colnames(dat) = c("CHROM", "POS", "REF", "FREQ", "EFF", "PVAL", "LOD")
  dat = dat[dat$CHROM != "Intercept", ]
}
### Calculate -log10(pvalues), and generate the streamlined GWAS output data frame
dat$LOG10PVAL = -log10(dat$PVAL + 1e-20)
GWAS = data.frame(CHROM=dat$CHROM, POS=dat$POS, LOG10PVAL=dat$LOG10PVAL)
# print(i)
# print(summary(dat$LOG10PVAL))
# }
rm(dat); gc()
### Load QTL information and merge with the genome information...
### ... to identify loci within the LD block of the QTL
QTL = aggregate(EFFECT ~ CHROM + POS, data=read.csv(fname_qtl), FUN=max) ### maximum effect of each QTL
QTL$SNP = paste0(QTL$CHROM, "_", QTL$POS)
GENOME = read.csv(fname_genome)
GENOME$SNP =paste0(GENOME$CHROM, "_", GENOME$POS)
GENOME$LOCI_IN_LD_WITH_QTL = rep(0, times=nrow(GENOME))
GENOME$QTL_SNP_ID = rep(NA, times=nrow(GENOME))
### Identify loci within LD_kb to the left and right of the QTL
for (snp in unique(QTL$SNP)) {
  # snp = unique(QTL$SNP)[1]
  chrom = QTL$CHROM[QTL$SNP==snp][1]
  pos = QTL$POS[QTL$SNP==snp][1]
  upper_lim_pos = pos + (1000*LD_kb)
  lower_lim_pos = max(c(1, pos - (1000*LD_kb)))
  idx =       (as.character(GENOME$CHROM) == as.character(chrom))
  idx = idx & (GENOME$POS <= upper_lim_pos)
  idx = idx & (GENOME$POS >= lower_lim_pos)
  GENOME$LOCI_IN_LD_WITH_QTL[idx] = 1
  GENOME$QTL_SNP_ID[idx] = snp
}
# table(GENOME$LOCI_IN_LD_WITH_QTL) ### to check if we really have the correct number of QTL
### OUTPUT1: Determine the proportion of QTL that were filtered out because they were fixed in the population
FRACTION_OF_FIXED_QTL = (nrow(QTL) - nrow(merge(QTL, GWAS, by=c("CHROM", "POS")))) / nrow(QTL)
### Merge the genome information with the GWAS results to generate the dataframe for ROC plotting
# ROC_DATA = merge(GENOME, GWAS, by=c("CHROM", "POS"), all=TRUE) ### include all loci in the genome
# ROC_DATA$LOG10PVAL[is.na(ROC_DATA$LOG10PVAL)] = 0
ROC_DATA = merge(GENOME, GWAS, by=c("CHROM", "POS"), all=FALSE) ### do not include loci that were filtered out
if (PLOT) {
  ### Manhattan plotting (for testing only because we're interested with the area under the ROC plot (AUC as the GWAS performance metric))
  png(paste0(prefix_output, "_MANHATTAN_PLOT.png"), width=2000, height=700)
  plot(x=1:nrow(ROC_DATA), y=ROC_DATA$LOG10PVAL, type="p", pch=20, col=rep(c("blue", "green"), times=ceiling(length(unique(ROC_DATA$CHROM))/2))[ROC_DATA$CHROM])
  points(x=c(1:nrow(ROC_DATA))[ROC_DATA$LOCI_IN_LD_WITH_QTL!=0], y=ROC_DATA$LOCI_IN_LD_WITH_QTL[ROC_DATA$LOCI_IN_LD_WITH_QTL!=0]*max(ROC_DATA$LOG10PVAL), pch=2, col="red")
  dev.off()
}
### Compute the area under the curve (AUC) of the true positive rate (TPR) as the false positive rate (FPR) increases
ROC_DATA$TRUE_POSITIVE = 1 * (ROC_DATA$LOCI_IN_LD_WITH_QTL==1) ### true positive loci == 1
ROC_DATA$FALSE_POSITIVE = 1 - ROC_DATA$TRUE_POSITIVE ### loci where if we were to continuously reduce the threshold will be false positves and so it's just all the loci less the actual QTL
ROC_DATA = ROC_DATA[order(ROC_DATA$LOG10PVAL, decreasing=TRUE), ] ### sort by decreasing -log10(pvalue) so that we are also sorting as we decrease the threshold (i.e. false positive rate should increase in this order)
ROC_DATA$FPR = cumsum(ROC_DATA$FALSE_POSITIVE) / sum(ROC_DATA$FALSE_POSITIVE) ### rate of increase of false positives conditional on the observed false positives (therefore may be zero which we will fix in the proceeding lines)
ROC_DATA$TPR = cumsum(ROC_DATA$TRUE_POSITIVE) / sum(ROC_DATA$TRUE_POSITIVE) ### rate of increase in true positives conditional on the observed true positives
ROC_DATA$FPR[is.na(ROC_DATA$FPR)] = 0.0 ### if no false positive loci were detected
ROC_DATA$TPR[is.na(ROC_DATA$TPR)] = 0.0 ### if no true positive loci were detected
TOTAL_FRACTION_OF_QTL_CAPTURED_IN_THE_GENOTYPE_DATA = length(unique(ROC_DATA$QTL_SNP_ID[ROC_DATA$TRUE_POSITIVE==1])) / length(unique(QTL$SNP)) ### Actual true positive rate: only one of the detected loci in LD for each of the QTL is counted
x = ROC_DATA$FPR
y = ROC_DATA$TPR
dx = x[2:length(x)] - x[1:length(x)-1]
mini_rectangles_areas = y[2:length(y)] * dx
### OUTPUT2: Area under the receiver operating characterstive plot (AUC ROC)
AUC = sum(mini_rectangles_areas)
if (PLOT) {
  ### ROC plotting (for testing only)
  png(paste0(prefix_output, "_ROC_PLOT.png"), width=700, height=700)
  plot(x=0:1, y=0:1, type="n", xlab="False Positive Rate", ylab="True Positive Rate")
  lines(x=0:1, y=0:1, lty=1, col="black")
  lines(x=x, y=y, lty=1, col="red")
  grid(lty=2, col="gray")
  Bonferroni_theshold = -log10(0.05/nrow(GWAS))
  abline(v=suppressWarnings(min(ROC_DATA$FPR[ROC_DATA$LOG10PVAL >= Bonferroni_theshold])), lty=3, col="blue")
  legend("top", legend=c("Null", paste0("ROC AUC = ", round(AUC,2)), paste0("Bonferroni Threshold (alpha=5%) = ", round(Bonferroni_theshold, 2))), lty=c(1,1,3), col=c("black", "red", "blue"))
  dev.off()
}
### ADD ON: true positive rate (TPR; TP/(TP+FN)); false discovery rate (FDR; FP/(FP+TP)); and the identities of QTL identified
###         at Bonferroni threshold (alpha=5%)
Bonferroni_theshold = -log10(0.05/nrow(GWAS))
DISCOVERIES = ROC_DATA[ROC_DATA$LOG10PVAL >= Bonferroni_theshold, ]
if (nrow(DISCOVERIES) > 0) {
  DISCOVERIES$SNP_ID = paste0(DISCOVERIES$CHROM, "_", DISCOVERIES$POS)
  BONFERRONI_5PERCENT_TRUE_POSITIVE = paste(DISCOVERIES$SNP_ID[DISCOVERIES$TRUE_POSITIVE==1], collapse=":")
  BONFERRONI_5PERCENT_FALSE_POSITIVE = paste(DISCOVERIES$SNP_ID[DISCOVERIES$TRUE_POSITIVE==0], collapse=":")
  BONFERRONI_5PERCENT_QTL_ID = paste(DISCOVERIES$QTL_SNP_ID[DISCOVERIES$TRUE_POSITIVE==1], collapse=":")
} else {
  BONFERRONI_5PERCENT_TRUE_POSITIVE = ""
  BONFERRONI_5PERCENT_FALSE_POSITIVE = ""
  BONFERRONI_5PERCENT_QTL_ID = ""
}
BONFERRONI_5PERCENT_TRUE_POSITIVE[BONFERRONI_5PERCENT_TRUE_POSITIVE==""] = NA
BONFERRONI_5PERCENT_FALSE_POSITIVE[BONFERRONI_5PERCENT_FALSE_POSITIVE==""] = NA
BONFERRONI_5PERCENT_QTL_ID[BONFERRONI_5PERCENT_QTL_ID==""] = NA

### Output
print(AUC)
print(FRACTION_OF_FIXED_QTL)
write.table(data.frame(AUC,
                       FRACTION_OF_FIXED_QTL,
                       BONFERRONI_5PERCENT_TRUE_POSITIVE,
                       BONFERRONI_5PERCENT_FALSE_POSITIVE,
                       BONFERRONI_5PERCENT_QTL_ID), file=paste0(prefix_output, ".auc"), row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
