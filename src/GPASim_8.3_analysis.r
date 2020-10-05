################
### Analysis ###
################
### Questions:
### (1) How many populations to sample?
### (2) When to use Indi-seq or Pool-seq? (Note that Indi-seq: 384 individuals/population, while Pool-seq: 5 pools/population)
### (3) Which populations to select?
### Input:
### (1) MERGED_GPAS_OUTPUT.rds ### parsed RDS data
### Outputs:
### (1) FIXEF*.csv
### (2) HISTOGRAMS-*.*
### (3) LOG_REGRESSION-*.svg
### (4) MEANS_COEF_HSD_*.csv
### (5) RANEF-*.csv
### (6) SCATTERPLOTS-metric_biplots.png
### (7) VIOLINPLOTS-*.svg
### (8) WEIGHTED_OPTIMISATION-NPOP-TPR-FPR.png
### (9) 
### (10) 

args = commandArgs(trailingOnly=TRUE)
# args = c("MERGED_GPAS_OUTPUT_NPOP_PERFORMANCE.rds", "MERGED_VARIATION_SUMMSTATS.csv")
fname_input = args[1]
fname_var_summstats = args[2]

### load libraries
library(violinplotter)
library(lme4)

### load data
dat = readRDS(fname_input)

### remove MODEL = c("GWAlpha_GLMNET", "GWAlpha_LASSO", "GWAlpha_RIDGE")
dat = droplevels(dat[!(as.character(dat$MODE) %in% c("GWAlpha_GLMNET", "GWAlpha_LASSO", "GWAlpha_RIDGE")), ])

### list metric names and labels
names_metrics = c("RMSE_MEAN", "AUC_MEAN", "TRUE_POSITIVE_RATE", "FALSE_POSITIVE_RATE")
labels_metrics = c("Root Mean Square Error",
                   "Area under the ROC Plot",
                   "True Positive Rate",
                   "False Positive Rate")

### list explanatory variable names and labels
dat$NQTL =                as.numeric(dat$NQTL)
dat$SELECTION_INTENSITY = as.numeric(1.00 - dat$SELECTION_INTENSITY) ### reverted since the original variable is the distance of phenotypic values from the maximu possible phenotype value, i.e. lower values means higher selecetion intensities
dat$MIGRATION_RATE =      as.numeric(dat$MIGRATION_RATE)
dat$QTL_GRADIENT =        as.factor(dat$QTL_GRADIENT)
dat$GENOTYPING_SCHEME =   as.factor(dat$GENOTYPING_SCHEME)
dat$GPAS_MODEL =          as.factor(paste0(dat$MODEL, "-", dat$COVARIATE))
dat$NPOP_MERGED =         as.numeric(dat$NPOP_MERGED)
names_explain = c("NQTL", "SELECTION_INTENSITY", "MIGRATION_RATE", "QTL_GRADIENT", "GENOTYPING_SCHEME", "GPAS_MODEL", "NPOP_MERGED")
labels_explain = c("Number of Simulated QTL",
                   "Selection Intensity",
                   "Migration Rate",
                   "QTL Diffusion Gradient",
                   "Genotyping Method",
                   "GPAS Model",
                   "Number of Populations Sampled")

### log transform the numeric explanatory variables
### to reflect some of the non-linear relationships with the metrics we will soon observe in the violinplots
dat$LOG_NQTL = log10(dat$NQTL)
dat$LOG_SELE = log10(dat$SELECTION_INTENSITY)
dat$LOG_MIGR = log10(dat$MIGRATION_RATE)
dat$LOG_NPOP = log10(dat$NPOP_MERGED)
names_explain_log = c("LOG_NQTL", "LOG_SELE", "LOG_MIGR", "LOG_NPOP")
labels_explain_log = c("log10(Number of Simulated QTL)",
                       "log10(Selection Intensity)",
                       "log10(Migration Rate)",
                       "log10(Number of Training Populations)")

### plot metric distributions
svg("0.1_HISTOGRAMS-metric_distributions.svg", width=10, height=10)
par(mfrow=c(2,2))
for (i in 1:length(names_metrics)){
    # i = 1
    metric = eval(parse(text=paste0("dat$", names_metrics[[i]])))
    name_metric = labels_metrics[i]
    hist(metric, xlab=name_metric, main=name_metric)
}
dev.off()

### plot metric scatterplots
n = length(names_metrics)
n_row = round(sqrt(n*(n-1)/2))
n_col = ceiling((n*(n-1)/2)/n_row)
png("0.2_SCATTERPLOTS-metric_biplots.png", width=n_col*400, height=n_row*400)
par(mfrow=c(n_row, n_col))
for (i in 1:(n-1)){
    for (j in (i+1):n){
        x = eval(parse(text=paste0("dat$", names_metrics[i])))
        y = eval(parse(text=paste0("dat$", names_metrics[j])))
        x_lab = labels_metrics[i]
        y_lab = labels_metrics[j]
        plot(x, y, xlab=x_lab, ylab=y_lab, main=paste0(y_lab, "\nX\n", x_lab), pch=19, col=rgb(0.2, 0.2, 0.2, alpha=0.2))
    }
}
dev.off()

### plot violin plots 
### NOTE: AUC_MEAN and RMSE_MEAN are expected not to vary across NPOP_MERGED since these were simply averaged across intra-population trained GPAS models
### violin plots of non-log-transformed numeric explanatory variables
for (i in 1:length(names_metrics)){
    # i = 1
    metric = names_metrics[i]
    name_metric = labels_metrics[i]
    fname_svg = paste0("0.3_VIOLINPLOTS-", metric, ".svg")
    svg(fname_svg, width=15, height=15)
    eval(parse(text=paste0("violinplotter(", metric, " ~ ", paste(names_explain, collapse=" + "), ", data=dat)")))
    dev.off()
}
### compare the linear fit of untransformed against log-transformed numeric explanatory variables
names_numeric_explain = c(names_explain[(names_explain=="NQTL")|(names_explain=="SELECTION_INTENSITY")|(names_explain=="MIGRATION_RATE")|(names_explain=="NPOP_MERGED")],
                          names_explain_log)
for (i in 1:length(names_metrics)){
    # i = 1
    metric = names_metrics[i]
    name_metric = labels_metrics[i]
    fname_svg = paste0("0.4_VIOLINPLOTS-LOGX-", metric, ".svg")
    svg(fname_svg, width=15, height=7)
    par(mfrow=c(2,4))
    for (j in 1:length(names_numeric_explain)){
        # j = 1
        explain_var = names_numeric_explain[j]
        eval(parse(text=paste0("violinplotter(", metric, " ~ ", explain_var, ", data=dat, REGRESSX=TRUE)")))
    }
    dev.off()
}
### Results:
### (1) Notice the improved fit (i.e. higher R2-adjusted) for the log-transformed explanatory variables for most of them.
###     - We will be using log10(MIGRATION_RATE) and log10(NPOP_MERGED) for the subsequest linear regressions so that the relationships of the metrics with migration rate and the number of populations will be linear instead of non-linear.
### (2) AUC_MEAN is unaffected by NQTL, SELECTION_INTENSITY, and NPOP_MERGED (i.e. the R2_adjusted is less than 0.01, and the HSD grouping and ranking do not align with the regression line)
### (3) RMSE_MEAN is unaffccted by SELECTION_INTENSITY and NPOP_MERGED (...)
### (4) TRUE_POSITIVE_RATE is unaffected by SELECTION_INTENSITY (...)
### (5) FALSE_POSITIVE_RATE is unaffected by NQTL, SELECTION_INTENSITY, and MIGRATION_RATE (...)
### These violinplot results will be validated with linear models below.

### Before we can address question 1, we need to identify which GPAS models to use.
### Specifically, we need to identify best performing representative of Indi-GPAS and Pool-GPAS models, respectively.
### These will be used to subset the full dataset to include only the best models for subsequent analyses.
### This is to remove the confounding effects of different GPAS models.
### To do this we need to extract only the intra-population dataset, because we are intereseted in per population GPAS performances and not as affected by the number of populations sampled.
### Specifically, we will be using AUC_MEAN and RMSE_MEAN to assess GPAS performance.
dat_indi_pop = droplevels(dat[dat$NPOP_MERGED==1, ])
### Calculate the means, perform ANOVA, F-test, HSD mean comparison, and also extract the intercept-adjusted coefficients.
func_model_Ftest_HSD_rank = function(formula, var_name, data, alpha=0.05){
    # ### test input
    # formula = AUC_MEAN ~ LOG_NQTL * LOG_SELE * LOG_MIGR * QTL_GRADIENT * GPAS_MODEL
    # var_name = "GPAS_MODEL"
    # data = dat_indi_pop
    # alpha = 0.05
    ### rename var_name levels by replacing "-" with "_" in order to not confound the HSD grouping names
    eval(parse(text=paste0('data$', var_name, ' = sub("-", "_", data$', var_name, ')')))
    eval(parse(text=paste0('data$', var_name, ' = as.factor(data$', var_name, ')')))
    ### model fitting
    mod_aov = aov(formula=formula, data=data)
    ### analysis of variance and F-test to confirm at least one of the GPAS models is significantly different from the others
    anova_table = as.data.frame(anova(mod_aov))
    anova_table = cbind(rownames(anova_table), anova_table)
    rownames(anova_table) = NULL
    colnames(anova_table) = c("SV", "df", "SS", "MS", "F_value", "Sig")
    eval(parse(text=paste0('
        if (anova_table$Sig[anova_table$SV=="', var_name, '"] < alpha){
            print("Significant difference/s exist/s between ', var_name, ' levels!")
        } else {
            print("No significant difference exists between ', var_name, ' levels!")
        }')))
    ### computate the means per explanatory variable level
    y_name = unlist(strsplit(as.character(formula), " "))[2]
    eval(parse(text=paste0('means = eval(parse(text=paste0("aggregate(",  y_name, "~ ', var_name, ', data=data, FUN=mean)")))')))
    colnames(means) = c("VAR", "MEANS")
    means = means[order(means$MEANS, decreasing=TRUE), ]
    ### compute HSD pairwise comparisons
    hsd =  eval(parse(text=paste0('as.data.frame(TukeyHSD(mod_aov, which="', var_name, '", ordered=TRUE, conf.level=0.95)$', var_name, ')')))
    ### add "VAR_" string to allow for explanatory variable that are originally numeric to be easily set as list names
    factor_labels = matrix(paste0("VAR_", unlist(strsplit(rownames(hsd), "-"))), ncol=2, byrow=TRUE)
    hsd$factor1 = factor_labels[,1]
    hsd$factor2 = factor_labels[,2]
    factors_all = paste0("VAR_", as.character(means$VAR))
    ### initialize the list of HSD grouping of each response variable level
    GROUPING_LIST = eval(parse(text=paste0("list('VAR_", paste(as.character(means$VAR), collapse="'=c(), 'VAR_"), "'=c())")))
    ### generate the vector of letters and numbers for grouping
    letters_vector = c(letters, LETTERS, 1:(nrow(hsd)^2))
    ### iterate across response variable level
    letter_counter = 1
    for (f in factors_all){
        # f = factors_all[1]
        ### subset the current factor level
        subhsd = hsd[(hsd$factor1==f) | (hsd$factor2==f), ]
        ### identify the factor levels that are not significantly from the current factor level: f
        nonsigfactors = unique(c(subhsd$factor1[subhsd$p > 0.05], subhsd$factor2[subhsd$p > 0.05]))
        nonsigfactors = nonsigfactors[!(nonsigfactors %in% f)]
        ### define the current letter grouping
        letter_add = letters_vector[letter_counter]
        new_letter_bool = 0 ### for testing if we need a new letter
        ### iterate across non-significantly different factor levels to the current factor
        for (g in nonsigfactors){
            # g = nonsigfactors[1]
            f_letters = eval(parse(text=paste0("GROUPING_LIST$`", f, "`"))) ### currect factor grouping
            g_letters = eval(parse(text=paste0("GROUPING_LIST$`", g, "`"))) ### grouping of the non-siginificantly different factor level
            ### if we have all significantly different means at the start
            if (is.na(g)){
            eval(parse(text=paste0("GROUPING_LIST$`", f, "` = c(", "GROUPING_LIST$`", g, "`, '", letter_add, "')")))
            new_letter_bool = new_letter_bool + 1
            } else if ( !((sum(f_letters %in% g_letters)>0) | (sum(g_letters %in% f_letters)>0)) | is.null(f_letters) ) {
            ### test if the current factor level is the same as the non-siginificantly different factor level or if we are at the start
            eval(parse(text=paste0("GROUPING_LIST$`", g, "` = c(", "GROUPING_LIST$`", g, "`, '", letter_add, "')")))
            new_letter_bool = new_letter_bool + 1
            }
        }
        ### add the current letter grouping
        if ((new_letter_bool>0) | (length(nonsigfactors)==0)){
            eval(parse(text=paste0("GROUPING_LIST$`", f, "` = c(", "GROUPING_LIST$`", f, "`, '", letter_add, "')")))
            letter_counter = letter_counter + 1
        }
    }
    ### convert grouping list into a data.frame
    GROUPING_LIST = as.matrix(lapply(GROUPING_LIST, FUN=paste, collapse=""))
    GROUPING_LIST = data.frame(VAR=gsub("VAR_", "", as.character(rownames(GROUPING_LIST))), HSD_GROUPING=as.character(GROUPING_LIST[,1]))
    ### looking at intercept-adjusted coeffients not just the means
    coef_table = as.data.frame(coef(mod_aov))
    coef_table = cbind(rownames(coef_table), coef_table)
    rownames(coef_table) = NULL
    colnames(coef_table) = c("VAR", "COEFFICIENT")
    Intercept = coef_table$COEFFICIENT[grepl("Intercept", coef_table$VAR)]
    eval(parse(text=paste0('coef_table = coef_table[grepl("', var_name, '", coef_table$VAR) & !grepl(":", coef_table$VAR), ]')))
    eval(parse(text=paste0('coef_table$VAR = sub("', var_name, '", "", coef_table$VAR)')))
    eval(parse(text=paste0('coef_table = merge(data.frame(VAR=levels(data$', var_name, ')), coef_table, by="VAR", all=TRUE)')))
    coef_table$COEF_ADJ = coef_table$COEFFICIENT + Intercept
    coef_table$COEF_ADJ[is.na(coef_table$COEF_ADJ)] = Intercept
    coef_table = coef_table[order(coef_table$COEF_ADJ, decreasing=TRUE), c(1,3)]
    ### merge means, coefficients, and HSD grouping
    OUTPUT = merge(means, merge(coef_table, GROUPING_LIST, by="VAR"), by="VAR")
    ### output
    return(OUTPUT)
}
### The explanatory variables we're using are the log-transformed numeric variables, the QTL diffusion gradient levels, and the GPAS models.
### Note that the genotyping scheme is not included because the GPAS models are genotyping scheme-specific.
MEANS_COEF_HSD_GWAS = func_model_Ftest_HSD_rank(formula=AUC_MEAN~NQTL*SELECTION_INTENSITY*LOG_MIGR*QTL_GRADIENT*GPAS_MODEL, var_name="GPAS_MODEL", data=dat_indi_pop)
MEANS_COEF_HSD_GP   = func_model_Ftest_HSD_rank(formula=RMSE_MEAN~NQTL*SELECTION_INTENSITY*LOG_MIGR*QTL_GRADIENT*GPAS_MODEL, var_name="GPAS_MODEL", data=dat_indi_pop)
### Sort and identify the best model for Indi- and Pool-GPAS, respectively.
MEANS_COEF_HSD_GWAS[order(MEANS_COEF_HSD_GWAS$MEANS, decreasing=TRUE), ] ### higher AUC_MEAN : better performance
MEANS_COEF_HSD_GP[order(MEANS_COEF_HSD_GP$MEANS, decreasing=FALSE), ]  ### lower RMSE_MEAN : better performance
write.table(MEANS_COEF_HSD_GWAS, file="0.5_MEANS_COEF_HSD_GWAS.csv", sep=",", row.names=FALSE, col.names=TRUE, quote=FALSE)
write.table(MEANS_COEF_HSD_GP, file="0.5_MEANS_COEF_HSD_GP.csv", sep=",", row.names=FALSE, col.names=TRUE, quote=FALSE)
### Results:
### In terms of means and the corresponding HSD groupings:
###     - For Indi-GPAS: GEMMA-STANDARDIZED and GEMMA-GRM for both GWAS and GP.
###     - For Pool-GPAS: GWAlpha_SNPwise for GWAS, and GWAlpha_MIXED_REML-{HIVERT & WEIRCOCK} for GP.
### Note that in terms of the model coefficients:
###     - GWAlpha_SNPwise outperformns the Indi-GPAS models for both GWAS and GP!
### We will be using GEMMA-STANDARDIZED and GWAlpha_SNPwise, GPAS models.
### Subset the main dataset to include only these two models.
dat = droplevels(dat[(dat$GPAS_MODEL=="GEMMA-STANDARDIZED") | (dat$GPAS_MODEL=="GWAlpha_SNPwise-NONE"), ])
### Re-draw the violin plots of non-log-transformed numeric explanatory variables using the trimmed dataset
for (i in 1:length(names_metrics)){
    # i = 1
    metric = names_metrics[i]
    name_metric = labels_metrics[i]
    fname_svg = paste0("0.6_VIOLINPLOTS-TRIMMED", metric, ".svg")
    svg(fname_svg, width=10, height=8)
    eval(parse(text=paste0("violinplotter(", metric, " ~ ", paste(names_explain[names_explain!="GPAS_MODEL"], collapse=" + "), ", data=dat)")))
    dev.off()
}
############################################################################################################
############################################################################################################
############################################################################################################
### ADDRESSING THE FOLLOWING 3 MAIN QUESTIONS: 
### QUESTION 1: How many populations do we need to represent the landscape?
### QUESTION 2: Under what population-specific circumstances should we use Pool-seq over Indi-seq and vice-versa?  (Note that Indi-seq: 384 individuals/population, while Pool-seq: 5 pools/population)
### QUESTION 3: Which populations best represent the landscape under different circumstances?
############################################################################################################
############################################################################################################
############################################################################################################

##################
### QUESTION 1 ### How many populations do we need to represent the landscape?
##################
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
### Metric and explanatory variables ###
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
### We need to identify which metrics and explanatory variable to use to address question 1.
### We plot the 4 metrics as a function of the number of populations sampled across the landscape.
svg("1.1_VIOLINPLOTS-METRICS_vs_NPOP.svg", width=15, height=10)
par(mfrow=c(2,2))
shaved_dat_for_violinplot = droplevels(dat[dat$NPOP_MERGED!=4, ])
for (i in 1:length(names_metrics)){
    # i = 1
    metric = names_metrics[i]
    eval(parse(text=paste0("violinplotter(", metric, " ~ NPOP_MERGED, data=shaved_dat_for_violinplot)")))
}
dev.off()
### True positive and false negative rates are the only metrics varying across different number of populations samples.
### For question 1 we will only be using GWAS performance metrics: TPR and FPR.
### The explanatory variables we'll use are the number of population samples, the landscape-specific variables, and the genotyping scheme - which are now just one GPAS model each (i.e. GEMMA-STANDARDIXED and GWAlpha_SNPwise-NONE)
question_1_metric_names = c("TRUE_POSITIVE_RATE", "FALSE_POSITIVE_RATE")
question_1_metric_labels = c("True Positive Rate\n(Sensitivity or Power)",
                             "False Positive Rate\n(1 - Specificity or Fall-out)")
question_1_explain_names = c("NPOP_MERGED", "NQTL", "SELECTION_INTENSITY", "MIGRATION_RATE", "QTL_GRADIENT", "GENOTYPING_SCHEME")
question_1_explain_labels = c("Number of Populations Sampled",
                              "Number of Simulated QTL",
                              "Selection Intensity",
                              "Migration Rate",
                              "QTL Diffusion Gradient",
                              "Genotyping Method")
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
### Uncovering the relationship between TPR and FPR with the number of populations sampled ###
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
### plot TPR x log10(npop) in general and under different levels of each explanatory variable
func_plot_metric_x_npop = function(x, y, y_range=c(0,1), x_lab, y_lab, main_text="", deriv=FALSE){
    # ### test input
    # # y = dat$TRUE_POSITIVE_RATE
    # y = dat$FALSE_POSITIVE_RATE
    # x = dat$NPOP_MERGED
    # # y_range=c(0,1)
    # y_range=c(min(y,na.rm=TRUE), max(y,na.rm=TRUE))
    # y_lab = "True Positive Rate"
    # x_lab = "Number of Sampled Populations"
    # main_text=""; deriv=FALSE
    ### calculate means and sd
    agg_mean = aggregate(y ~ x, FUN=mean)
    agg_sd = aggregate(y ~ x, FUN=sd)
    ### scatterplot
    plot(x=c(min(x), max(x)), y=y_range, xlab=x_lab, ylab=y_lab, main=main_text, type="n")
    grid(lty=3, col="gray")
    ### remove duplicates for ease of plotting
    x_y = unique(paste0(x, "-", y))
    points(x=as.numeric(unlist(lapply(strsplit(x_y, "-"), FUN=function(x){x[1]}))),
           y=as.numeric(unlist(lapply(strsplit(x_y, "-"), FUN=function(x){x[2]}))),
           pch=20, col="#ccebc5")
    ### plot the standard deviation bars and the means
    arrows(x0=agg_mean$x, y0=agg_mean$y+agg_sd$y, y1=agg_mean$y-agg_sd$y, angle=90, code=3, lwd=2, length=0.1, col="#80b1d3")
    points(x=agg_mean$x, y=agg_mean$y, pch=19, col="#fb8072")
    ### log-regression line
    mod_log = lm(y ~ log10(x))
    newx = seq(min(agg_mean$x), max(agg_mean$x), length=100)
    newy = predict(mod_log, newdata=data.frame(x=newx))
    lines(x=newx, y=newy, lwd=2)
    points(x=agg_mean$x, y=agg_mean$y, pch=19, col="#fb8072") ### replot the means on top of the log-regression line
    ### plot first and second derivatives and the legend
    ### Note: R2_adjusted are based on mean fitted model (i.e. means were used instead of individual values)
    if (deriv==TRUE){
        ### first derivative of the log-regression (dy/dx = 1/ln(10)x for our y~log10(x) model)
        dy_dx = 1/(log(10)*newx)
        lines(x=newx, y = dy_dx, lty=1, lwd=2, col="#fdb462")
        ### second derivative of the log-regression (d2y/dx2 = -(1/ln(10)) * ((ln(x)+1)/(ln(x)^2 * x^2)))
        d2y_dx2 = -(1/log(10)) * ((log(newx)+1)/(log(newx)^2 * newx^2))
        lines(x=newx, y = d2y_dx2, lty=1, lwd=2, col="#e78ac3")
        legend("topright", legend=c((paste0("log10 regression (R2-adj=", round(summary(mod_log)$adj.r.squared*100), "%)")), "+/-1 standard deviation", "First derivative", "Second derivative"), lty=1, lwd=2, col=c("black", "#80b1d3", "#fdb462", "#e78ac3"))
    } else {
        legend("topright", legend=c((paste0("log10 regression (R2-adj=", round(summary(mod_log)$adj.r.squared*100), "%)")), "+/-1 standard deviation"), lty=1, lwd=2, col=c("black", "#80b1d3"))
    }
    return(0)
    # ### instead of a logarithmic fit, use a logistic fit
    # fx = function(par, data){
    #     y_pred = par[1] + ( par[2] / (par[3] + exp(data$x)) )
    #     return(y_pred)
    # }
    # ll = function(par, data){
    #     y_pred = fx(par, data)
    #     loss = sum((data$y - y_pred)^2)
    #     return(loss)
    # }
    # par = optim(par=c(0,0,0), fn=ll, data=data.frame(x=newx, y=newy))$par
    # ### plot line
    # y_logit = fx(par=par, data=data.frame(x=newx, y=newy))
    # lines(x=newx, y=y_logit, col="blue")
}
for (i in 1:length(question_1_metric_names)){
    # i = 1
    ### plot across all landscape-specific variables
    y_name = question_1_metric_names[i]
    x_name = "NPOP_MERGED"
    y_lab = question_1_metric_labels[i]
    x_lab = "Number of Sampled Populations"
    y = eval(parse(text=paste0("dat$", y_name)))
    x = eval(parse(text=paste0("dat$", x_name)))
    y_range = c(min(y,na.rm=TRUE), max(y,na.rm=TRUE))
    svg(paste0("1.2_LOG_REGRESSION-", y_name, "_x_NPOP-ALL.svg"), width=6, height=4)
    par(mar=c(7,7,5,2))
    func_plot_metric_x_npop(x=x, y=y, y_range=y_range, x_lab=x_lab, y_lab=y_lab, main_text="ALL", deriv=FALSE)
    dev.off()
    ### plot per landscape-specific explanatory variable
    for (j in 1:length(question_1_explain_names)){
        # j = 1
        explain = question_1_explain_names[j]
        print(explain)
        n_row = round(sqrt(nlevels(eval(parse(text=paste0("as.factor(dat$", explain, ")"))))+1))
        n_col = ceiling((nlevels(eval(parse(text=paste0("as.factor(dat$", explain, ")"))))+1)/n_row)
        svg(paste0("1.2_LOG_REGRESSION-", y_name, "_x_NPOP-", explain,".svg"), width=n_col*6, height=n_row*4)
        par(mfrow=c(n_row,n_col), mar=c(7,7,5,2))
        ### subset by the different levels of each explanatory variable
        for (k in 1:nlevels(eval(parse(text=paste0("as.factor(dat$", explain, ")"))))){
            # k = 1
            level = levels(eval(parse(text=paste0("as.factor(dat$", explain, ")"))))[k]
            sub_dat = droplevels(dat[eval(parse(text=paste0("as.factor(dat$", explain, ")")))==level, ])
            y = eval(parse(text=paste0("sub_dat$", y_name)))
            x = eval(parse(text=paste0("sub_dat$", x_name)))
            func_plot_metric_x_npop(x=x, y=y, y_range=y_range, x_lab=x_lab, y_lab=y_lab, main_text=paste0(explain, " = ", level), deriv=FALSE)
        }
        dev.off()
    }
}
### plot TPR and FPR vs the number of populations sampled for Indi-GPAS and Pool-GPAS in one graph
dat$FRAC_LAND = dat$NPOP_MERGED/max(dat$NPOP_MERGED, na.rm=TRUE)
dat_Indi = droplevels(dat[dat$GENOTYPING_SCHEME=="INDIVIDUAL",])
dat_Pool = droplevels(dat[dat$GENOTYPING_SCHEME=="POOL",])
vec_dat = list(dat_Indi, dat_Pool)
vec_colors = c("#ca0020", "#0571b0")

### test relationship under different landscape scenarios
dat_old = dat
idx_all = rep(TRUE, nrow(dat_old))
idx_qtl_010 = dat$NQTL == 010
idx_qtl_050 = dat$NQTL == 050
idx_qtl_100 = dat$NQTL == 100
idx_sel_0.50 = dat$SELECTION_INTENSITY == 0.50
idx_sel_0.90 = dat$SELECTION_INTENSITY == 0.90
idx_sel_0.95 = dat$SELECTION_INTENSITY == 0.95
idx_mig_1en4 = dat$MIGRATION_RATE == 1e-4
idx_mig_1en3 = dat$MIGRATION_RATE == 1e-3
idx_mig_1en2 = dat$MIGRATION_RATE == 1e-2
idx_grad_0 = dat$QTL_GRADIENT =="0"
idx_grad_1 = dat$QTL_GRADIENT =="1"
idx_grad_2 = dat$QTL_GRADIENT =="2"
list_idx = list(idx_all, idx_qtl_010, idx_qtl_050, idx_qtl_100, idx_sel_0.50, idx_sel_0.90, idx_sel_0.95, idx_mig_1en4, idx_mig_1en3, idx_mig_1en2, idx_grad_0, idx_grad_1, idx_grad_2)
list_lab = c("", "_qtl_010", "_qtl_050", "_qtl_100", "_sel_0.50", "_sel_0.90", "_sel_0.95", "_mig_1en4", "_mig_1en3", "_mig_1en2", "_grad_0", "_grad_1", "_grad_2")
for (i in 1:length(list_idx)){
    # i = 2
    idx = unlist(list_idx[i])
    lab = unlist(list_lab[i])
    print(lab)
    ### subset dat_old
    dat = dat_old[idx, ]
    svg(paste0("1.3_TPR_FPR_vs_NPOP", lab, ".svg"), width=7, height=10)
    par(mfrow=c(2,1), mar=c(5,5,2,2))
    xlim=c(min(c(0, dat$FRAC_LAND, na.rm=TRUE)), max(dat$FRAC_LAND, na.rm=TRUE))
    for (i in 1:length(question_1_metric_names)){
        # i = 1
        metric = question_1_metric_names[i]
        metric_label = question_1_metric_labels[i]
        agg_mean = eval(parse(text=paste0("aggregate(", metric, " ~ FRAC_LAND, FUN=mean, data=dat, na.rm=TRUE)")))
        agg_sd = eval(parse(text=paste0("aggregate(", metric, " ~ FRAC_LAND, FUN=sd, data=dat, na.rm=TRUE)")))
        ylim=c(max(0, min(agg_mean[,2]-2*max(agg_sd[,2]))), max(agg_mean[,2]+2*max(agg_sd[,2])))
        plot(x=xlim, y=ylim, xlab="Proportion of the Landscape Sampled", ylab=metric_label, type="n")
        grid(lty=2, col="gray")
        for (j in 1:length(vec_dat)){
            # j = 1
            d = vec_dat[[j]]
            color = vec_colors[j]
            agg_mean = eval(parse(text=paste0("aggregate(", metric, " ~ FRAC_LAND, FUN=mean, data=d, na.rm=TRUE)")))
            agg_sd = eval(parse(text=paste0("aggregate(", metric, " ~ FRAC_LAND, FUN=sd, data=d, na.rm=TRUE)")))
            points(x=agg_mean[,1], y=agg_mean[,2], pch=19, col=color)
            lines(x=agg_mean[,1], y=agg_mean[,2], lty=3, lwd=2, col=color)
            arrows(x0=agg_mean[,1], y0=agg_mean[,2]+agg_sd[,2], y1=agg_mean[,2]-agg_sd[,2], angle=90, code=3, lty=1, lwd=1, length=0.1, col=color)
        }
        legend("topleft", legend=c("Indi-GPAS", "Pool-GPAS"), col=vec_colors, pch=20, lty=1, lwd=2)
    }
    dev.off()

    # ### visualising the ratio between TPR and FPR ### not too informative - the separate plots above seem better
    # dat$RATIO = dat$TRUE_POSITIVE_RATE / dat$FALSE_POSITIVE_RATE
    # violinplotter(RATIO ~ FRAC_LAND, data=dat[dat$GENOTYPING_SCHEME=="POOL", ])

    ### TPR and FPR both increase logarithmically as the number of populations sampled increases among and within explanatory variables.
    ### We will disrearding the landscape-specific variations and initially just focus on the general trends and effects.
    ### This relationship is monotonous and lacks optimal points.
    ### It is worthwhile to note here that Pool-GPAS performs qualitatively better than Indi-GPAS in terms of minimising false positive rate!
    ### We want to have high TPR and low FPR.
    ### In order to find an optimum, we will be using weighted optimisation, i.e. y_weighted = w1*(TPR) - w2*(FPR); where 1 = w1 + w2
    ### HOWEVER!!
    ### FPR is minimised at the origin, i.e. we won't get any false positives if we don't perform any GWAS experiments in the first place! Hahaha!
    ### And so should we just ditch FPR as metric and just focus on TPR?
    ### Or come up with a better metric which accounts for both TPR and FPR?
    ### Use 3 different weights, i.e y_weighted = (w1*TPR) - (w2*FPR) - (w3*NPOP)?
    func_weighted_optimisation_npop_tpr_fpr = function(data){
        # ### test input
        # data = dat
        ### define the linear relationships of TPR with NPOP, and FPR with NPOP
        mod_tpr = lm(TRUE_POSITIVE_RATE ~ log10(NPOP_MERGED), data=data)
        mod_fpr = lm(FALSE_POSITIVE_RATE ~ log10(NPOP_MERGED), data=data)
        ### map out these relationships as the number of populations increase at regular intervals (as opposed to the irregular intervals we have, i.e. 1, 5, 10, 20, 25, ..., 100)
        x_pop = seq(from=min(data$NPOP_MERGED,na.rm=TRUE), to=max(data$NPOP_MERGED,na.rm=TRUE), length=100)
        x_tpr = predict(mod_tpr, newdata=data.frame(NPOP_MERGED=x_pop))
        x_fpr = predict(mod_fpr, newdata=data.frame(NPOP_MERGED=x_pop))
        ### map the resulting values into the 0 to range so that the three variables NPOP, TPR, and FPR are in the same range and the weightings are unbiased
        x_pop_binmap = (x_pop-min(x_pop,na.rm=TRUE))/(max(x_pop,na.rm=TRUE)-min(x_pop,na.rm=TRUE))
        x_tpr_binmap = (x_tpr-min(x_tpr,na.rm=TRUE))/(max(x_tpr,na.rm=TRUE)-min(x_tpr,na.rm=TRUE))
        x_fpr_binmap = (x_fpr-min(x_fpr,na.rm=TRUE))/(max(x_fpr,na.rm=TRUE)-min(x_fpr,na.rm=TRUE))
        ### define the weights such that the weights for NPOP, TPR and FPR sum up to one
        w_pop = c()
        w_tpr = c()
        w_fpr = c()
        for (i in seq(0,1,length=100)){
            for (j in seq(0,1,length=100)){
                for (k in seq(0,1,length=100)){
                    if (i + j + k == 1){
                        w_pop = c(w_pop, i)
                        w_tpr = c(w_tpr, j)
                        w_fpr = c(w_fpr, k)
                    }
                }
            }
        }
        ### calculate the weighted sum of NPOP, TPR and FPR
        ### find the maximum, and the corresponding values for NPOP, TPR and FPR
        y_max = c()
        optim_pop = c()
        optim_tpr = c()
        optim_fpr = c()
        for (i in 1:length(w_pop)){
            # i = 150
            y = (w_tpr[i]*x_tpr_binmap) - (w_fpr[i]*x_fpr_binmap) - (w_pop[i]*x_pop_binmap)
            # par(mfrow=c(2,2))
            # plot(x_tpr, y)
            # plot(x_fpr, y)
            # plot(x_pop, y)
            # hist(y)
            # max(y, na.rm=TRUE)
            y_max = c(y_max, max(y, na.rm=TRUE))
            optim_pop = c(optim_pop, x_pop[y==max(y, na.rm=TRUE)])
            optim_tpr = c(optim_tpr, x_tpr[y==max(y, na.rm=TRUE)])
            optim_fpr = c(optim_fpr, x_fpr[y==max(y, na.rm=TRUE)])
        }
        ### merge the weights, maximum weighted sum, and the correspoding optimum values for NPOP< TPR and FPR
        OPTIM = data.frame(COST_NPOP  = w_pop,
                        WEIGHT_TPR = w_tpr,
                        COST_FPR   = w_fpr,
                        Y_WEIGHTED_MAX = y_max,
                        OPTIM_NPOP = optim_pop,
                        OPTIM_TPR  = optim_tpr,
                        OPTIM_FPR  = optim_fpr)
        return(OPTIM)
    }
    ### Determine the optimum NPOP, TPR and FPR at different combinations of weights or costs to these variables.
    ### Perform these optimisations using the combined as well as seprately for the Indi-GPAS, and Pool-GPAS models.
    # OPTIM = func_weighted_optimisation_npop_tpr_fpr(data=dat)
    OPTIM_IndiGPAS = func_weighted_optimisation_npop_tpr_fpr(data=droplevels(dat[dat$GENOTYPING_SCHEME=="INDIVIDUAL",]))
    OPTIM_PoolGPAS = func_weighted_optimisation_npop_tpr_fpr(data=droplevels(dat[dat$GENOTYPING_SCHEME=="POOL",]))
    # write.table(OPTIM, file="1.4_WEIGHTED_OPTIMISATION-AVG.csv", sep=",", quote=FALSE, row.names=FALSE, col.names=TRUE)
    write.table(OPTIM_IndiGPAS, file=paste0("1.4_WEIGHTED_OPTIMISATION-IndiGPAS", lab, ".csv"), sep=",", quote=FALSE, row.names=FALSE, col.names=TRUE)
    write.table(OPTIM_PoolGPAS, file=paste0("1.4_WEIGHTED_OPTIMISATION-PoolGPAS", lab, ".csv"), sep=",", quote=FALSE, row.names=FALSE, col.names=TRUE)
    ### Results:
    ### [1,1] At 0.00 to 0.20 cost to sampling populations, we get the most flexibility in the number of population samples.
    ###       At this cost range, the weight towards TPR and cost towards FPR, are the determining factors in terms of how many populations we want to sample.
    ### [1,2] The monotonic upward sloping linear relationship between TPR and FPR suggests that a subjective decision needs to be made to find the optimum number of populations,
    ###       while having high TPR and low FPR.
    ### [2,1] The optimum TPR increases logarithmically as the number of population samples increases.
    ### [2,2] The same relationship for optimum FPR is observed. The only difference is the range of y-axis values, i.e. smaller for FPR comapred to TPR.
    ### These graphs serve as guides for selecting the number of populations that a GPAS project can afford and
    ### the corresponding expected TPR and FPR, averaged across and for each of the best Indi-GPAS and Pool-GPAS models,
    ### as well as across the different landscape circumstances we simulated - which we hope have realistic range.
    ### Comparing Indi-GPAS and Pool-GPAS models, we can see that Indi-GPAS outperforms Pool-GPAS in terms of TPR but at the cost of higher FPR!
    ### In line with this, we can ask the second question.
    # ### But before moving to the second question, let us assess the results of
    # ### regressing TPR and FPR across increasing population samples 
    # ### under each level of the landscape-specific explanatory variables.
    # ### i.e paste0("LOG_REGRESSION-", metric, "_x_NPOP-", explain,".svg")
    # ### Different levels of the landscape-specific variables affect the rate at which TPR and FPR increase with increasing number of populations sampled.
    # ### Let us assess this with linear mixed modelling.
    # # mod_LMM = lmer(TRUE_POSITIVE_RATE ~ LOG_NQTL + LOG_SELE + LOG_MIGR + QTL_GRADIENT + GENOTYPING_SCHEME + LOG_NPOP +
    # for (metric in question_1_metric_names){
    #     # metric = question_1_metric_names[1]
    #     eval(parse(text=paste0("mod_LMM = lmer(", metric, " ~ NQTL + SELECTION_INTENSITY + LOG_MIGR + QTL_GRADIENT + GENOTYPING_SCHEME + LOG_NPOP +
    #                 (0+LOG_NPOP|NQTL)+
    #                 (0+LOG_NPOP|SELECTION_INTENSITY)+
    #                 (0+LOG_NPOP|MIGRATION_RATE) +
    #                 (0+LOG_NPOP|GENOTYPING_SCHEME), data=dat)")))
    #     # fixef(mod_LMM)
    #     # ranef(mod_LMM)
    #     ### Fixed effects
    #     FIXEF = data.frame(fixef(mod_LMM))
    #     FIXEF = cbind(rownames(FIXEF), FIXEF); rownames(FIXEF) = NULL
    #     colnames(FIXEF) = c("Variable", "Fixed Effects")
    #     write.table(FIXEF, file=paste0("1.4_FIXEF-", metric, ".csv"), sep=",", quote=TRUE, row.names=FALSE, col.names=TRUE)
    #     ### Random slopes, TPR/NPOP under the different levels of the explanatory variables
    #     ### and plot line graphs of npop X metric where each line correspond to a factor level
    #     svg(paste0("1.4_RANDOM_SLOPES-", metric, ".svg"), width=10, height=10)
    #     par(mfrow=c(2,2))
    #     for (explain in names(ranef(mod_LMM))){
    #         # explain = names(ranef(mod_LMM))[1]
    #         RANEF = eval(parse(text=paste0("data.frame(ranef(mod_LMM)$", explain, ")")))
    #         RANEF = cbind(rownames(RANEF), RANEF); rownames(RANEF) = NULL
    #         colnames(RANEF) = c(explain, "SLOPE_TPR_OVER_NPOP")
    #         write.table(RANEF, file=paste0("1.4_RANEF-", metric, "-", explain, ".csv"), sep=",", quote=TRUE, row.names=FALSE, col.names=TRUE)
    #         plot(x=c(0,2), y=c(min(c(0, min(RANEF[,2]))),max(RANEF[,2])), type="n", xlab="log10(Number of Populations Sampled)", ylab=paste0("Effect on ", metric), main=explain, las=2)
    #         grid()
    #         for (i in 1:nrow(RANEF)){
    #             # i = 1
    #             lines(x=c(0,2), y=c(0,RANEF[i,2]), lty=1, lwd=2, col=c("#d73027", "#fdae61", "#31a354")[i])
    #         }
    #         legend("topleft", legend=RANEF[,1], lty=1, lwd=2, col=c("#d73027", "#fdae61", "#31a354"))
    #     }
    #     dev.off()
    # }
    # ### Reuslts:
    # ### NQTL: slope TPR/NPOP decreases as the number of QTL increases, i.e. it is less important to increase the number of samples if the number of QTL are predicted to be large since we expect diminishing returns because of the property of the simulated polygenic traits, i.e. some loci can never be detected because the effects are too small at the given heritability of 50% and fixed samples sizes of 384 (4X96-well plates for Indi-GPAS and 5 pools for Pool-GPAS.
    # ### SELECTION_INTENSITY: slope TPR/NPOP increases as the selection intensity increases, i.e. increasing the number of samples populations increases in importance if we expect the selection pressures to be high
    # ### MIGRATION_RATE: slope TPR/NPOP decreases as the migration rate increases, i.e. increasing the number of population samples is less important if we expect migration rate to be high
    # ### GENOTYPING_SCHEME: slope TPR/NPOP is lower in Pool-GPAS than Indi-GPAS, i.e. increasing the number of population samples is less important when using Pool-GPAS compared with Indi-GPAS, however, Indi-GPAS performs better than Pool-GPAS

    ### "Back of the envelope calculations" to test different Indi-seq:Pool-seq populations per sequencing capacity ratio
    # OPTIM_IndiGPAS = read.csv("1.3_WEIGHTED_OPTIMISATION-IndiGPAS.csv")
    # OPTIM_PoolGPAS = read.csv("1.3_WEIGHTED_OPTIMISATION-PoolGPAS.csv")
    ### list of ratios of the number of populations that can be characterised with Indi-seq
    ### for every population characterised with Pool-seq (assumes more populations can be characterised with Pool-seq than Indi-seq)
    tpr_indi = c()
    tpr_pool = c()
    fpr_indi = c()
    fpr_pool = c()
    vec_ratio = c()
    vec_npop_indi = c()
    vec_npop_pool = c()
    for (i in 1:100){
        for (j in 1:100){
            if ((i+j <= 100) & (i<=j)) {
                vec_npop_indi = c(vec_npop_indi, i)
                vec_npop_pool = c(vec_npop_pool, j)
            }
        }
    }
    for (i in 1:length(vec_npop_indi)){
        # i = 345
        npop_indi = vec_npop_indi[i]
        npop_pool = vec_npop_pool[i]
        vec_ratio = c(vec_ratio, npop_indi/npop_pool)
        tpr_indi = c(tpr_indi, mean(OPTIM_IndiGPAS$OPTIM_TPR[OPTIM_IndiGPAS$OPTIM_NPOP==npop_indi]))
        tpr_pool = c(tpr_pool, mean(OPTIM_PoolGPAS$OPTIM_TPR[OPTIM_PoolGPAS$OPTIM_NPOP==npop_pool]))
        fpr_indi = c(fpr_indi, mean(OPTIM_IndiGPAS$OPTIM_FPR[OPTIM_IndiGPAS$OPTIM_NPOP==npop_indi]))
        fpr_pool = c(fpr_pool, mean(OPTIM_PoolGPAS$OPTIM_FPR[OPTIM_PoolGPAS$OPTIM_NPOP==npop_pool]))
    }
    RATIO_DF = data.frame(RATIO=vec_ratio, NPOP_INDI=vec_npop_indi, NPOP_POOL=vec_npop_pool, TPR_INDI=tpr_indi, TPR_POOL=tpr_pool, FPR_INDI=fpr_indi, FPR_POOL=fpr_pool)
    RATIO_DF = RATIO_DF[order(RATIO_DF$RATIO, decreasing=FALSE), ]
    RATIO_DF$DIFF_TPR_Indi_Less_Pool = RATIO_DF$TPR_INDI - RATIO_DF$TPR_POOL
    RATIO_DF$DIFF_FPR_Indi_Less_Pool = RATIO_DF$FPR_INDI - RATIO_DF$FPR_POOL
    Selected_range = RATIO_DF[(RATIO_DF$DIFF_TPR_Indi_Less_Pool<0) & (RATIO_DF$DIFF_FPR_Indi_Less_Pool>0),]
    min_thresh = min(Selected_range$RATIO)
    max_thresh = max(Selected_range$RATIO)
    svg(paste0("1.5_TPR_FPR_DIFFERENCES_IndiGPAS_vs_PoolGPAS", lab, ".svg"), width=7, height=9)
        par(mfrow=c(3,1), mar=c(5,5,2,2))
        # plot(x=RATIO_DF$RATIO, y=RATIO_DF$DIFF_TPR_Indi_Less_Pool, type="l", lwd=3, col=rgb(0.5,0.5,0.5),
        #      xlab="Ratio of the Number of Populations Characterised with Indi-seq\nfor every Population Characterised with Pool-seq",
        #      ylab="TPR of Indi-GPAS Minus FPR of Pool-GPAS", las=2)
        plot(x=RATIO_DF$RATIO, y=RATIO_DF$DIFF_TPR_Indi_Less_Pool, type="n",
            xlab="Ratio of the Number of Populations Characterised with Indi-seq\nfor every Population Characterised with Pool-seq",
            ylab="TPR of Indi-GPAS Minus TPR of Pool-GPAS", las=2)
        col_vec = c("#253494", "#74add1", "#edf8b1", "#fdae61", "#d73027")
        for (i in 1:nrow(RATIO_DF)){
            if ((RATIO_DF$NPOP_POOL[i] > 0) & (RATIO_DF$NPOP_POOL[i] <= 20)){
                colour = col_vec[1]
            } else if ((RATIO_DF$NPOP_POOL[i] > 20) & (RATIO_DF$NPOP_POOL[i] <= 40)){
                colour = col_vec[2]
            } else if ((RATIO_DF$NPOP_POOL[i] > 40) & (RATIO_DF$NPOP_POOL[i] <= 60)){
                colour = col_vec[3]
            } else if ((RATIO_DF$NPOP_POOL[i] > 60) & (RATIO_DF$NPOP_POOL[i] <= 80)){
                colour = col_vec[4]
            } else {
                colour = col_vec[5]
            }
            points(x=RATIO_DF$RATIO[i], y=RATIO_DF$DIFF_TPR_Indi_Less_Pool[i], pch=20, col=colour)
        }
        abline(h=0.0, lty=1)
        abline(v=min_thresh, lwd=2, col="red")
        abline(v=max_thresh, lwd=2, col="red")
        grid()
        legend("bottomright",
            legend=c("Fraction of the Landscape\nCharacterised with Pool-seq:",
                    " 0 - 20%",
                    "21 - 40%",
                    "41 - 60%",
                    "61 - 80%",
                    "81 - 100%"),
            pch=19, col=c("white", col_vec)
            )
        # plot(RATIO_DF$RATIO, y=RATIO_DF$DIFF_FPR_Indi_Less_Pool, type="l", lwd=3, col=rgb(0.5,0.5,0.5),
        #     xlab="Ratio of the Number of Populations Characterised with Indi-seq\nfor every Population Characterised with Pool-seq",
        #     ylab="FPR of Indi-GPAS Minus FPR of Pool-GPAS", las=2)
        plot(x=RATIO_DF$RATIO, y=RATIO_DF$DIFF_FPR_Indi_Less_Pool, type="n",
            xlab="Ratio of the Number of Populations Characterised with Indi-seq\nfor every Population Characterised with Pool-seq",
            ylab="FPR of Indi-GPAS Minus FPR of Pool-GPAS", las=2)
        col_vec = c("#253494", "#74add1", "#edf8b1", "#fdae61", "#d73027")
        for (i in 1:nrow(RATIO_DF)){
            if ((RATIO_DF$NPOP_POOL[i] > 0) & (RATIO_DF$NPOP_POOL[i] <= 20)){
                colour = col_vec[1]
            } else if ((RATIO_DF$NPOP_POOL[i] > 20) & (RATIO_DF$NPOP_POOL[i] <= 40)){
                colour = col_vec[2]
            } else if ((RATIO_DF$NPOP_POOL[i] > 40) & (RATIO_DF$NPOP_POOL[i] <= 60)){
                colour = col_vec[3]
            } else if ((RATIO_DF$NPOP_POOL[i] > 60) & (RATIO_DF$NPOP_POOL[i] <= 80)){
                colour = col_vec[4]
            } else {
                colour = col_vec[5]
            }
            points(x=RATIO_DF$RATIO[i], y=RATIO_DF$DIFF_FPR_Indi_Less_Pool[i], pch=20, col=colour)
        }
        abline(h=0.0, lty=1)
        abline(v=min_thresh, lwd=2, col="red")
        abline(v=max_thresh, lwd=2, col="red")
        grid()
        legend("bottomright", legend=c(paste0("Min. Ratio = ", round(min_thresh,2)), paste0("Max. Ratio = ", round(max_thresh,2))))
        plot(x=Selected_range$NPOP_INDI, y=Selected_range$NPOP_POOL,
            xlab="Number of Populations Characterised with Indi-seq",
            ylab="Number of Populations Characterised with Pool-seq",
            pch=19, col=rgb(0.3,0.3,0.3,0.5))
        grid()
        legend("bottomright", legend=c(paste0("Min. Indi = ", min(Selected_range$NPOP_INDI)), paste0("Max. Indi = ", max(Selected_range$NPOP_INDI)),
                                    paste0("Min. Pool = ", min(Selected_range$NPOP_POOL)), paste0("Max. Pool = ", max(Selected_range$NPOP_POOL)) ))
        # par(mfrow=c(2,1), mar=c(5,5,2,2))
        # ### prep df
        # RATIO_DF = RATIO_DF[order(RATIO_DF$RATIO, decreasing=FALSE), ]
        # RATIO_DF$FOR_RIBBON = rep(RATIO_DF$RATIO[round(seq(1, nrow(RATIO_DF), length=10))], each=nrow(RATIO_DF)/10)
        # ### TPR
        # df_mean = aggregate(DIFF_TPR_Indi_Less_Pool ~ RATIO, data=RATIO_DF, FUN=mean)
        # df_min = aggregate(DIFF_TPR_Indi_Less_Pool ~ RATIO, data=RATIO_DF, FUN=min)
        # df_max = aggregate(DIFF_TPR_Indi_Less_Pool ~ RATIO, data=RATIO_DF, FUN=max)
        # x = c(df_max[,1], rev(df_min[,1]))
        # y = c(df_max[,2], rev(df_min[,2]))
        # plot(x, y, type="n", xlab="Ratio of the Number of Populations Characterised with Indi-seq\nfor every Population Characterised with Pool-seq",
        #                     ylab="FPR of Indi-GPAS Minus FPR of Pool-GPAS", las=2)
        # polygon(x, y, border="NA", col="gray")
        # lines(x=df_mean[,1], y=df_mean[,2])
        # abline(h=0.0, lty=1)
        # abline(v=min_thresh, lwd=2, col="red")
        # abline(v=max_thresh, lwd=2, col="red")
        # grid()
        # legend("bottomright", legend=c(paste0("Min. Ratio = ", round(min_thresh,2)), paste0("Max. Ratio = ", round(max_thresh,2))))
        # ### FPR
        # df_mean = aggregate(DIFF_FPR_Indi_Less_Pool ~ RATIO, data=RATIO_DF, FUN=mean)
        # df_min = aggregate(DIFF_FPR_Indi_Less_Pool ~ RATIO, data=RATIO_DF, FUN=min)
        # df_max = aggregate(DIFF_FPR_Indi_Less_Pool ~ RATIO, data=RATIO_DF, FUN=max)
        # x = c(df_max[,1], rev(df_min[,1]))
        # y = c(df_max[,2], rev(df_min[,2]))
        # plot(x, y, type="n", xlab="Ratio of the Number of Populations Characterised with Indi-seq\nfor every Population Characterised with Pool-seq",
        #                     ylab="FPR of Indi-GPAS Minus FPR of Pool-GPAS", las=2)
        # polygon(x, y, border="NA", col="gray")
        # lines(x=df_mean[,1], y=df_mean[,2])
        # abline(h=0.0, lty=1)
        # abline(v=min_thresh, lwd=2, col="red")
        # abline(v=max_thresh, lwd=2, col="red")
        # grid()
    dev.off()
    ### Results:
    ### Pool-GPAS is more preferable than Indi-GPAS when the ratio between the populations that can be characterised with Indi-seq for every population that can be characterised with Pool-seq
    ### is at most 0.40 (10 Indi-seq : 25 Pool-seq) for TPR, and
    ### is at most 0.67 (67 Indi-seq : 100 Pool-seq) for FPR.

} ### end of loop testing different landscape scenarios
### reset dat:
dat = dat_old

# ### Power analysis to test what is the minimum allelic effect we could ever hope to detect with Indi-GPAS and Pool-GPAS
# ### Based on Appendix A of Visscher et al, 2017 (10 Years of GWAS Discovery:Biology, Function, and Translation)
# ### Assumes the absolute values of the QTL effects are chi-square distributed
# n = list(indi=384, pool=5)      ### sample sizes
# m = 10000                       ### number of loci
# h2 = seq(0, 1, length=1000)      ### heritability of the QTL
# beta = c(0.0, 0.0001, 0.001, 0.01, 0.1, 0.5, 0.75, 1.00)    ### effect size of the QTL
# Phet = c(0.0, 1/unlist(n), 0.50, 1-(1/unlist(rev(n))), 1.00)            ### frequency of heterozygotes assuming HWE
# pval = 0.05/m          ### Bonferroni threshold
# ### Power as a function of QTL heritability, number of samplesm and p-value
# power_h2 <- function(h2, n, pval) {
#     ### minimum effet of the QTL ot the quantile value on the chi-square distribution at df=1 corresponding to the p-value
#     q = qchisq(pval, df=1, lower.tail=F)
#     ### initialise the vector of powers
#     power = matrix(NA, ncol=1, nrow=length(h2))
#     ### calculate the power for each QTL heritability
#     for (i in 1:length(h2)) {
#         ### the non-centrality parameter of the chi-square distribution
#         ncp = (n * h2[i]) / (1 - h2[i])
#         # calculate the power as the area under the chi-square distribution from the minimum effect of the QTL based on the p-value up to infinity
#         if (is.finite(ncp) && ncp>=0) { ### but first test if ncp is within the domain of the chi-square function
#            power[i, 1] = pchisq(q=q, df=1, lower.tail=F, ncp=ncp)
#         }
#     }
#     return(power)
# }
# ### Power as a function of QTL effect, heterozygote frequency, number of samplesm and p-value
# power_bph <- function(beta, Phet, n, pval) {
#     ### minimum effet of the QTL ot the quantile value on the chi-square distribution at df=1 corresponding to the p-value
#     q = qchisq(pval, df=1, lower.tail=F)
#     ### initialise the matrix of powers: length(beta) rows x length(Phet) columns
#     POWER = matrix(NA, nrow=length(beta), ncol=length(Phet))
#     for (i in 1:length(beta)) {
#         for (j in 1:length(Phet)) {
#             # i = 2; j=2
#             ### calculate the variance explained by the QTL
#             h2 = Phet[j] * (beta[i]^2)
#             ### the non-centrality parameter of the chi-square distribution
#             ncp = (n * h2) / (1 - h2)
#             # calculate the power as the area under the chi-square distribution from the minimum effect of the QTL based on the p-value up to infinity
#             if (is.finite(ncp) && ncp>=0) {
#                 POWER[i,j] = pchisq(q, df=1, lower.tail=F, ncp=ncp)
#             }
#         }
#     }
#     return(POWER)
# }

# ### Calculate power as a function of QTL heritability
# p1_indi = power_h2(h2=h2, n=n$indi, pval=pval)
# p1_pool = power_h2(h2=h2, n=n$pool, pval=pval)
# # ### Calculate power as a function of the QTL effec and heterozygote frequency
# # P2_indi = power_bph(beta, Phet, n=n$indi, pval=pval)
# # P2_pool = power_bph(beta, Phet, n=n$pool, pval=pval)

# ### QTL heritability (assuminbg HWE across loci and p=q=0.5)
# p = 0.5; q = 0.5
# vec_nQTL = c(10, 50, 100)
# for (i in 1:length(vec_nQTL)){
#     b = round(rchisq(n=vec_nQTL[i], df=1),2)
#     var_b = var(b)
# }


# ### Plot powers
# ### Find h2 at which power is at 50%
# h2_0.5_indi = min(h2[p1_indi>=0.5], na.rm=TRUE)
# h2_0.5_pool = min(h2[p1_pool>=0.5], na.rm=TRUE)
# svg("1.6_Power_calculations.svg", width=10, height=10)
#     par(mfrow=c(2,1))
#     plot(x=h2, y=p1_indi, xlab="Variance Explained by the QTL", ylab="Power", type="l", lwd=2, las=2, main="Indi-GPAS")
#     abline(h=0.5, col="red")
#     abline(v=h2_0.5_indi, col="red")
#     grid()
#     legend("bottomright", legend=c(paste0("n = ", n$indi), paste0("Bonferroni threshold = ", pval), paste0("h2_QTL at 50% power = ", round(h2_0.5_indi,2))))
#     plot(x=h2, y=p1_pool, xlab="Variance Explained by the QTL", ylab="Power", type="l", lwd=2, las=2, main="Pool-GPAS (Improper Calculations!!!)\n[assumes allele count data and not allele frequencies]")
#     abline(h=0.5, col="red")
#     abline(v=h2_0.5_pool, col="red")
#     grid()
#     legend("bottomright", legend=c(paste0("n = ", n$pool), paste0("Bonferroni threshold = ", pval), paste0("h2_QTL at 50% power = ", round(h2_0.5_pool,2))))
# dev.off()

##################
### QUESTION 2 ### Under which population-specific circumstances should we use Pool-seq over Indi-seq?
################## Note that Indi-seq: 384 individuals/population, while Pool-seq: 5 pools/population.
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
### Dataset, metrics and explanatory variables ###
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
### We will use the intra-population dataset and only the two best GPAS models per genotryping scheme.
### This is to avoid the confounding effefs of multiple populations and the different models in each genotyping scheme.
dat_intra = droplevels(dat[dat$NPOP_MERGED==1, ])
### We will be using both GWAS and GP performance metrics, specifically AUC and RMSE.
### For the explanatory variables, we will use the all variables except for the number of populations samples and GPAS model, since we will only be comparing the best Indi- and Pool-GPAS models.
question_2_metric_names = c("AUC_MEAN", "RMSE_MEAN")
question_2_metric_labels = c("Area under the ROC plot (AUC)",
                             "Root Mean Square Error (RMSE)")
question_2_explain_names = c("NQTL", "SELECTION_INTENSITY", "MIGRATION_RATE", "QTL_GRADIENT", "GENOTYPING_SCHEME")
question_2_explain_labels = c("Number of Simulated QTL",
                              "Selection Intensity",
                              "Migration Rate",
                              "QTL Diffusion Gradient",
                              "Genotyping Method")
### Visualise the distributions and relationship of the two metrics we will be using.
png("2.1_HISTOGRAMS-metrics-subset-intra-2GPAS_models.png", width=500, height=900)
par(mfrow=c(3,1))
hist(dat_intra$AUC_MEAN, xlab="AUC", main="")
hist(dat_intra$RMSE_MEAN, xlab="RMSE", main="")
plot(dat_intra$AUC_MEAN, dat_intra$RMSE_MEAN, xlab="AUC", ylab="RMSE", pch=20, col=rgb(0.2,0.2,0.2,alpha=0.2))
dev.off()
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
### Mean compariasons under different levels of each explanatory variable and linear mixed modelling ###
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
### Violin plotting and simple mean comparisons between Indi-GPAS and Pool-GPAS under the different levels of the lanscape-specific variables.
vec_metric = c()
vec_explain = c()
vec_explain_level = c()
vec_difference = c()
vec_significance = c()
svg("2.2_VIOLINPLOTS-INTRAPOP-METRICS-ACROSS-XLEVELS.svg", width=30, height=8)
par(mfrow=c(2,12))
for (metric in c("AUC_MEAN", "RMSE_MEAN")){
    # metric = "AUC_MEAN"
    for (explain in names_explain[1:4]){
        # explain = names_explain[1:4][1]
        for (explain_level in eval(parse(text=paste0("unique(dat_intra$", explain, ")")))){
            # explain_level = eval(parse(text=paste0("unique(dat_intra$", explain, ")")))[1]
            sub_dat = eval(parse(text=paste0("droplevels(dat_intra[dat_intra$", explain, "==explain_level, ])")))
            ### violin plot
            out = eval(parse(text=paste0("violinplotter(", metric, "~ GENOTYPING_SCHEME, data=sub_dat, TITLE=paste0(explain, ':\n', explain_level))")))
            ### means and parsing HSD output
            agg = eval(parse(text=paste0("aggregate(", metric, " ~ GENOTYPING_SCHEME, data=sub_dat, FUN=mean)")))
            difference = agg[1,2] - agg[2,2]
            if (length(unique(out[[1]]$HSD_out.GROUPING))==2){
                significance = TRUE
            } else {
                significance = FALSE
            }
            vec_metric = c(vec_metric, metric)
            vec_explain = c(vec_explain, explain)
            vec_explain_level = c(vec_explain_level, explain_level)
            vec_difference = c(vec_difference, difference)
            vec_significance = c(vec_significance, significance)
        }
    }
}
dev.off()
### Barplotting the significant and non-significant differences
BARPLOT = data.frame(METRIC=vec_metric, EXPLANATORY_VAR=vec_explain, LEVEL=as.numeric(vec_explain_level), DIFF=vec_difference, SIG=vec_significance)
BARPLOT = BARPLOT[order(BARPLOT$LEVEL, decreasing=FALSE), ]
BARPLOT = BARPLOT[order(BARPLOT$EXPLANATORY_VAR, decreasing=FALSE), ]
for (metric in unique(BARPLOT$METRIC)){
    # metric = unique(BARPLOT$METRIC)[1]
    sub_BARPLOT = droplevels(BARPLOT[(BARPLOT$METRIC==metric)&(BARPLOT$EXPLANATORY_VAR!="QTL_GRADIENT"), ])
    if (metric == "AUC_MEAN"){
        y.lim=c(0, max(sub_BARPLOT$DIFF)+(sd(sub_BARPLOT$DIFF)/3))
        text_adj = +sd(sub_BARPLOT$DIFF)/3
    } else {
        y.lim = c(min(sub_BARPLOT$DIFF)+(sd(sub_BARPLOT$DIFF)/3), max(sub_BARPLOT$DIFF)+(sd(sub_BARPLOT$DIFF)/3))
        text_adj = -sd(sub_BARPLOT$DIFF)/3
    }
    svg(paste0("2.3_BARPLOTS_DIFF_Indi_Pool-", metric, ".svg"), width=10, height=6)
    par(mar=c(7,7,5,2))
    bp = barplot(sub_BARPLOT$DIFF,
                 ylim=y.lim,
                 names=paste0(substr(sub_BARPLOT$EXPLANATORY_VAR,1,4), "-", sub_BARPLOT$LEVEL),
                 ylab=paste0(metric, " Difference\n(Indi-GPAS - Pool-GPAS)"),
                 bord=FALSE,
                 col=c("#fcae91", "#fb6a4a", "#cb181d",
                       "#a1dab4", "#41b6c4", "#225ea8",
                       "#bae4b3", "#74c476", "#238b45"),
                 las=2)
    sig_lab = rep("*", 9)
    sig_lab[sub_BARPLOT$SIG==FALSE] = "ns"
    text(bp, sub_BARPLOT$DIFF+text_adj, lab=sig_lab)
    dev.off()
}
### Results of simple mean comparisons:
### Indi-GPAS performs better than Pool-GPAS across the different landscape variables except for when selection intensity and migration rates are high.
### At high selection intensities (SELECTION_INTENSITY={0.90, 0.95}) and high gene flow (migration rate = 0.01), Indi-GPAS and Pool-GPAS performances were not significantly different.
# ### We attempt to confirm these findings using linear mixed modelling.
# ### Instantiate these vectors for plotting the line plots below.
# vec_metric = c()
# vec_explain = c()
# vec_explain_level = c()
# vec_difference = c()
# for (metric in question_2_metric_names){
#     # metric = "RMSE_MEAN"
#     # mod_LMM = eval(parse(text=paste0("lmer(", metric, " ~ LOG_NQTL + LOG_SELE + LOG_MIGR + QTL_GRADIENT + GENOTYPING_SCHEME +
#     mod_LMM = eval(parse(text=paste0("lmer(", metric, " ~ NQTL + SELECTION_INTENSITY + LOG_MIGR + QTL_GRADIENT + GENOTYPING_SCHEME +
#                                     (0 + GENOTYPING_SCHEME|NQTL) +
#                                     (0 + GENOTYPING_SCHEME|SELECTION_INTENSITY) +
#                                     (0 + GENOTYPING_SCHEME|MIGRATION_RATE), data=dat_intra, 
#                                     control=lmerControl(check.conv.singular = .makeCC(action = 'ignore',  tol = 1e-4)))")))
#     # fixef(mod_LMM)
#     # ranef(mod_LMM)
#     FIXEF = data.frame(fixef(mod_LMM))
#     FIXEF = cbind(rownames(FIXEF), FIXEF); rownames(FIXEF) = NULL
#     colnames(FIXEF) = c("Variable", paste0("Fixed Effects on ", metric))
#     write.table(FIXEF, file=paste0("2.4_FIXEF-intrapop-", metric, ".csv"), sep=",", quote=TRUE, row.names=FALSE, col.names=TRUE)
#     ### Random slopes, TPR/NPOP under the different levels of the explanatory variables
#     for (explain in names(ranef(mod_LMM))){
#         # explain = names(ranef(mod_LMM))[1]
#         RANEF = eval(parse(text=paste0("data.frame(ranef(mod_LMM)$", explain, ")")))
#         RANEF = cbind(rownames(RANEF), RANEF); rownames(RANEF) = NULL
#         colnames(RANEF) = c(explain, paste0("SLOPE_", metric, "_OVER_INDI-GPAS"), paste0("SLOPE_", metric, "_OVER_POOL-GPAS"))
#         write.table(RANEF, file=paste0("2.4_RANEF-intrapop-", metric, "-", explain, ".csv"), sep=",", quote=TRUE, row.names=FALSE, col.names=TRUE)
#     }
#     ### Populate the vectors in preparation for the line plots below.
#     for (explain in names(ranef(mod_LMM))){
#         # explain = names(ranef(mod_LMM))[1]
#         RANEF = eval(parse(text=paste0("data.frame(ranef(mod_LMM)$", explain, ")")))
#         vec_metric = c(vec_metric, rep(metric, nrow(RANEF)))
#         vec_explain = c(vec_explain, rep(explain, nrow(RANEF)))
#         vec_explain_level = c(vec_explain_level, rownames(RANEF))
#         vec_difference = c(vec_difference, RANEF[,1] - RANEF[,2])
#     }
# }
# ### plot line graphs of npop X metric where each line correspond to a factor level
# LINEPLOT = data.frame(METRIC=vec_metric, EXPLANATORY_VAR=vec_explain, LEVEL=vec_explain_level, DIFF=vec_difference)
# LINEPLOT$LEVEL = as.numeric(as.character(LINEPLOT$LEVEL))
# svg(paste0("2.5_RANDOM_SLOPE_DIFF.svg"), width=12, height=7)
# par(mfrow=c(2,3))
# for (metric in unique(LINEPLOT$METRIC)){
#     # metric = unique(LINEPLOT$METRIC)[1]
#     for (explain in unique(LINEPLOT$EXPLANATORY_VAR)){
#         # explain = unique(LINEPLOT$EXPLANATORY_VAR)[1]
#         sub_LINEPLOT = droplevels(LINEPLOT[(LINEPLOT$METRIC==metric)&(LINEPLOT$EXPLANATORY_VAR==explain), ])
#         y.lim = c(min(LINEPLOT$DIFF[(LINEPLOT$METRIC==metric)]), max(LINEPLOT$DIFF[(LINEPLOT$METRIC==metric)]))
#         par(mar=c(5,6,2,2))
#         plot(x=sub_LINEPLOT$LEVEL, y=sub_LINEPLOT$DIFF,
#              ylim=y.lim,
#              type="b", pch=19, lty=2,
#              xlab=explain, ylab=paste0(metric, " Random Slope Difference\n(Indi-GPAS - Pool-GPAS)"),
#              las=2)
#         abline(h=0.0, lty=2, col="#fb6a4a", lwd=2)
#         grid(lty=2,col="gray",lwd=0.5)
#     }
# }
# dev.off()
# ### Results:
# ### For AUC (GWAS performance): Pool-GPAS are better than Indi-GPAS when selection intensity, migration rate and NQTL are high, (i.e. random slopes are higher for Pool-GPAS).
# ### For RMSE (GP performance): Pool-GPAS are better than Indi-GPAS when selection intensity, migration rate and NQTL are high, (i.e. random slopes are lower for Pool-GPAS).
# ### Consolidating the results:
# ### When we expect the number of QTL, selection intensity, and gene flow to be high in the landscape of interest,
# ### then we can use Pool-GPAS instead of Indi-GPAS especially when budget is tight.
# ### Just be sure to keep the expectations of TPR and FPR aligned with the previous results of the weighted optimisation: "WEIGHTED_OPTIMISATION-NPOP-TPR-FPR.png"

##################
### QUESTION 3 ### Which populations best represent the landscape under different circumstances?
##################
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
### Dataset, metrics and explanatory variables ###
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
### The same with question 2, we will use the intra-population dataset and only the two best GPAS models per genotryping scheme.
dat_intra = droplevels(dat[dat$NPOP_MERGED==1, ])
### We will also be using the same metrics and explanatory variables as question 2.
question_3_metric_names = c("AUC_MEAN", "RMSE_MEAN")
question_3_metric_labels = c("Area under the ROC plot (AUC)",
                             "Root Mean Square Error (RMSE)")
question_3_explain_names = c("NQTL", "SELECTION_INTENSITY", "MIGRATION_RATE", "QTL_GRADIENT", "GENOTYPING_SCHEME")
question_3_explain_labels = c("Number of Simulated QTL",
                              "Selection Intensity",
                              "Migration Rate",
                              "QTL Diffusion Gradient",
                              "Genotyping Method")
### However, we will be dividing our intra-population dataset by QTL diffusion gradients.
### This is to allow us to group populations by rows as a function of these gradient.
### This will help in interpretability of the results for applications in real-world sampling endeavors.
dat_GRAD0 = droplevels(dat_intra[dat_intra$QTL_GRADIENT==0, ])
dat_GRAD1 = droplevels(dat_intra[dat_intra$QTL_GRADIENT==1, ])
dat_GRAD2 = droplevels(dat_intra[dat_intra$QTL_GRADIENT==2, ])
### Is row-based grouping sensible? Let's look at the metric distributions per population per gradient dataset
### Plot performance distributions per population using each QTl diffusion gradient dataset.
names_dataset = c("dat_GRAD0", "dat_GRAD1", "dat_GRAD2")
for (metric in question_3_metric_names){
    # metric = question_3_metric_names[1]
    svg(paste0("3.1_VIOLINPLOTS-POP_GRAD-", metric, ".svg"), width=9, height=12)
    par(mfrow=c(3,1))
    for (dataset in names_dataset){
        # dataset = "dat_GRAD1"
        eval(parse(text=paste0("violinplotter(", metric, " ~ POP_NAMES, data=", dataset, ", TITLE='", dataset, "', HSDX=FALSE)")))
    }
    dev.off()
}
### For each QTL diffusion gradient dataset, we will map each population to belong to 1 of 10 rows perpedicular to the gradient.
for (dataset in names_dataset){
    eval(parse(text=paste0(dataset, "$ROW_GROUP = rep(1, times = nrow(", dataset, "))")))
    vec_row_group = c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10")
    vec_start = seq(from=1, to=100, by=10)
    vec_end = seq(from=10, to=100, by=10)
    for (i in 1:length(vec_row_group)){
        row_group = vec_row_group[i]
        start = vec_start[i]
        end = vec_end[i]
        eval(parse(text=paste0("vec_pop_num = as.numeric(sub('POP_', '', as.character(", dataset, "$POP_NAMES)))")))
        eval(parse(text=paste0(dataset, "$ROW_GROUP[(vec_pop_num >= start) & (vec_pop_num <= end)] = row_group")))
    }
    eval(parse(text=paste0(dataset, "$ROW_GROUP = as.factor(", dataset, "$ROW_GROUP)")))
}
### Now, plot the performance distributions per row grouping this time.
for (metric in question_3_metric_names){
    # metric = question_3_metric_names[1]
    svg(paste0("3.1_VIOLINPLOTS-ROW_GRAD-", metric, ".svg"), width=11, height=4)
    par(mfrow=c(1,3))
    for (dataset in names_dataset){
        # dataset = "dat_GRAD1"
        eval(parse(text=paste0("violinplotter(", metric, " ~ ROW_GROUP, data=", dataset, ", TITLE='", dataset, "', HSDX=TRUE)")))
    }
    dev.off()
}
### Results:
### Yes, grouping by row is informative.
### Notice the progresively decreasing performances from POP001 to POP100 in the GRAD1 dataset, as well as
### the slight dip in performances in the middle populations in the GRAD2 dataset.
### And finally the constant performances for GRAD0.
### These patterns reflect the QTL diffusion gradients (see Appendix 1 for details).

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
### Elucidate the relationships between GPAS performance and the relative position or identity of populations in the landscape ###
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
### From the violiplots and mean comparisons above, we can see that the regions (i.e. rows) from which the QTL originate have higher GPAS performances.
### We summarize these results using a full factorial linear model for mean comparisons and to extract the linear model-derived coefficients.
### We will re-use the function defined during the identification of the best GPAS models (before question 1).
formula_GWAS = AUC_MEAN  ~ LOG_NQTL*LOG_SELE*LOG_MIGR*GENOTYPING_SCHEME*ROW_GROUP
formula_GP   = RMSE_MEAN ~ LOG_NQTL*LOG_SELE*LOG_MIGR*GENOTYPING_SCHEME*ROW_GROUP
for (dataset in names_dataset){
    # dataset = names_dataset[2]
    MEANS_COEF_HSD_GWAS = func_model_Ftest_HSD_rank(formula=formula_GWAS, var_name="ROW_GROUP", data=eval(parse(text=dataset)), alpha=0.05)
    MEANS_COEF_HSD_GP   = func_model_Ftest_HSD_rank(formula=formula_GP,   var_name="ROW_GROUP", data=eval(parse(text=dataset)), alpha=0.05)
    write.table(MEANS_COEF_HSD_GWAS, file=paste0("3.2_MEANS_COEF_HSD_GWAS-", dataset, ".csv"), sep=",", row.names=FALSE, col.names=TRUE, quote=FALSE)
    write.table(MEANS_COEF_HSD_GP,   file=paste0("3.2_MEANS_COEF_HSD_GP-",   dataset, ".csv"), sep=",", row.names=FALSE, col.names=TRUE, quote=FALSE)
}
### Results:
### Regions (i.e. rows) from which the QTL originated from have higher GPAS performances.
### Simply put, collect populations from regions with genetic variability for the trait of interest.
### Stress on genetic variabily, since it should be differentiated from phenotypic variability which may just be caused by environmental variability.
### We want to sample from regions where we expect to have had history of selection while having genetic variability for the trait.
### How is this finding affected by the other explanatory variables?
### To answer this, we plotted violinplots and fitted linear mixed models.

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
### Mean compariasons under different levels of each explanatory variable and linear mixed modelling ###
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
### Plot distributions, perform mean comparisons and pseudo-regressions of metrics across the row groupings.
for (dataset in names_dataset){
    # dataset = "dat_GRAD1"
    svg(paste0("3.3_VIOLINPLOTS-ROW_GRAD-METRICS-ACROSS-XLEVELS-", dataset, ".svg"), width=30, height=8)
    par(mfrow=c(2,11))
    for (metric in question_3_metric_names){
        # metric = "AUC_MEAN"
        for (explain in question_3_explain_names[!grepl("QTL_GRADIENT", question_3_explain_names)]){
            # explain = question_3_explain_names[!grepl("QTL_GRADIENT", question_3_explain_names)][1]
            for (explain_level in eval(parse(text=paste0("unique(", dataset, "$", explain, ")")))){
                # explain_level = eval(parse(text=paste0("unique(", dataset, "$", explain, ")")))[1]
                sub_dat = eval(parse(text=paste0("droplevels(", dataset, "[", dataset, "$", explain, "==explain_level, ])")))
                eval(parse(text=paste0("violinplotter(", metric, "~ ROW_GROUP, data=sub_dat, TITLE=paste0(explain, ':\n', explain_level), REGRESSX=TRUE)")))
            }
        }
    }
    dev.off()
}
### Results:
### Similar with previous findings across the other explanatory variables:
### -   Generally, in the absence of a gradient [GRAD=0], no significant row seems to be favourable to sample,
###     except for a bit of favourability for populations in the middle.
### -   Generally, in the presence of a single diffusion front [GRAD=1],
###     in terms of GWAS performance, the populations from where the QTL originated from are the most favourable to be sampled, however
###     in terms of GP performance, the populations in the middle of the ladnscape seem to be more preferable.
### -   Generally, in the presence of two diffusion fronts from opposite sides of the landscape [GRAD=2],
###     in terms of both GWAS and GP performances, the populations from both sides from where the QTL originated are the most favorable to be sampled.
### Specifically in terms of GWAS performance:
### (1) NQTL: sampling farther away from the diffusion fronts becomes less important as the number of QTL increases (i.e. absolute value of slopes decrease).
### (2) SELECTION_INTENSITY: sampling farther away from the diffusion fronts becomes less important as selection intensity decreases (i.e. absolute value of slopes decrease).
### (3) MIGRATION_RATE: sampling farther away from the diffusion fronts becomes less important as migration rate increases (i.e. absolute value of slopes decrease).
### (4) GENOTYPING_SCHEME: sampling farther away from the diffusion fronts is more important for Pool-GPAS than Indi-GPAS.
### Finally, in terms of GP performance:
### The middle region is favourable under a single diffusion front (GRAD=1), but unfavourable under 2 diffusion fronts (GRAD=2).
### Under different levels of the other explanatory variables, the magnitudes of these favourable/ or unfavourableness are as follows:
### (1) NQTL: decreases as the number of QTL increases
### (2) SELECTION_INTENSITY: increases as selection intensity increases
### (3) MIGRATION_RATE: decreases as migration rate increases
### (4) GENOTYPING_SCHEME: more important for Pool-GPAS than Indi-GPAS

### We will try to confirm these finding with linear mixed models and polynomial (second degree) models based on the relationships observed from the violinplots and mean comparisons above.
for (dataset in names_dataset){
    # dataset = names_dataset[2]
    for (metric in question_3_metric_names){
        # metric = "AUC_MEAN"
        ### convert ROW_GROUP into numeric
        eval(parse(text=paste0(dataset, "$ROW_GROUP = as.numeric(as.character(", dataset, "$ROW_GROUP))")))
        ###@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        ### Mixed models without second degree polynomial random effects of row grouping.
        mod_LMM = eval(parse(text=paste0("lmer(", metric, " ~ LOG_NQTL + LOG_SELE + LOG_MIGR + GENOTYPING_SCHEME +
                                    (0 + ROW_GROUP|NQTL) +
                                    (0 + ROW_GROUP|SELECTION_INTENSITY) +
                                    (0 + ROW_GROUP|MIGRATION_RATE) +
                                    (0 + ROW_GROUP|GENOTYPING_SCHEME), data=", dataset, ",
                                    control=lmerControl(check.conv.singular=.makeCC(action='ignore', tol=1e-4)))")))
        # summary(rePCA(mod_LMM))
        # fixef(mod_LMM)
        # ranef(mod_LMM)
        FIXEF = data.frame(fixef(mod_LMM))
        FIXEF = cbind(rownames(FIXEF), FIXEF); rownames(FIXEF) = NULL
        colnames(FIXEF) = c("Variable", paste0("Fixed Effects on ", metric))
        write.table(FIXEF, file=paste0("3.4_FIXEF-ROW_GRAD-", metric, "-", dataset, "-NOPOLY.csv"), sep=",", quote=TRUE, row.names=FALSE, col.names=TRUE)
        ### Random slopes, TPR/NPOP under the different levels of the explanatory variables
        ### and plot the slopes (i.e random effects nested within each explanatory factor)
        svg(paste0("3.4_LINEAR_FIT-", dataset, "-", metric, ".svg"), width=12, height=3)
        par(mfrow=c(1,4))
        for (explain in names(ranef(mod_LMM))){
            # explain = names(ranef(mod_LMM))[1]
            RANEF = eval(parse(text=paste0("data.frame(ranef(mod_LMM)$", explain, ")")))
            RANEF = cbind(rownames(RANEF), RANEF); rownames(RANEF) = NULL
            colnames(RANEF) = c(explain, paste0("SLOPE_", metric, "_OVER_ROW_GROUP"))
            write.table(RANEF, file=paste0("3.4_RANEF-ROW_GRAD-", metric, "-", dataset, "-", explain, "-NOPOLY.csv"), sep=",", quote=TRUE, row.names=FALSE, col.names=TRUE)
            plot(x=c(1,10), y=c(min(c(0, min(RANEF[,2]))),max(RANEF[,2])), type="n", xlab="Row Group", ylab=metric, main=explain, las=2)
            grid()
            for (i in 1:nrow(RANEF)){
                # i = 1
                lines(x=c(0,10), y=c(0,RANEF[i,2]), lty=1, lwd=2, col=c("#d73027", "#fdae61", "#31a354")[i])
            }
            legend("bottomleft", legend=RANEF[,1], lty=1, lwd=2, col=c("#d73027", "#fdae61", "#31a354"))
        }
        dev.off()
        ### Test if the fit was singular using the PCA of the random covariance matrix. (NAs and zero standard deviations are evidence of singularities. This implies that the random components are too complex, and ommiting one to several or even all the random effects may be better.)
        PCA = matrix(unlist(summary(rePCA(mod_LMM))), ncol=7, byrow=TRUE)
        colnames(PCA) = c("Standard Deviation", "Rotation", "Center", "Scale", "Importance1", "Importance2", "Importance3")
        rownames(PCA) = names(summary(rePCA(mod_LMM)))
        write.table(RANEF, file=paste0("3.4_RANEF-ROW_GRAD-", metric, "-", dataset, "-PCA-NOPOLY.csv"), sep=",", quote=TRUE, row.names=FALSE, col.names=TRUE)
        ###@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        ### Mixed models with second degree polynomial random effects of row grouping. (Note that the random effects of poly(ROW_GROUP,2) indicate the direction and magnitude of the trough or crest of the relationship)
        mod_POLY = eval(parse(text=paste0("lmer(", metric, " ~ LOG_NQTL + LOG_SELE + LOG_MIGR + GENOTYPING_SCHEME +
                                    (poly(ROW_GROUP,2, raw=TRUE)|NQTL) +
                                    (poly(ROW_GROUP,2, raw=TRUE)|SELECTION_INTENSITY) +
                                    (poly(ROW_GROUP,2, raw=TRUE)|MIGRATION_RATE) +
                                    (poly(ROW_GROUP,2, raw=TRUE)|GENOTYPING_SCHEME), data=", dataset, ",
                                    control=lmerControl(check.conv.singular=.makeCC(action='ignore', tol=1e-4)))")))
        # summary(rePCA(mod_POLY))
        # fixef(mod_POLY)
        # ranef(mod_POLY)
        FIXEF = data.frame(fixef(mod_POLY))
        FIXEF = cbind(rownames(FIXEF), FIXEF); rownames(FIXEF) = NULL
        colnames(FIXEF) = c("Variable", paste0("Fixed Effects on ", metric))
        write.table(FIXEF, file=paste0("3.5_FIXEF-ROW_GRAD-", metric, "-", dataset, "-POLY.csv"), sep=",", quote=TRUE, row.names=FALSE, col.names=TRUE)
        ### Solve fo the vertex of the parabola (y = ax^2 - 2ahx + ah^2 + k; where vertex=(h,k)) and plot
        svg(paste0("3.5_POLY_FIT-", dataset, "-", metric, ".svg"), width=12, height=3)
        par(mfrow=c(1,4))
        for (explain in names(ranef(mod_POLY))){
            # explain = names(ranef(mod_POLY))[3]
            RANEF = eval(parse(text=paste0("data.frame(ranef(mod_POLY)$", explain, ")")))
            RANEF = cbind(rownames(RANEF), RANEF); rownames(RANEF) = NULL
            colnames(RANEF) = c(explain, "c", "b", "a")
            a = RANEF$a
            b = RANEF$b
            c = RANEF$c ### intercept
            # c = 0
            RANEF$h = -b/(2*a)
            RANEF$k = c - (a*RANEF$h^2)
            write.table(RANEF, file=paste0("3.5_RANEF-ROW_GRAD-", metric, "-", dataset, "-", explain, "-POLY.csv"), sep=",", quote=TRUE, row.names=FALSE, col.names=TRUE)
            ### plot parabolas/parabolae?
            x = seq(1, 10, length=100)
            # y_extreme = max(c(abs(RANEF$k), 1))
            y_extreme = max(abs(a*(rbind(x,x,x) - RANEF$h)^2 + RANEF$k))
            plot(x=c(1,10), y=c(-y_extreme, +y_extreme), type="n", xlab="Row Group", ylab=metric, main=explain, las=2)
            grid()
            for (i in 1:nrow(RANEF)){
                # i = 1
                y = a[i]*(x - RANEF$h[i])^2 + RANEF$k[i]
                lines(x, y, lwd=2, col=c("#d73027", "#fdae61", "#31a354")[i])
            }
            legend("top", legend=as.character(RANEF[,1]), lty=1, lwd=2, col=c("#d73027", "#fdae61", "#31a354"))
        }
        dev.off()
        ### Test if the fit was singular using the PCA of the random covariance matrix. (NAs and zero standard deviations are evidence of singularities. This implies that the random components are too complex, and ommiting one to several or even all the random effects may be better.)
        PCA = matrix(unlist(summary(rePCA(mod_POLY))), ncol=23, byrow=TRUE)
        colnames(PCA) = c(paste0("Standard Deviation", 1:3), paste0("Rotation", 1:9), "Center", "Scale", paste0("Importance", 1:9))
        rownames(PCA) = names(summary(rePCA(mod_POLY)))
        write.table(RANEF, file=paste0("3.5_RANEF-ROW_GRAD-", metric, "-", dataset, "-PCA-POLY.csv"), sep=",", quote=TRUE, row.names=FALSE, col.names=TRUE)
        ###@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    }
}
### Results:
### Observations from the violin plots are confirmed.
### Note that polynomial and non-polynomial fits were chosen to confirm the violin plot results according to the shapes of the relationships based on the same plots.

### Assessing if these favourable regions in the landscape coincide with regions of high variability
var_summstats = read.csv(fname_var_summstats)
### compute the fraction of the causal loci captured per population
var_summstats$FRAC_QTL = var_summstats$NQTL_CAPTURED/var_summstats$NQTL
### express selection intensity in an intuitive way
var_summstats$SEL= 1 - var_summstats$SEL
### divide into the 3 QTL diffusion gradients
varSS_GRAD0 = droplevels(var_summstats[var_summstats$GRAD==0, ])
varSS_GRAD1 = droplevels(var_summstats[var_summstats$GRAD==1, ])
varSS_GRAD2 = droplevels(var_summstats[var_summstats$GRAD==2, ])
### Is row-based grouping sensible? Let's look at the metric distributions per population per gradient dataset
### Plot performance distributions per population using each QTl diffusion gradient dataset.
names_dataset = c("varSS_GRAD0", "varSS_GRAD1", "varSS_GRAD2")
### For each QTL diffusion gradient dataset, we will map each population to belong to 1 of 10 rows perpedicular to the gradient.
for (dataset in c("var_summstats", names_dataset)){
    eval(parse(text=paste0(dataset, "$ROW_GROUP = rep(1, times = nrow(", dataset, "))")))
    vec_row_group = c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10")
    vec_start = seq(from=1, to=100, by=10)
    vec_end = seq(from=10, to=100, by=10)
    for (i in 1:length(vec_row_group)){
        row_group = vec_row_group[i]
        start = vec_start[i]
        end = vec_end[i]
        eval(parse(text=paste0("vec_pop_num = as.numeric(sub('POP_', '', as.character(", dataset, "$POP)))")))
        eval(parse(text=paste0(dataset, "$ROW_GROUP[(vec_pop_num >= start) & (vec_pop_num <= end)] = row_group")))
    }
    eval(parse(text=paste0(dataset, "$ROW_GROUP = as.factor(", dataset, "$ROW_GROUP)")))
}
### Plot distributions per row
for (metric in c("PHENO_MEAN", "PHENO_VAR", "NQTL_CAPTURED", "FRAC_QTL", "EXPECTED_HET", "OBSERVED_HET", "QTL_2pq", "CAUSAL_ALLELE_FREQ")){
    # metric = "PHENO_MEAN"
    svg(paste0("3.6_VIOLINPLOTS-VARSS-GRAD-", metric, ".svg"), width=11, height=4)
    par(mfrow=c(1,3))
    for (dataset in names_dataset){
        # dataset = "dat_GRAD1"
        eval(parse(text=paste0("violinplotter(", metric, " ~ ROW_GROUP, data=", dataset, ", TITLE='", dataset, "', HSDX=TRUE)")))
    }
    dev.off()
}


### Assessing changes in QTL capture in each population across different levels of the landscape parameters
# svg(paste0("0.7_VIOLINPLOTS-QTL_CAPTURED_vs_NQTL_SEL_MIG.svg"), width=11, height=11)
# par(mfrow=c(3,3))
# violinplotter(FRAC_QTL ~ NQTL, data=var_summstats)
# violinplotter(FRAC_QTL ~ SEL, data=var_summstats)
# violinplotter(FRAC_QTL ~ MIG, data=var_summstats)
# violinplotter(NQTL_CAPTURED ~ NQTL, data=var_summstats)
# violinplotter(NQTL_CAPTURED ~ SEL, data=var_summstats)
# violinplotter(NQTL_CAPTURED ~ MIG, data=var_summstats)
# violinplotter(CAUSAL_ALLELE_FREQ ~ NQTL, data=var_summstats)
# violinplotter(CAUSAL_ALLELE_FREQ ~ SEL, data=var_summstats)
# violinplotter(CAUSAL_ALLELE_FREQ ~ MIG, data=var_summstats)
# dev.off()
for (i in 1:2){
    # dataset = "varSS_GRAD1"
    svg(paste0("3.7_VIOLINPLOTS-ROW_GRAD-VAR-ACROSS-XLEVELS-", c("GRAD1", "GRAD2")[i], ".svg"), width=30, height=10)
    par(mfrow=c(3,9))
    ### AUC & RMSE
    dataset = c("dat_GRAD1", "dat_GRAD2")[i]
    for (metric in question_3_metric_names){
        # metric = "AUC_MEAN"
        for (explain in c("NQTL", "SELECTION_INTENSITY", "MIGRATION_RATE")){
            # explain = question_3_explain_names[!grepl("QTL_GRADIENT", question_3_explain_names)][1]
            for (explain_level in eval(parse(text=paste0("unique(", dataset, "$", explain, ")")))){
                # explain_level = eval(parse(text=paste0("unique(", dataset, "$", explain, ")")))[1]
                sub_dat = eval(parse(text=paste0("droplevels(", dataset, "[", dataset, "$", explain, "==explain_level, ])")))
                # eval(parse(text=paste0("violinplotter(", metric, "~ ROW_GROUP, data=sub_dat, TITLE=paste0(explain, ':\n', explain_level), REGRESSX=TRUE)")))
                eval(parse(text=paste0("violinplotter(", metric, "~ ROW_GROUP, data=sub_dat, TITLE=paste0(explain, ':\n', explain_level), REGRESSX=FALSE)")))
            }
        }
    }
    ### Causal allele freq
    dataset = c("varSS_GRAD1", "varSS_GRAD2")[i]
    for (metric in c("CAUSAL_ALLELE_FREQ")){
        # metric = "PHENO_MEAN"
        for (explain in c("NQTL", "SEL", "MIG")){
            # explain = question_3_explain_names[!grepl("QTL_GRADIENT", question_3_explain_names)][1]
            for (explain_level in eval(parse(text=paste0("unique(", dataset, "$", explain, ")")))){
                # explain_level = eval(parse(text=paste0("unique(", dataset, "$", explain, ")")))[2]
                sub_dat = eval(parse(text=paste0("droplevels(", dataset, "[", dataset, "$", explain, "==explain_level, ])")))
                # eval(parse(text=paste0("violinplotter(", metric, "~ ROW_GROUP, data=sub_dat, TITLE=paste0(explain, ':\n', explain_level), REGRESSX=TRUE)")))
                eval(parse(text=paste0("violinplotter(", metric, "~ ROW_GROUP, data=sub_dat, TITLE=paste0(explain, ':\n', explain_level), REGRESSX=FALSE)")))
            }
        }
    }
    dev.off()
}
### These trends in genetic diversity explain the GPAS performance trends in the 3.3 violin plots and the mixed modelling results

# ####################################################
# ### ADDITIONAL QUESTION OF ACROSS VS WITHIN GPAS ###
# #################################################### 2020-09-01
# library(violinplotter)

# # ### load merged streamlined dataset containing intra-(WITHIN) anf inter-(ACROSS) population datasets
# dat = readRDS("MERGED_GPAS_OUTPUT.rds")
# ### list explanatory variable names and labels
# dat$NQTL =                as.numeric(dat$NQTL)
# dat$SELECTION_INTENSITY = as.numeric(1.00 - dat$SELECTION_INTENSITY) ### reverted since the original variable is the distance of phenotypic values from the maximu possible phenotype value, i.e. lower values means higher selecetion intensities
# dat$MIGRATION_RATE =      as.numeric(dat$MIGRATION_RATE)
# dat$QTL_GRADIENT =        as.factor(dat$QTL_GRADIENT)
# dat$GENOTYPING_SCHEME =   as.factor(dat$GENOTYPING_SCHEME)
# dat$GPAS_MODEL =          as.factor(paste0(dat$MODEL, "-", dat$COVARIATE))
# dat$NPOP_TRAIN = unlist(lapply(strsplit(as.character(dat$TRAIN_POP), ":"), FUN=length))
# ### Subset the main dataset to include only GEMMA-STANDARDIZED and GWAlpha_SNPwise to represent Indi-GPAS and Pool-GPAS models, respectively
# dat = droplevels(dat[(dat$GPAS_MODEL=="GEMMA-STANDARDIZED") | (dat$GPAS_MODEL=="GWAlpha_SNPwise-NONE"), ])
# # dat = droplevels(dat[(dat$GPAS_MODEL=="GEMMA-STANDARDIZED") | (dat$GPAS_MODEL=="GWAlpha_MIXED_REML-HIVERT"), ])
# gc()
# saveRDS(dat, "MERGED_GPAS_OUTPUT_STREAMLINED.rds")
# dat = readRDS("MERGED_GPAS_OUTPUT_STREAMLINED.rds")

# aggregate(AUC ~ NPOP_TRAIN, data=dat, FUN=mean)
# aggregate(RMSE ~ NPOP_TRAIN, data=dat, FUN=mean)


# ### aggregate data across training sets
# agg_AUC = aggregate(AUC ~ NQTL + SELECTION_INTENSITY + MIGRATION_RATE + QTL_GRADIENT + GENOTYPING_SCHEME + TRAIN_POP + NPOP_TRAIN, data=dat, FUN=mean)
# agg_RMSE = aggregate(RMSE ~ NQTL + SELECTION_INTENSITY + MIGRATION_RATE + QTL_GRADIENT + GENOTYPING_SCHEME + TRAIN_POP + NPOP_TRAIN, data=dat, FUN=mean)
# ### remove NPOP_TRAIN == 5 to make the plots clearer
# agg_AUC = droplevels(agg_AUC[agg_AUC$NPOP_TRAIN != 5, ])
# agg_RMSE = droplevels(agg_RMSE[agg_RMSE$NPOP_TRAIN != 5, ])

# svg("test_intra_vs_inter.svg", width=12, height=12)
# # svg("test_intra_vs_inter-GEMAA_STD_vs_MIXED_HIVERT.svg", width=12, height=12)
# par(mfrow=c(2,2))
# violinplotter(AUC ~ NPOP_TRAIN, data=agg_AUC[agg_AUC$GENOTYPING_SCHEME=="INDIVIDUAL", ], TITLE="INDIVIDUAL")
# violinplotter(AUC ~ NPOP_TRAIN, data=agg_AUC[agg_AUC$GENOTYPING_SCHEME=="POOL", ], TITLE="POOL")
# violinplotter(RMSE ~ NPOP_TRAIN, data=agg_RMSE[agg_RMSE$GENOTYPING_SCHEME=="INDIVIDUAL", ], TITLE="INDIVIDUAL")
# violinplotter(RMSE ~ NPOP_TRAIN, data=agg_RMSE[agg_RMSE$GENOTYPING_SCHEME=="POOL", ], TITLE="POOL")
# dev.off()


