#########################################################################################################
### CALCULATE GPAS PERFORMANCE AS THE NUMBER OF INDEPENDENT GPAS EXPERIMENTS PER POPULATION INCREASES ###
#########################################################################################################
### Input:
### (1) GPAS_OUTPUT.csv (full path)
### Output:
### (1) GPAS_OUTPUT_NPOP_PERFORMANCE.csv

args = commandArgs(trailingOnly=TRUE)
# args = c("GPASim_1REP_10NQTL_100NBGS_0.0001MIG_0.50SEL_0.10BGSEL_1GRAD/OUTPUT/GPAS_OUTPUT.csv")
fname_input = args[1]
print(fname_input)
### load data
dat = read.csv(fname_input)
### restricting maximum RMSE to 1 for computational ease (since an RMSE of 1 simply says that the genomic prediction error maximized, i.e. predicting 1 when the truth is zero and vice versa)
dat$RMSE[dat$RMSE>1] = 1
### extract population names
POP_NAMES = unique(sort(as.character(dat$TRAIN_POP)))
### remove "ACROSS" SAMPLING SCHEME in TRAIN_POP and TEST_POP, becuase we will be merging only the intra-population data
dat = droplevels(dat[(dat$SAMPLING_SCHEME!="ACROSS") & (nchar(as.character(dat$TEST_POP))<10) , ])
### remove auto-cross-validation
dat = droplevels(dat[dat$TRAIN_POP != dat$TEST_POP, ])
### for each POP_NAMES dataset, extract GPAS performances, i.e. mean RMSE across all cross validations and recalculated percent QTL detected and percent false positives
for (i in 1:length(POP_NAMES)){
    # i = 29
    # print(i)
    list_pop = unlist(strsplit(POP_NAMES[i], ":"))
    sub_dat = droplevels(dat[as.character(dat$TRAIN_POP) %in% list_pop, ])
    ######################
    ### GP PERFORMANCE ###
    ######################
    for (func in c(mean, min, max, var)){
        # func = mean
        rmse_agg = aggregate(RMSE ~ GENOTYPING_SCHEME + MODEL + COVARIATE, FUN=func, data=sub_dat)
        if (exists("rmse_merged")==FALSE){
            rmse_merged = rmse_agg
        } else {
            rmse_merged = merge(rmse_merged, rmse_agg, by=c("GENOTYPING_SCHEME", "MODEL", "COVARIATE"), all=TRUE)
        }
    }
    colnames(rmse_merged)[4:7] = c("RMSE_MEAN", "RMSE_MIN", "RMSE_MAX", "RMSE_VAR")

    ########################
    ### GWAS PERFORMANCE ###
    ########################
    ### TPR, FPR, AUC summstats, and P(QTL fixed) summstats
    ### convert NA in the detected QTL as "NADA" to indicate no QTL were detected
    sub_dat$BONFERRONI_5PERCENT_QTL_ID = as.character(sub_dat$BONFERRONI_5PERCENT_QTL_ID)
    sub_dat$BONFERRONI_5PERCENT_QTL_ID[is.na(sub_dat$BONFERRONI_5PERCENT_QTL_ID)] = "NADA"
    ### define the function to calculate true positive rate
    func_tpr = function(x, nloc=1){
        # x = sub_dat$BONFERRONI_5PERCENT_QTL_ID[1:10]
        loc_names = unique(unlist(strsplit(as.character(x), ":")))
        counts = length(loc_names[loc_names != "NADA"])### "NADA" when no QTL were detected
        return(counts/nloc)
    }
    ### define the function to calculate false positive rate (with LD at 1kb)
    func_fpr = function(x, LD_kb=1, nloc=10000){
        # x = sub_dat$BONFERRONI_5PERCENT_FALSE_POSITIVE[1:10]
        loc_names = unique(unlist(strsplit(as.character(x), ":")))
        if(length(loc_names)>0){ ### if no false positives were detected
            # print(length(loc_names))
            ### group into LD blocks
            split_names = unlist(strsplit(loc_names, "_"))
            loc_coor = matrix(as.numeric(split_names[!is.na(split_names)]), ncol=2, byrow=TRUE)
            loc_coor = loc_coor[!is.na(loc_coor[,1]), ]
            if(length(loc_coor)>2){ ### single locus (chrom,pos)
                counts = 0
                for (chrom in unique(loc_coor[,1])){
                    # chrom=1
                    sub_loc = loc_coor[loc_coor[,1]==chrom, ]
                    i = 1
                    idx_out = c(i)
                    # print(sub_loc)
                    if(length(sub_loc)>2){ ### if only one locus was extracted
                        while (i < nrow(sub_loc)){
                            # i = 1
                            x1 = sub_loc[i,2]
                            i = i + 1
                            for (j in i:nrow(sub_loc)){
                                x2 = sub_loc[j,2]
                                if (abs(x1-x2)<=LD_kb){
                                    i = j
                                    break
                                }
                            }
                            idx_out = c(idx_out, i)
                        }
                    }
                    counts = counts + length(unique(idx_out))
                }
            } else {
                counts = 1 ### single locus
            }
        } else {
            counts = 0 ### no false positive
        }
        # out = counts
        out = counts/nloc
        return(out)
    }
    NQTL = unique(sub_dat$NQTL)
    NLOCI = unique(sub_dat$NLOCI)
    tpr_agg = aggregate(BONFERRONI_5PERCENT_QTL_ID         ~ GENOTYPING_SCHEME + MODEL + COVARIATE, FUN=func_tpr, nloc=NQTL,       data=sub_dat)
    fpr_agg = tryCatch(aggregate(BONFERRONI_5PERCENT_FALSE_POSITIVE ~ GENOTYPING_SCHEME + MODEL + COVARIATE, FUN=func_fpr, nloc=NLOCI-NQTL, data=sub_dat),
                       error=function(e){x = aggregate(REPLICATE ~ GENOTYPING_SCHEME + MODEL + COVARIATE, FUN=mean, data=sub_dat); x[,4]=NA; return(x)})
    colnames(tpr_agg)[4] = "TRUE_POSITIVE_RATE"
    colnames(fpr_agg)[4] = "FALSE_POSITIVE_RATE"
    for (func in c(mean, min, max, var)){
        # func = mean
        auc_agg = aggregate(AUC ~ GENOTYPING_SCHEME + MODEL + COVARIATE, FUN=func, data=sub_dat)
        if (exists("auc_merged")==FALSE){
            auc_merged = auc_agg
        } else {
            auc_merged = merge(auc_merged, auc_agg, by=c("GENOTYPING_SCHEME", "MODEL", "COVARIATE"), all=TRUE)
        }
        # fixedqtl_agg = aggregate(FRACTION_QTL_FIXED ~ GENOTYPING_SCHEME + MODEL + COVARIATE, FUN=func, data=sub_dat)
        # if (exists("fixedqtl_merged")==FALSE){
        #     fixedqtl_merged = fixedqtl_agg
        # } else {
        #     fixedqtl_merged = merge(fixedqtl_merged, fixedqtl_agg, by=c("GENOTYPING_SCHEME", "MODEL", "COVARIATE"), all=TRUE)
        # }
    }
    colnames(auc_merged)[4:7] = c("AUC_MEAN", "AUC_MIN", "AUC_MAX", "AUC_VAR")
    # colnames(fixedqtl_merged)[4:7] = c("FIXEDQTL_MEAN", "FIXEDQTL_MIN", "FIXEDQTL_MAX", "FIXEDQTL_VAR")


    #############
    ### MERGE ###
    #############
    merged = merge(
                    merge(
                            merge(rmse_merged, auc_merged, by=c("GENOTYPING_SCHEME", "MODEL", "COVARIATE"), all=TRUE), 
                    tpr_agg, by=c("GENOTYPING_SCHEME", "MODEL", "COVARIATE"), all=TRUE),
                fpr_agg, by=c("GENOTYPING_SCHEME", "MODEL", "COVARIATE"), all=TRUE)
    if (exists("out")==FALSE){
        out = data.frame(POP_NAMES = rep(POP_NAMES[i],nrow(merged)), NPOP_MERGED = length(list_pop), merged)
    } else {
        out = rbind(out, data.frame(POP_NAMES = rep(POP_NAMES[i],nrow(merged)), NPOP_MERGED = length(list_pop), merged))
    }
    ### clean-up
    rm("rmse_merged")
    rm("tpr_agg")
    rm("fpr_agg")
    rm("auc_merged")
    rm("fixedqtl_merged")
}
### output
OUT = cbind(dat[1:nrow(out), 1:19], out)
fname_output = paste0(dirname(fname_input), "/GPAS_OUTPUT_NPOP_PERFORMANCE.csv")
write.table(OUT, file=fname_output, sep=",", quote=FALSE, row.names=FALSE)
