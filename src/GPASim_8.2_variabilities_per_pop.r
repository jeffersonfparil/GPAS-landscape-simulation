### Extract genetic variability per population per landscape
args = commandArgs(trailingOnly=TRUE)
# args = c("/data/Lolium/Quantitative_Genetics/LOLSIM_2020/GPASim_5REP_50NQTL_100NBGS_0.01MIG_0.50SEL_0.10BGSEL_2GRAD")
# args = c("/data/Lolium/Quantitative_Genetics/LOLSIM_2020/GPASim_1REP_100NQTL_100NBGS_0.01MIG_0.50SEL_0.10BGSEL_1GRAD")
DIR = args[1]
print(DIR)
QTL_spec = read.csv(paste0(DIR, "/OUTPUT/QTL_SPEC.csv"))
QTL_freq = read.delim(paste0(DIR, "/RESOURCES/OUT.frq.csv"), header=TRUE, sep=",")

list_fnames_fam = system(paste0("ls ", DIR, "/RESOURCES/POP_*.fam"), inter=TRUE)
list_fnames_hwe = system(paste0("ls ", DIR, "/RESOURCES/POP_*.hwe.csv"), inter=TRUE)
### initialise output vectors
fbasenames = sub(".hwe.csv", "", basename(list_fnames_hwe))
pop_names = c()
pheno_mean = c()
pheno_var = c()
n_qtl_captured = c()
expected_hetero = c()
observed_hetero = c()
qtl_2pq = c()
qtl_maf = c()
for (i in 1:length(fbasenames)){
    # i = 1
    # print(i)
    fam = read.delim(paste0(DIR, "/RESOURCES/", fbasenames[i], ".fam"), header=FALSE, sep=" ")
    if (nrow(fam)>1){
        hwe = read.delim(paste0(DIR, "/RESOURCES/", fbasenames[i], ".hwe.csv"), header=TRUE, sep=",")
        frq = droplevels(QTL_freq[QTL_freq$CLST==as.numeric(sub("POP_", "", fbasenames[i])), ])
        ### use only the QTL
        hwe = hwe[as.character(hwe$SNP) %in% unique(paste0(as.character(QTL_spec$CHROM), "_", as.character(QTL_spec$POS))), ]
        frq = frq[as.character(frq$SNP) %in% paste0(as.character(QTL_spec$CHROM), "_", as.character(QTL_spec$POS)), ]
        ### output
        pop_names = c(pop_names, fbasenames[i])
        pheno_mean = c(pheno_mean, mean(fam$V6, na.rm=TRUE))
        pheno_var = c(pheno_var, var(fam$V6, na.rm=TRUE))
        n_qtl_captured = c(n_qtl_captured, nrow(hwe))
        expected_hetero = c(expected_hetero, mean(hwe$E.HET., na.rm=TRUE))
        observed_hetero = c(observed_hetero, mean(hwe$O.HET., na.rm=TRUE))
        ### expected heterozygosity of all the QTL not just the ones capture by genotyping after filtering
        qtl_2pq = c(qtl_2pq, mean(2*frq$MAF*(1-frq$MAF), na.rm=TRUE)) ### use unfiltered minor allele frequencies from frq instead of hwe
        ### mean minor allele freq at QTL
        causal_qtl_spec = QTL_spec[QTL_spec$EFFEC!=0, ]
        causal_qtl_spec$SNP = paste0(causal_qtl_spec$CHROM, "_", causal_qtl_spec$POS)
        causal_qtl_maf = merge(frq, causal_qtl_spec, by="SNP")
        idx_rev = (causal_qtl_maf$ALLELE=="A" & causal_qtl_maf$A1==2) | 
                  (causal_qtl_maf$ALLELE=="T" & causal_qtl_maf$A1==1) |
                  (causal_qtl_maf$ALLELE=="A" & causal_qtl_maf$A1==0 & causal_qtl_maf$A2==1) |
                  (causal_qtl_maf$ALLELE=="T" & causal_qtl_maf$A1==0 & causal_qtl_maf$A2==2)
        causal_qtl_maf$MAF[idx_rev] = 1 - causal_qtl_maf$MAF[idx_rev]
        qtl_maf = c(qtl_maf, mean(causal_qtl_maf$MAF, na.rm=TRUE))
    }
}
REP = rep(sub("REP", "", unlist(strsplit(basename(DIR), "_"))[2]), length(pop_names))
NQTL = rep(sub("NQTL", "", unlist(strsplit(basename(DIR), "_"))[3]), length(pop_names))
MIG = rep(sub("MIG", "", unlist(strsplit(basename(DIR), "_"))[5]), length(pop_names))
SEL = rep(sub("SEL", "", unlist(strsplit(basename(DIR), "_"))[6]), length(pop_names))
GRAD = rep(sub("GRAD", "", unlist(strsplit(basename(DIR), "_"))[8]), length(pop_names))
OUT = data.frame(REP, NQTL, MIG, SEL, GRAD, 
                 POP=pop_names,
                 PHENO_MEAN=pheno_mean,
                 PHENO_VAR=pheno_var,
                 NQTL_CAPTURED=n_qtl_captured,
                 EXPECTED_HET=expected_hetero,
                 OBSERVED_HET=observed_hetero,
                 QTL_2pq=qtl_2pq,
                 CAUSAL_ALLELE_FREQ=qtl_maf
                )
# par(mfrow=c(2,4))
# plot(OUT$PHENO_MEAN)
# plot(OUT$PHENO_VAR)
# plot(OUT$NQTL_CAPTURED)
# plot(OUT$OBSERVED_HET)
# plot(OUT$EXPECTED_HET)
# plot(OUT$QTL_2pq)
# plot(OUT$CAUSAL_ALLELE_FREQ)
write.table(OUT, file=paste0(DIR, "/OUTPUT/VARIATION_PER_POP.csv"), sep=",", quote=FALSE, row.names=FALSE, col.names=TRUE)