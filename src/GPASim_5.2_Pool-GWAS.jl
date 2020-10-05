##############################################
### GENOME-WIDE ASSOCIATION STUDIES (GWAS) ###
### using Pool sequencing (Pool-seq) data  ###
##############################################
### INPUT:
### (1) prefix of the sync file of genotype and phenotype data (POP_* or MULTIPOP_*_*)
### OUTPUTS:
### (1) string(prefix, "_GWAlpha_SNPwise.gwas") (HEADER: |col1:CHROM|col2:POS|col3:ALLELES|col4:FREQ|col5:ALPHA|col6:PVALUES|col7:LOD|)
### (2) string(prefix, "_GWAlpha_REML_LS_", fst_id, ".gwas")
### (3) string(prefix, "_GWAlpha_REML_LS_", fst_id, ".ranef")
### (2) string(prefix, "_GWAlpha_REML_RR_", fst_id, ".gwas")
### (3) string(prefix, "_GWAlpha_REML_RR_", fst_id, ".ranef")
### (2) string(prefix, "_GWAlpha_REML_GELMNET_", fst_id, ".gwas")
### (3) string(prefix, "_GWAlpha_REML_GELMNET_", fst_id, ".ranef")
### (2) string(prefix, "_GWAlpha_REML_LASSO_", fst_id, ".gwas")
### (3) string(prefix, "_GWAlpha_REML_LASSO_", fst_id, ".ranef")
### where:
# xiferp=$(echo $prefix | rev); tsf=$(echo ${fst%.*} | rev)
# fst_id=$(echo ${tsf%_$(echo $xiferp)*} | rev) ### (HIVERT or WEIRCOCK)
# *.gwas: HEADER: |col1:CHROM|col2:POS|col3:ALLELES|col4:FREQ|col5:ALPHA|col6:PVALUES|col7:LOD|
# *.ranef: HEADERLESS: |col1:random effects|

# define the GWAlpha script to run GWAlpha-SNPwise, GWAlpha mixed model ML_LS, REML_RR, REML_GLMNET, and REML_LASSO
using GWAlpha
using DelimitedFiles
# ARGS = ["POP_01", "HIVERT"]
prefix = ARGS[1]
filename_sync = string(prefix, ".sync")
filename_phen_py = string(prefix, ".py")
filename_phen_csv = string(prefix, ".csv")
### Calculate MAF and set DEPTH=1 for simplicity
MAF = 1.00/(2*sum(DelimitedFiles.readdlm(filename_phen_csv, ',')[:,1]))
DEPTH = 1
println("############################")
println("### GWAlpha: SNP-wise ML ###")
println("############################")
if isfile(string(prefix, "_GWAlpha_SNPwise.gwas"))==false
  @time GWAlpha.PoolGPAS(filename_sync=filename_sync, filename_phen=filename_phen_py, model="GWAlpha", maf=MAF, depth=DEPTH, fpr=0.01)
  # output: comma-delimited; HEADER: |col1:CHROM|col2:POS|col3:ALLELES|col4:FREQ|col5:ALPHA|col6:PVALUES|col7:LOD|
  mv(string(prefix, "-GWAlpha-OUTPUT.csv"), string(prefix, "_GWAlpha_SNPwise.gwas"))
end
println("###########################")
println("### GWAlpha: MIXED REML ###")
println("###########################")
### Fst file to use (headerless, comma-delimited, square-symmetric matrix) (HIVERT or WEIRCOCK)
for fst_id in ["WEIRCOCK", "HIVERT"]
  filename_random_covariate = string(prefix, "_", fst_id, ".fst")
  @time GWAlpha.PoolGPAS(filename_sync=filename_sync, filename_phen=filename_phen_csv, model="MIXED", filename_random_covariate=filename_random_covariate, varcomp_est="REML", maf=MAF, depth=DEPTH, fpr=0.01)
  # fixed effects output: comma-delimited; HEADER: |col1:CHROM|col2:POS|col3:ALLELES|col4:FREQ|col5:ALPHA|col6:PVALUES|col7:LOD|
  mv(string(prefix, "_MAF", MAF,"_DEPTH", DEPTH, "-MIXEDREML_USER_DEFINED_COVARIATE-OUTPUT.csv"), string(prefix, "_GWAlpha_MIXED_REML_", fst_id, ".gwas"), force=true)
  # random effects output: comma-delimited; HEADERLESS: |col1:random effects|
  mv(string(prefix, "_MAF", MAF,"_DEPTH", DEPTH, "-MIXEDREML_USER_DEFINED_COVARIATE-RANEF-OUTPUT.csv"), string(prefix, "_GWAlpha_MIXED_REML_", fst_id, ".ranef"), force=true)
end
println("#######################")
println("### GWAlpha: RIDGE ###")
println("######################")
ALPHA = 0.00
@time GWAlpha.PoolGPAS(filename_sync=filename_sync, filename_phen=filename_phen_csv, model="GLMNET", glmnet_alpha=ALPHA, maf=MAF, depth=DEPTH, fpr=0.01)
# fixed effects output: comma-delimited; HEADER: |col1:CHROM|col2:POS|col3:ALLELES|col4:FREQ|col5:ALPHA|col6:PVALUES|col7:LOD|
mv(string(prefix, "_MAF", MAF,"_DEPTH", DEPTH, "-GLMNET_ALPHA", ALPHA, "-OUTPUT.csv"), string(prefix, "_GWAlpha_RIDGE.gwas"), force=true)
println("#######################")
println("### GWAlpha: LASSO ###")
println("######################")
ALPHA = 1.00
@time GWAlpha.PoolGPAS(filename_sync=filename_sync, filename_phen=filename_phen_csv, model="GLMNET", glmnet_alpha=ALPHA, maf=MAF, depth=DEPTH, fpr=0.01)
# fixed effects output: comma-delimited; HEADER: |col1:CHROM|col2:POS|col3:ALLELES|col4:FREQ|col5:ALPHA|col6:PVALUES|col7:LOD|
mv(string(prefix, "_MAF", MAF,"_DEPTH", DEPTH, "-GLMNET_ALPHA", ALPHA, "-OUTPUT.csv"), string(prefix, "_GWAlpha_LASSO.gwas"), force=true)
println("#######################")
println("### GWAlpha: GLMNET ###")
println("#######################")
ALPHA = 0.50
@time GWAlpha.PoolGPAS(filename_sync=filename_sync, filename_phen=filename_phen_csv, model="GLMNET", glmnet_alpha=ALPHA, maf=MAF, depth=DEPTH, fpr=0.01)
# fixed effects output: comma-delimited; HEADER: |col1:CHROM|col2:POS|col3:ALLELES|col4:FREQ|col5:ALPHA|col6:PVALUES|col7:LOD|
mv(string(prefix, "_MAF", MAF,"_DEPTH", DEPTH, "-GLMNET_ALPHA", ALPHA, "-OUTPUT.csv"), string(prefix, "_GWAlpha_GLMNET.gwas"), force=true)
