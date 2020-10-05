#!/bin/bash
###############################################
### ESTIMATE FST MATRICES FOR POOL-SEQ DATA ###
###############################################

### INPUTS:
prefix=$1
# prefix=POP_01
SRC_DIR=$2 ### directory where plink and genomic_prediction/ are located
PARALLEL_FST=$3 ### perform parallel Fst computation?

### OUTPUTS:
### MAF-filtered sync files (old files were renamed with the prefix: FULLSNPSET_*):
### (1) ${prefix}.sync
### Headerless square symmetric Fst (relatedness) matrices:
### (2) ${prefix}_HIVERT.fst
### (3) ${prefix}_WEIRCOCK.fst

### Define julia script that will call GWAlha.jl filter_sync() function
echo -e "
  ### INPUTS:
  fname_sync = ARGS[1] ### genotype file Pool-seq
  fname_phen_csv = ARGS[2] ### phenotype file (headerless: |col1:pool_sizes|col2:mean_phenotype|)
  ### OUTPUT:
  ### (1) string(join(split(filename_sync, '.')[1:(end-1)], '.'), '_MAF', MAF, '_DEPTH', DEPTH, '.sync')
  ### Load libraries
  using DelimitedFiles
  using GWAlpha
  ### Calculate MAF and set DEPTH=1 for simplicity
  MAF = 1.00/(2*sum(DelimitedFiles.readdlm(fname_phen_csv, ',')[:,1]))
  DEPTH = 1
  ### Filter sync
  GWAlpha.sync_processing_module.sync_filter(filename_sync=fname_sync, MAF=MAF, DEPTH=DEPTH)
" > ${prefix}_filter_sync.jl

### Define julia script that will call GWAlha.jl Fst_pairwise() function
if [ $PARALLEL_FST == TRUE ]
then
echo -e "
  ### INPUTS:
  fname_sync = ARGS[1] ### genotype file Pool-seq
  fname_csv = ARGS[2] ### phenotype file Pool-seq (csv format)
  METHOD = ARGS[3] ### Fst estimation method: Hivert or WeirCock or Relatedness
  ### OUTPUT:
  ### (1) string(join(split(fname_sync, '.')[1:(end-1)], '.'), '_COVARIATE_FST.csv')
  ### Load libraries
  using DelimitedFiles
  using Distributed
  Distributed.addprocs(length(Sys.cpu_info()))
  @everywhere using GWAlpha
  ### Load pool sizes
  pool_sizes = convert(Array{Int64,1}, DelimitedFiles.readdlm(fname_csv, ',')[:,1])
  ### Estimate Fst matrix
  if METHOD != \"Relatedness\"
    GWAlpha.relatedness_module.Fst_pairwise(filename_sync=fname_sync, window_size=100000, pool_sizes=pool_sizes, METHOD=METHOD)
  else
    GWAlpha.relatedness_module.standardized_relatedness(fname_sync)
  end
" > ${prefix}_compute_poolFST.jl
else
  echo -e "
    ### INPUTS:
    fname_sync = ARGS[1] ### genotype file Pool-seq
    fname_csv = ARGS[2] ### phenotype file Pool-seq (csv format)
    METHOD = ARGS[3] ### Fst estimation method: Hivert or WeirCock or Relatedness
    ### OUTPUT:
    ### (1) string(join(split(fname_sync, '.')[1:(end-1)], '.'), '_COVARIATE_FST.csv')
    ### Load libraries
    using DelimitedFiles
    using GWAlpha
    ### Load pool sizes
    pool_sizes = convert(Array{Int64,1}, DelimitedFiles.readdlm(fname_csv, ',')[:,1])
    ### Estimate Fst matrix
    if METHOD != \"Relatedness\"
      GWAlpha.relatedness_module.Fst_pairwise(filename_sync=fname_sync, window_size=100000, pool_sizes=pool_sizes, METHOD=METHOD)
    else
      GWAlpha.relatedness_module.standardized_relatedness(fname_sync)
    end
  " > ${prefix}_compute_poolFST.jl
fi

### Filter sync files by MAF=1/pop_size and DEPTH=1 (for simplicity)
julia ${prefix}_filter_sync.jl ${prefix}.sync ${prefix}.csv
### rename original dataset (to be used in genomic prediction cross-validation)
mv ${prefix}.sync FULLSNPSET_${prefix}.sync
### rename filtered dataset
mv ${prefix}_MAF*_DEPTH*.sync ${prefix}.sync

### Compute Fst using Hivert's method
julia ${prefix}_compute_poolFST.jl ${prefix}.sync ${prefix}.csv Hivert
mv ${prefix}_COVARIATE_FST.csv ${prefix}_HIVERT.fst

### Compute Fst using Weir & Cockerham's method
julia ${prefix}_compute_poolFST.jl ${prefix}.sync ${prefix}.csv WeirCock
mv ${prefix}_COVARIATE_FST.csv ${prefix}_WEIRCOCK.fst

# ### Compute simple standardized relatedness matrix
# julia ${prefix}_compute_poolFST.jl ${prefix}.sync ${prefix}.csv Relatedness
# mv ${prefix}_COVARIATE_RELATEDNESS.csv ${prefix}_RELATEDNESS.fst

### Clean-up
rm ${prefix}_filter_sync.jl ${prefix}_compute_poolFST.jl