#!/bin/bash

#############################################
### FILTER GENOTYPES AND BUILD COVARIATES ###
#############################################

### Input and output specifications:
while [ ! -z "$1" ]; do
  case "$1" in
    --dir|-d) shift
      if [ $(ls $1 | wc -l) == 0 ]; then
        echo "Directory is empty: $1"
      else
        DIR=$1
      fi;;
    --src-dir|-s) shift
      if [ $(ls $1 | wc -l) == 0 ]; then
        echo "Directory is empty: $1"
      else
        SRC_DIR=$1
      fi;;
    --n-threads|-t) shift
      NTHREADS=$1;;
    --help|-h) shift
      echo "#############################################"
      echo "### FILTER GENOTYPES AND BUILD COVARIATES ###"
      echo "#############################################"
      echo ""
      echo ""
      echo '##############'
      echo '### INPUTS ###'
      echo '##############'
      echo "--dir|-d          directory containing all the ouputs of GPASim_3.0_across_population_sampling.sh"
      echo "--src-dir|-s      directory where plink, gemma, gcta and genomic_prediction/ are located"
      echo "--n-threads|-t    number of threads to use in parallel"
      echo "--help|-h         help documentation"
      echo ""
      echo ""
      echo '###############'
      echo '### OUTPUTS ###'
      echo '###############'
      echo 'MAF-filtered bed, bim, & fam files (old files were over-written):'
      echo '(1) ${prefix}.bed'
      echo '(2) ${prefix}.bim'
      echo '(3) ${prefix}.fam'
      echo 'Headerless square symmetric kinship matrices and GCTA GRM files:'
      echo '(4) ${prefix}_STANDARDIZED.kin (GCTA-specific GRM files: ${prefix}_STANDARDIZED.grm.bin,${prefix}_STANDARDIZED.grm.id, & ${prefix}_STANDARDIZED.grm.N.bin)'
      echo '(5) ${prefix}_GRM.kin (GCTA-specific GRM files: ${prefix}_GRM.grm.bin,${prefix}_GRM.grm.id, & ${prefix}_GRM.grm.N.bin)'
      echo 'MAF-filtered sync files (old files were over-written):'
      echo '(6) ${prefix}.sync'
      echo 'Headerless square symmetric Fst matrices:'
      echo '(7) ${prefix}_HIVERT.fst'
      echo '(8) ${prefix}_WEIRCOCK.fst'
      echo 'These are headerless square symmetric Fst matrices,'
      echo 'where: prefix=$(echo $f | sed "s/.sync//g"); f={$(ls POP_*.sync)}.'
      echo ""
      echo ""
      echo '###############'
      echo '### EXAMPLE ###'
      echo '###############'
      echo 'rep=1'
      echo 'nQTL=10'
      echo 'migration=0.00'
      echo 'selection=0.25'
      echo 'bg_selection=0.00'
      echo 'GRADIENT=1'
      echo 'DIR=/data/Lolium/Quantitative_Genetics/LOLSIM2020/LOLSIM_${rep}rep_${nQTL}qtl_${migration}mr_${selection}fgs_${bg_selection}bgs_${GRADIENT}grad'
      echo 'SRC_DIR=/data/Lolium/Softwares'
      echo 'NTHREADS=12'
      echo 'time \'
      echo '${SRC_DIR}/genomic_prediction/src/GPASim_4.0_filter_build_covariates.sh \'
      echo '      -d $DIR \'
      echo '      -s $SRC_DIR \'
      echo '      -t $NTHREADS'
      echo ""
      echo ""
      exit 0
      ;;
    *)
      echo "What is this?! $1"
      exit 1
      ;;
  esac
shift
done
### Input parameter check
echo -e "--dir\n--src-dir\n--n-threads" > input_parameter_names.temp
echo -e "$DIR\n$SRC_DIR\n$NTHREADS" > input_parameter_values.temp
for i in $(seq 1 $(cat input_parameter_names.temp | wc -l))
do
  input_name=$(head -n${i} input_parameter_names.temp | tail -n1)
  input_value=$(head -n${i} input_parameter_values.temp | tail -n1)
  if [ $(echo $input_value | wc -w) == 0 ]
  then
    echo "Please provide an input for: $input_name"
    exit 1
  fi
done
rm input_parameter_names.temp input_parameter_values.temp
### Set working directories
cd $DIR
#####################
### Indi-seq data ###
#####################
# txt STANDARDIZED to bin STANDARDIZED (subfunction within "GPASim_4.1_build_kinship_Indiseq.sh")
echo '
  args = commandArgs(trailingOnly=TRUE)
  prefix = args[1]
  # prefix = "POP_01"
  id = read.table(paste0(prefix, ".fam"), colClasse=c("character", "character"))[,1:2]
  n = nrow(id)
  N = as.vector(rep(nrow(read.table(paste0(prefix, ".bim"))), time=(((n^2)-n)/2)+n))
  K = as.matrix(read.table(paste0(prefix, "_STANDARDIZED.kin")))
  linearK = as.vector(rep(0, times=(((n^2)-n)/2)+n))
  i = cumsum(1:n)
  linearK[i] = diag(K)
  linearK[-i] = K[upper.tri(K)]
  BinIO = file(paste0(prefix,"_STANDARDIZED.grm.bin"), "wb")
  NIO = file(paste0(prefix,"_STANDARDIZED.grm.N.bin"), "wb")
  IDFileName=paste(prefix,"_STANDARDIZED.grm.id",sep="")
  size=4
  writeBin(linearK, BinIO, size=size)
  writeBin(N, NIO, size=size)
  write.table(id, file=IDFileName, quote=FALSE, row.names=FALSE, col.names=FALSE)
  close(BinIO)
  close(NIO)
' > txt2bin.r
# bin GRM to txt GRM (subfunction within "GPASim_4.1_build_kinship_Indiseq.sh")
echo '
  args = commandArgs(trailingOnly=TRUE)
  prefix = paste0(args[1], "_GRM")
  # prefix = "POP_01_GRM"
  # prefix = "POP_01_STANDARDIZED"
  ReadGRMBin=function(prefix, AllN=F, size=4){
    BinFileName=paste(prefix,".grm.bin",sep="")
    NFileName=paste(prefix,".grm.N.bin",sep="")
    IDFileName=paste(prefix,".grm.id",sep="")
    id = read.table(IDFileName)
    n=dim(id)[1]
    BinFile=file(BinFileName, "rb");
    grm=readBin(BinFile, n=n*(n+1)/2, what=numeric(0), size=size)
    NFile=file(NFileName, "rb");
    if(AllN==T){
      N=readBin(NFile, n=n*(n+1)/2, what=numeric(0), size=size)
    }
    else N=readBin(NFile, n=1, what=numeric(0), size=size)
    i=cumsum(1:n)
    return(list(diag=grm[i], off=grm[-i], id=id, N=N))
  }
  k = ReadGRMBin(prefix)
  n = length(k$diag)
  K = matrix(0, n, n)
  diag(K) = k$diag
  counter = 1
  for (i in 2:n){
    for (j in 1:(i-1)){
      K[i,j] = k$off[counter]
      K[j,i] = k$off[counter]
      counter = counter + 1
    }
  }
  write.table(K, paste0(prefix, ".kin"), row.names=FALSE, col.names=FALSE)
  # gemmaK = read.table(paste0(sub("_GRM", "", prefix), "_STANDARDIZED.kin"), F)
  # cor(K[,1], gemmaK[,1])
' > bin2txt.r
### parallel execution
ls POP_*.bed MULTIPOP_*_*.bed  | sort | uniq | sed 's/.bed//g' > INDISEQ_PREFIXES.txt
### Set the amount of memory to allocate for each instance of plink execution
PLINK_MEM=$(echo $(head /proc/meminfo | grep "MemAvailable" | cut -d: -f2 | sed "s/[ ]//g" | sed "s/kB//g") / ${NTHREADS}000 | bc)
time \
parallel -j ${NTHREADS} ${SRC_DIR}/genomic_prediction/src/GPASim_4.1_build_kinship_Indiseq.sh \
                                    {} \
                                    ${SRC_DIR} \
                                    ${PLINK_MEM} \
                                    ::: $(cat INDISEQ_PREFIXES.txt)
rm txt2bin.r bin2txt.r
### OUTPUTS:
### MAF-filtered bed, bim, & fam files (old files were over-written):
### (1) ${prefix}.bed
### (2) ${prefix}.bim
### (3) ${prefix}.fam
### Headerless square symmetric kinship matrices and GCTA GRM files:
### (4) ${prefix}_STANDARDIZED.kin
### (5) ${prefix}_GRM.kin (GCTA-specific GRM files: ${prefix}_GRM.grm.bin,${prefix}_GRM.grm.id, & ${prefix}_GRM.grm.N.bin)

#####################
### Pool-seq data ###
#####################
ls POP_*.sync MULTIPOP_*_*.sync  | sort | uniq | sed 's/.sync//g' > POOLSEQ_PREFIXES.txt
### parallel computation for within population Pool-seq
time \
parallel -j ${NTHREADS} ${SRC_DIR}/genomic_prediction/src/GPASim_4.2_build_Fst_Poolseq.sh \
                                    {} \
                                    ${SRC_DIR} \
                                    FALSE \
                                    ::: $(grep "^POP_" POOLSEQ_PREFIXES.txt)
### iterative computation for across population Pool-seq to enable parallel Fst computation
if [ $(ls POP_*.sync | wc -l) -lt 20 ]
then
  time \
  parallel -j ${NTHREADS} ${SRC_DIR}/genomic_prediction/src/GPASim_4.2_build_Fst_Poolseq.sh \
                                      {} \
                                      ${SRC_DIR} \
                                      FALSE \
                                      ::: $(grep "^MULTIPOP_" POOLSEQ_PREFIXES.txt)
else
  time \
  for i in $(grep "^MULTIPOP_" POOLSEQ_PREFIXES.txt)
  do
    ${SRC_DIR}/genomic_prediction/src/GPASim_4.2_build_Fst_Poolseq.sh \
                                      ${i} \
                                      ${SRC_DIR} \
                                      TRUE
  done
fi
### OUTPUTS:
### MAF-filtered sync files (old files were over-written):
### (1) ${prefix}.sync
### Headerless square symmetric Fst matrices:
### (2) ${prefix}_HIVERT.fst
### (3) ${prefix}_WEIRCOCK.fst
echo "#########################################################################################"
