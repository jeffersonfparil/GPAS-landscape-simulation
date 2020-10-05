#!/bin/bash

##################################
### ACROSS POPULATION SAMPLING ###
##################################

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
    --n-libraries|-l) shift
      NLIB=$1;;
    --help|-h) shift
      echo "##################################"
      echo "### ACROSS POPULATION SAMPLING ###"
      echo "##################################"
      echo ""
      echo ""
      echo '##############'
      echo '### INPUTS ###'
      echo '##############'
      echo "--dir|-d          directory containing all the ouputs of GPASim_2.0_parse.sh"
      echo "--src-dir|-s      directory where plink and genomic_prediction/ are located"
      echo "--n-threads|-t    number of threads to use in parallel"
      echo "--n-libraries|-l  maximum number of individually barcoded libraries for Indi-seq or Pool-seq"
      echo "--help|-h         help documentation"
      echo ""
      echo ""
      echo '###############'
      echo '### OUTPUTS ###'
      echo '###############'
      echo '(1) string("MULTIPOP_", population_count, "_", hash(sample), ".idx") (population ID and sample size per population: headerless, tab-delimited, |col1:pop_ID|col2:sample_sizes|)'
      echo '(2) string("MULTIPOP_", population_count, "_", hash(sample), ".bed")'
      echo '(3) string("MULTIPOP_", population_count, "_", hash(sample), ".bim")'
      echo '(4) string("MULTIPOP_", population_count, "_", hash(sample), ".fam")'
      echo '(5) string("MULTIPOP_", population_count, "_", hash(sample), ".syn")'
      echo '(6) string("MULTIPOP_", population_count, "_", hash(sample), ".py")'
      echo '(7) string("MULTIPOP_", population_count, "_", hash(sample), ".csv")'
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
      echo 'NLIB=1000'
      echo 'time \'
      echo '${SRC_DIR}/genomic_prediction/src/GPASim_3.0_across_population_sampling.sh \'
      echo '      -d $DIR \'
      echo '      -s $SRC_DIR \'
      echo '      -t $NTHREADS \'
      echo '      -l $NLIB'
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
echo -e "--dir\n--src-dir\n--n-threads\n--n-libraries" > input_parameter_names.temp
echo -e "$DIR\n$SRC_DIR\n$NTHREADS\n$NLIB" > input_parameter_values.temp
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
### Set the amount of memory to allocate for each instance of plink execution
PLINK_MEM=$(echo $(head /proc/meminfo | grep "MemAvailable" | cut -d: -f2 | sed "s/[ ]//g" | sed "s/kB//g") / ${NTHREADS}000 | bc)
### Set working directories
cd $DIR
### List population names
ls POP_*.bed | sort -g | uniq | sed 's/.bed//g' > POP_PREFIXES.temp
#########################################
### Layout in a 2-D array given that: ###
#########################################
###   - quantiNemo2 simulation specifications defined in GPASim_1.0_simulate.sh generates a square landscape, and
###   - assuming the populations are arranged by row across this square landscape
### NOTE: divide the landscape into equally-sized rectangles and sample the central population in a jittered way
###       then list the sample sizes for each population constrained by the number of libraries (nLib)
time \
julia ${SRC_DIR}/genomic_prediction/src/GPASim_3.1_sample_populations.jl \
    POP_PREFIXES.temp \
    ${NLIB}
### OUTPUTS:
### (1n) string("MULTIPOP_", population_count, "_", hash(sample), ".idx") (population ID and sample size per population: headerless, tab-delimited, |col1:pop_ID|col2:sample_sizes|)
#######################################
### Across population Indi-seq data ###
#######################################
### List across population sampling names
ls MULTIPOP_*_*.idx | sort -g | uniq | sed 's/.idx//g' > MULTIPOP_PREFIXES.temp
### Iterate per sample
rm filter.temp
for prefix in $(cat MULTIPOP_PREFIXES.temp)
do
  # prefix=POP_01
  touch filter.temp
  # random sampling of individuals per population
  while read l
  do
    # l=$(head -n1 $f)
    pop=$(echo $l | cut -d' ' -f1)
    size=$(echo $l | cut -d' ' -f2)
    shuf -n${size} ${pop}.fam | cut -d' ' -f1-2 | sort >> filter.temp
  done < ${prefix}.idx
  # generate bed, bim, fam files
  ${SRC_DIR}/plink --bed OUT.bed \
                    --bim OUT.bim \
                    --fam OUT.fam \
                    --keep filter.temp \
                    --make-bed \
                    --threads ${NTHREADS} \
                    --memory ${PLINK_MEM} \
                    --out ${prefix}
  rm filter.temp *.log
done
### OUTPUTS:
### (1n) string("MULTIPOP_", population_count, "_", hash(sample), ".bed")
### (2n) string("MULTIPOP_", population_count, "_", hash(sample), ".bim")
### (3n) string("MULTIPOP_", population_count, "_", hash(sample), ".fam")
#######################################
### Across population Pool-seq data ###
#######################################
### build phenotype files and sort pop_pooling_info.temp.${prefix} (subfunction within "/GPASim_3.2_merge_onepool_sync_py_csv.sh")
echo '
  args = commandArgs(trailingOnly=TRUE)
  fname = args[1]
  out_prefix = args[2]
  dat = read.delim(fname, header=FALSE)
  dat = dat[order(dat$V3, decreasing=FALSE), ]
  total_variance = 0
  for (i in dat$V1){
    total_variance = total_variance + var(read.delim(paste0(sub("ONEPOOL_", "", i), ".fam"), sep=" ", header=FALSE)$V6)
  }
  perc = cumsum(dat$V2) / sum(dat$V2)
  line1 = "Pheno_name=\"pheno\";"
  line2 = paste0("sig=", sqrt(total_variance), ";")
  line3 = paste0("MIN=", min(dat$V3), ";")
  line4 = paste0("MAX=", max(dat$V3), ";")
  line5 = paste0("perc=[", paste(perc[1:length(perc)-1], collapse=","), "];")
  line6 = paste0("q=[", paste(dat$V3[1:nrow(dat)-1], collapse=","), "];")
  PY_PHEN = rbind(line1, line2, line3, line4, line5, line6)
  CSV_PHEN = dat[, 2:3]
  write.table(PY_PHEN, paste0(out_prefix, ".py"), quote=FALSE, row.names=FALSE, col.names=FALSE)
  write.table(CSV_PHEN, paste0(out_prefix, ".csv"), sep=",", quote=FALSE, row.names=FALSE, col.names=FALSE)
  write.table(dat, fname, quote=FALSE, row.names=FALSE, col.names=FALSE)
' > build_pheno_and_sort.r
### parallel execution
time \
parallel ${SRC_DIR}/genomic_prediction/src/GPASim_3.2_merge_onepool_sync_py_csv.sh {} ::: $(cat MULTIPOP_PREFIXES.temp)
rm build_pheno_and_sort.r
### OUTPUTS:
### (1n) string("MULTIPOP_", population_count, "_", hash(sample), ".syn")
### (2n) string("MULTIPOP_", population_count, "_", hash(sample), ".py")
### (3n) string("MULTIPOP_", population_count, "_", hash(sample), ".csv")
