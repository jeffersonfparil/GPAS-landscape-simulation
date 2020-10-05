#!/bin/bash

##################################
### PARSING QUANTINEMO2 OUTPUT ###
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
    --n-pools|-p) shift
      NPOOLS=$1;;
    --help|-h) shift
      echo "##################################"
      echo "### PARSING QUANTINEMO2 OUTPUT ###"
      echo "##################################"
      echo ""
      echo ""
      echo '##############'
      echo '### INPUTS ###'
      echo '##############'
      echo "--dir|-d        directory containing all the ouputs of GPASim_1.0_simulate.sh"
      echo "--src-dir|-s    directory where plink and genomic_prediction/ are located"
      echo "--n-threads|-t  number of threads to use in parallel"
      echo "--n-pools|-p    number of pools per population"
      echo "--help|-h       help documentation"
      echo ""
      echo ""
      echo '###############'
      echo '### OUTPUTS ###'
      echo '###############'
      echo '(01) POP_*.bed (genotype data in binary file format per population; NOTE: biallelic only!)'
      echo '(02) POP_*.bim (loci identities per population; headerless; tab-delimited; |col1:chrom|col2:varID|col3:pos_cM|col4:pos_bp|col5:ref|col6:alt|)'
      echo '(03) POP_*.fam (individual identities per population; headerless; tab-delimited; |col1:family_ID|col2:within_family_ID|col3:within_family_ID_male|col4:within_family_ID_female|col5:sex|col6:phenotype_0to1_range|)'
      echo '(04) POP_*.sync (Pool-seq genotype data with NPOOLS pools per population; headerless; tab-delimited; |col1:chrom|col2:pos_bp|col3:ref|col4-colN:allele_counts_delimited_by_colon_ATCGNDel)'
      echo '(05) POP_*.py (Pool-seq phenotype data with NPOOLS pools per population; headerless; 1 column; |line1:phenotype_name|line2:min_pheno|line3:max_pheno|line4:nth_percentiles|line5:percentile_values|)'
      echo '(06) POP_*.csv (Pool-seq phenotype data with NPOOLS pools per population; headerless; comma-separated; |col1:pool_sizes|col2:phenotype_means|)'
      echo '(07) ONEPOOL_POP_*.sync (Pool-seq genotype data with 1 pool per population; headerless; tab-delimited; |col1:chrom|col2:pos_bp|col3:ref|col4:allele_counts_delimited_by_colon_ATCGNDel)'
      echo '(08) ONEPOOL_POP_*.py (Pool-seq phenotype data with 1 pool per population; headerless; 1 column; |line1:phenotype_name|line2:min_pheno|line3:max_pheno|line4:nth_percentiles|line5:percentile_values|)'
      echo '(09) ONEPOOL_POP_*.csv (Pool-seq phenotype data with 1 pool per population; headerless; comma-separated; |col1:pool_sizes|col2:phenotype_means|)'
      echo '(10) OUT.bed (genotype data in binary file format for the whole landscape; NOTE: biallelic only!)'
      echo '(11) OUT.bim (loci identities for the whole landscape; headerless; tab-delimited; |col1:chrom|col2:varID|col3:pos_cM|col4:pos_bp|col5:ref|col6:alt|)'
      echo '(12) OUT.fam (individual identities for the whole landscape; headerless; tab-delimited; |col1:family_ID|col2:within_family_ID|col3:within_family_ID_male|col4:within_family_ID_female|col5:sex|col6:phenotype_0to1_range|)'
      echo '(13) OUT.pheno (phenotype data for the whole landscape; first 2 rows are comments; tab-delimited; |col1:family_ID|col2:individual_ID|col3:foreground_selection_phenotypes|col4:background_selection_phenotypes|)'
      echo '(14) LANDSCAPE.stat (phenotypic summary statistics per population; comma-separated; |col1:POP|col2:SIZE|col3:MEAN|col4:VAR|)'
      echo '(15) LANDSCAPE.fst (pairwise fixation indices; comma-separated; headerless)'
      echo '(16) LANDSCAPE.png (heatmaps of phenotype means, variances and Fst; portable network graphics (.png) image format)'
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
      echo 'NPOOLS=5'
      echo 'time \'
      echo '${SRC_DIR}/genomic_prediction/src/GPASim_2.0_parse.sh \'
      echo '      -d $DIR \'
      echo '      -s $SRC_DIR \'
      echo '      -t $NTHREADS \'
      echo '      -p $NPOOLS'
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
echo -e "--dir\n--src-dir\n--n-threads\n--n-pools" > input_parameter_names.temp
echo -e "$DIR\n$SRC_DIR\n$NTHREADS\n$NPOOLS" > input_parameter_values.temp
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
SUB_DIR=$(ls -d $(basename $DIR)_*/)
mv ${SUB_DIR}/* .
rm -R ${SUB_DIR}/
### Convert ped, map, & pheno files into bed, bim, & fam files
# renaming for consistency
mv *.ped OUT.ped
mv *.pheno OUT.pheno
# correcting map file bug (BUG: includes only the QTL for first trait if simulating multiple traits)
tail -n+2 GENOME_SPEC.csv | sed 's/,/\t/g' | cut -f1 > chrom.temp
tail -n+2 GENOME_SPEC.csv | sed 's/,/\t/g' | cut -f3 > pos.temp
tail -n+2 GENOME_SPEC.csv | sed 's/,/\t/g' | cut -f2-3 > cMpos.temp
paste -d_ chrom.temp pos.temp > chrompos.temp
paste chrom.temp chrompos.temp cMpos.temp > OUT.map
# convert ped & map files into bed, bim & fam files
${SRC_DIR}/plink --file OUT --make-bed --threads ${NTHREADS} --memory ${PLINK_MEM} --out OUT
# correct phenotype data:
# - replace fitness values from the ped file into the 1st quantitative phenotype data, and
# - restrict between zero and one
cut -d' ' -f1-5 OUT.fam > OUT.fam.col1-5.temp
tail -n+3 OUT.pheno | cut -d' ' -f4 > OUT.fam.col6.temp
echo '
  dat = read.delim("OUT.fam.col6.temp", header=FALSE)
  min_max = read.delim("QTL_min_max_GEBV.spec")
  true_min = min(c(dat$V1, min_max$MIN_GEBV))
  true_max = max(c(dat$V1, min_max$MAX_GEBV))
  min_max$MIN_GEBV = true_min
  min_max$MAX_GEBV = true_max
  dat$V1 = (dat$V1 - true_min) / (true_max - true_min)
  write.table(min_max, "QTL_min_max_GEBV.spec", quote=FALSE, row.names=FALSE)
  write.table(dat, "OUT.fam.col6.temp", quote=FALSE, col.names=FALSE, row.names=FALSE)
' > map_pheno_into_0_to_1.r
Rscript map_pheno_into_0_to_1.r
rm map_pheno_into_0_to_1.r
paste -d' ' OUT.fam.col1-5.temp OUT.fam.col6.temp > OUT.fam
# clean-up
rm *.temp *.log *.ped *.map

### separate data per population (NOTE: An issue here is that bim files are actually duplicated since it's all the same across populations may hog-up space)
echo "#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#"
echo "Parsing Indi-seq data per population"
echo "#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#"
for fam_id in $(cut -d' ' -f1 OUT.fam | sort -g | uniq)
do
  # fam_id=01
  echo $fam_id > filter.temp
  ${SRC_DIR}/plink --bed OUT.bed \
                    --bim OUT.bim \
                    --fam OUT.fam \
                    --keep-fam filter.temp \
                    --make-bed \
                    --threads ${NTHREADS} \
                    --memory ${PLINK_MEM} \
                    --out POP_${fam_id}
done
# clean-up
rm *.log

### convert bed-bim-fam files into vcf
echo "#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#"
echo "Parsing Pool-seq data per population"
echo "#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#"
for pop in $(ls POP_*.bed | sort -g | uniq | sed 's/.bed//g')
do
  # pop=POP_01
  ${SRC_DIR}/plink  --bed ${pop}.bed \
                    --bim ${pop}.bim \
                    --fam ${pop}.fam \
                    --recode vcf \
                    --threads ${NTHREADS} \
                    --memory ${PLINK_MEM} \
                    --out ${pop}
done

### convert vcf into sync
echo '#!/bin/bash' > vcf2sync.sh
echo "
  i=\$1
  NPOOLS=\$2
  if [ \$NPOOLS == 1 ]
  then
    ### add prefix to single pool per population
    cp \${i}.vcf ONEPOOL_\${i}.vcf
    julia ${SRC_DIR}/genomic_prediction/src/GPASim_2.1_vcf2sync.jl ONEPOOL_\${i}.vcf \${i}.fam \${NPOOLS}
    rm ONEPOOL_\${i}.vcf
  else
    julia ${SRC_DIR}/genomic_prediction/src/GPASim_2.1_vcf2sync.jl \${i}.vcf \${i}.fam \${NPOOLS}
  fi
" >> vcf2sync.sh; chmod +x vcf2sync.sh
# 5 pools per population for within population sampling with Pool-seq genotyping
time parallel ./vcf2sync.sh {} $NPOOLS ::: $(ls POP_*.vcf | sort -g | uniq | sed 's/.vcf//g')
# 1 pool per population for across population sampling with Pool-seq genotyping
# NOTE: We are assuming we can accurately measure allele frequencies with whole population Pool-seq!
time parallel ./vcf2sync.sh {} 1 ::: $(ls POP_*.vcf | sort -g | uniq | sed 's/.vcf//g')
# OUTPUTS (all headerless):
# (1) POP_*.sync
# (2) POP_*.py
# (3) POP_*.csv
# (4) ONEPOOL_POP_*.sync
# (5) ONEPOOL_POP_*.py
# (6) ONEPOOL_POP_*.csv
# clean-up
rm *.vcf

echo "##########################"
echo "Landscape characterization"
echo "##########################"
### Inputs:
### filename of the Pool-seq data where each population is a single pool
### total number of populations in the landscape
### Ouputs:
### (1) Phenotypic summary statistics per population ("LANDSCAPE.stat"; header: |col1:POP|col2:SIZE|col3:MEAN|col4:VAR|)
### (2) Pairwise fixation indices ("LANDSCAPE.fst"; headerless)
### (3) Heatmaps of phenotype means, variances and Fst ("LANDSCAPE.png"; portable network graphics format)
# build sync file
cut -f1-3 $(ls POP_*.sync | head -n1) > LANDSCAPE.sync.temp
for f in $(ls POP_*.sync)
do
  # f=$(ls POP_*.sync | head -n1)
  cut -f4 ${f} > colADD.temp
  paste LANDSCAPE.sync.temp colADD.temp > LANDSCAPE.sync
  cp LANDSCAPE.sync LANDSCAPE.sync.temp
done
rm LANDSCAPE.sync.temp
# calculate phenotype statistics, pairwise fixation indices and plot heatmaps
nPop=$(echo $(head -n1 LANDSCAPE.sync | awk '{print NF}') - 3 | bc)
julia ${SRC_DIR}/genomic_prediction/src/GPASim_2.2_landscape_statistics.jl LANDSCAPE.sync ${nPop}
