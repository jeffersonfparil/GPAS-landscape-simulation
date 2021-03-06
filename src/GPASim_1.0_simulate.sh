#!/bin/bash

#####################################################################
### SIMULATE GENOTYPE, QTL, AND PHENOTYPE DATA USIONG QUANTINEMO2 ###
#####################################################################

### Input and output specifications:
while [ ! -z "$1" ]; do
  case "$1" in
    --softwares-dir|-s) shift
      if [ $(ls $1 | wc -l) == 0 ]; then
        echo "Folder containg the genomic_prediction git repository does not exist: $1"
      else
        SRC_DIR=$1
      fi;;
    --output-dir|-d) shift
      if [ $(grep $(echo $(basename $1)) <(echo $(ls $(dirname $1))) | wc -l) == 0 ]; then
        echo "Output directory does not exist: $1"
        exit 1
      else
        OUTDIR=$1
      fi;;
    --output-prefix|-p) shift
      OUTPREFIX=$1;;
    --n-individuals|-n) shift
      nIndividuals=$1;;
    --n-loci|-l) shift
      nLoci=$1;;
    --n-QTL|-q) shift
      nQTL=$1;;
    --n-BGS|-b) shift
      nBGS=$1;;
    --n-alleles|-a) shift
      nAlleles=$1;;
    --alleles-dist|-r) shift
      allele_eff_model=$1;;
    --n-generations|-t) shift
      nGen=$1;;
    --n-populations|-P) shift
      nPop=$1;;
    --migration-rate|-m) shift
      migration=$1;;
    --selection-intensity-QTL|-sq) shift
      selection=$1;;
    --selection-intensity-BGS|-sb) shift
      bg_selection=$1;;
    --diffusion-gradient|-g) shift
      GRADIENT=$1;;
    --help|-h) shift
      echo "#####################################################################"
      echo "### SIMULATE GENOTYPE, QTL, AND PHENOTYPE DATA USIONG QUANTINEMO2 ###"
      echo "#####################################################################"
      echo ""
      echo ""
      echo '##############'
      echo '### INPUTS ###'
      echo '##############'
      echo "--softwares-dir|-s              full path to the softwares directory cointaining the quantinemo2 folder and the genomic_prediction folder (the cloned git repository)"
      echo "--output-dir|-d                 existing output directory to dump all the output files and folder"
      echo "--output-prefix|-p              prefix for the quantiNemo2 initiation (*.ini) file and the output folder"
      echo "--n-individuals|-n              number of individuals to simulate"
      echo "--n-loci|-l                     number of loci to simulate (neutral + QTL)"
      echo "--n-QTL|-q                      number of QTL among the all the loci simulated"
      echo "--n-BGS|-b                      number of background selection loci (should be highly polygenic >= 100)"
      echo "--n-alleles|-a                  number of alleles per loci, e.g. 5 for A,T,C,G, and DEL"
      echo "--alleles-dist|-r               probability distribution to sample the allele effects from (i.e. CHISQ, NORM, & UNIF)"
      echo "--n-generations|-t              number of generations to simulate"
      echo "--n-populations|-P              number of populations/subpopulations to simulate (NOTE: must have a natural number square-root)"
      echo "--migration-rate|-m             migration rate across the populations (surrently using the 2D stepping stone model see line 117) [@Wang2017]"
      echo "--selection-intensity-QTL|-sq   proportion of individual selected based on trait of interest (nQTL) defined as proportion of the most fit individuals selected, where fitness is defined using the generalized logistic curve (Richards, 1959): ranges from 0.00 to 1.00; lower s means more intense selection"
      echo "--selection-intensity-BGS|-sb   proportion of individual selected based on background selection trait (nBGS) defined as proportion of the most fit individuals selected, where fitness is defined using the generalized logistic curve (Richards, 1959): ranges from 0.00 to 1.00; lower s means more intense selection"
      echo "--diffusion-gradient|-g         gradient of QTL (foreground and background selection) across the landscape: 0 for uniform, 1 for first row of populations, and 2 for first and last rows of populations are the sources of the the non-zero (non-wild type) alleles"
      echo "--help|-h                       help documentation"
      echo ""
      echo ""
      echo '###############'
      echo '### OUTPUTS ###'
      echo '###############'
      echo '(01) ${OUTPREFIX}.ini (quantiNemo2 initialization or configuration file see: https://github.com/lucaferretti/npstat/blob/master/NPStat-manual-v1.pdf)'
      echo '(02) Allele_counts_dist.spec (Multi-allelic distributions; hard-coded; headerless; comma-separated; |col1:ALLELE_COUNT_*|col2:frequency|)'
      echo '(03) GENOME_SPEC.csv (simulated genome; HEADER:[CHROM, r, POS]; comma-separated)'
      echo '(04) QTL_SPEC.csv (simulated QTL for foreground selection; HEADER:[CHROM, POS, ALLELE, EFFECT]; comma-separated)'
      echo '(05) BGS.spec (quantiNemo2 allele specificification file for background selection QTL see: https://github.com/lucaferretti/npstat/blob/master/NPStat-manual-v1.pdf)'
      echo '(06) BGS_idx_map_to_genome.spec (indices of each of the locus_allele combinations (i.e. each row) in the BGS.spec file that maps correctly to the simulated genome: GENOME_SPEC.csv)'
      echo '(07) QTL.spec (quantiNemo2 allele specificification file for foreground selection QTL see: https://github.com/lucaferretti/npstat/blob/master/NPStat-manual-v1.pdf)'
      echo '(08) QTL_idx_map_to_genome.spec (indices of each of the locus_allele combinations (i.e. each row) in the QTL.spec file that maps correctly to the simulated genome: GENOME_SPEC.csv)'
      echo '(09) QTL_min_max_GEBV.spec (expected minimum and maximum phenotypic values; will be corrected during parsing to account for the observed phenotypes; HEADER:[MIN_GEBV, MAX_GEBV]; tab-delimited)'
      echo '(10) Lperenne_genome.spec (Lolium perenne specific chromosome lengths in bp (@Ansari2016) and cM (@Pfeifer2013); hard-coded; HEADER:[Chrom, Mbp, cM]; tab-delimited)'
      echo '(11) ${OUTPREFIX}_${TIMESTAMP}/${OUTPREFIX}_g${nGen}.ped (biallelic plink genotype file in pedigree format; first 2 lines are comments; tab-delimited; |col1:family_ID|col2:within_family_ID|col3:within_family_ID_male|col4:within_family_ID_female|col5:sex|col6:phenotype_0to1_range|col7-col2*M:genotype_string_eg_{1 2}|)'
      echo '(12) ${OUTPREFIX}_${TIMESTAMP}/${OUTPREFIX}_g${nGen}.pheno (phenotype file; first 2 lines are comments; tab-delimited; |col1:family_ID|col2:individual_ID|col3:foreground_selection_phenotypes|col4:background_selection_phenotypes|)'
      echo '(13) ${OUTPREFIX}_${TIMESTAMP}/${OUTPREFIX}.map (improper loci information:ID and positional info are wrong and missing background selection loci; probably a bug in quantiNemo2; will be replaced by GENOME_SPEC.csv during parsing)'
      echo ""
      echo ""
      echo '###############'
      echo '### EXAMPLE ###'
      echo '###############'
      echo 'DIR=/data/Lolium'
      echo 'SRC_DIR=${DIR}/Softwares'
      echo 'rep=1'
      echo 'nIndividuals=1000'
      echo 'nLoci=1000'
      echo 'nQTL=10'
      echo 'nBGS=100'
      echo 'nAlleles=2'
      echo 'allele_eff_model=CHISQ'
      echo 'nGen=100'
      echo 'nPop=16'
      echo 'migration=0.00'
      echo 'selection=0.25'
      echo 'bg_selection=0.00'
      echo 'GRADIENT=1'
      echo 'OUTPREFIX=GPASim-${rep}_REP-${nQTL}_NQTL-${nBGS}_NBGS-${migration}_MIG-${selection}_SEL-${bg_selection}_BGSEL-${GRADIENT}_GRAD'
      echo 'OUTDIR=${DIR}/Quantitative_Genetics/LOLSIM2020/${OUTPREFIX}'
      echo 'mkdir ${OUTDIR}'
      echo 'time \'
      echo '${SRC_DIR}/genomic_prediction/src/GPASim_1.0_simulate.sh \'
      echo '      -s $SRC_DIR \'
      echo '      -o $OUTDIR \'
      echo '      -p $OUTPREFIX \'
      echo '      -n $nIndividuals \'
      echo '      -l $nLoci \'
      echo '      -q $nQTL \'
      echo '      -b $nBGS \'
      echo '      -a $nAlleles \'
      echo '      -r $allele_eff_model \'
      echo '      -t $nGen \'
      echo '      -P $nPop \'
      echo '      -m $migration \'
      echo '      -sq $selection \'
      echo '      -sb $bg_selection \'
      echo '      -g $GRADIENT'
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

### SET WORKING DIRECTORY
cd $OUTDIR

### INPUT PARAMETER CHECK
echo -e "--softwares-dir\n--output-dir\n--output-prefix\n--n-individuals\n--n-loci\n--n-QTL\n--n-BGS\n--n-alleles\n--alleles-dist\n--n-generations\n--n-populations\n--migration-rate\n--selection-intensity\n--selection-intensity\n--diffusion-gradient" > input_parameter_names.temp
echo -e "$SRC_DIR\n$OUTDIR\n$OUTPREFIX\n$nIndividuals\n$nLoci\n$nQTL\n$nBGS\n$nAlleles\n$allele_eff_model\n$nGen\n$nPop\n$migration\n$selection\n$bg_selection\n$GRADIENT" > input_parameter_values.temp
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

### SANITY CHECK OF THE NUMBERS OF LOCI, FOREGROUND SELECTION QTL AND BACKGROUND SELECTION QTL
if [ $(echo $nQTL + $nBGS | bc) -gt  $nLoci ]
then
  echo "More QTL than there are loci!"
  echo "Number of foreground selection QTL (n-QTL|-q) + number of background selection QTL (n-BGS|-b) must be less than or equal to the total number of loci (n-loci|-l)."
  echo "Number of foreground selection QTL: $nQTL"
  echo "Number of background selection QTL: $nBGS"
  echo "Total number of loci: $nLoci"
  exit 1
fi
if [ $nQTL -lt 1 ] || [ $nBGS -lt 1 ]
then
  echo "The number of foreground (n-QTL|-q) and background (n-BGS|-b) QTL must range from 1 to the number of loci minus 1."
  exit 1
fi

### PREPARE THE QUANTINEMO2 INITIALIZATION FILE
### Setup mating system
echo -e "mating_system 0" > ${OUTPREFIX}.ini                    #hermaphrodite random mating
### Demography
echo -e "patch_number $nPop" >> ${OUTPREFIX}.ini                # Number of populations or subpopulations to simulate: must have a natural number square-root
echo -e "generations $nGen" >> ${OUTPREFIX}.ini                 #Lolium has been introduced to Australia around 1880
echo -e "patch_capacity $nIndividuals" >> ${OUTPREFIX}.ini      #carrying capacity
echo -e "regulation_model_offspring 1" >> ${OUTPREFIX}.ini      #regulation of population size to carrying capacity (patch_capacity) via random culling of offsprings
echo -e "regulation_model_adults 1" >> ${OUTPREFIX}.ini         #regulation of population size to carrying capacity (patch_capacity) via random culling of adults
# echo -e "mating_nb_offspring_model 3" >> ${OUTPREFIX}.ini       #simple fecundity by rouding the number of offsrpings) with fecundity rate of...
echo -e "mating_nb_offspring_model 0" >> ${OUTPREFIX}.ini       #total number of offspring is set to carrying capacity
echo -e "mean_fecundity 1" >> ${OUTPREFIX}.ini                  #... 1 which means a constant population size which is just an approximation though I don't know how bad of an approximation it is
echo -e "dispersal_rate $migration" >> ${OUTPREFIX}.ini         #migration rate across adjacent patches [@Wang2016]
# echo -e "dispersal_model 2" >> ${OUTPREFIX}.ini                 #1-dimensional stepping-stone model
# echo -e "dispersal_border_model 2" >> ${OUTPREFIX}.ini          #migrants from the border gets lost for the 1D stepping-stone model
echo -e "dispersal_model 3" >> ${OUTPREFIX}.ini                 #2D stepping-stone
echo -e "dispersal_lattice_range 1" >> ${OUTPREFIX}.ini         #2D stepping-stone: dispersal range of 8 adjacent patches (vertical, horizontal and diagonal; m/8)
echo -e "dispersal_border_model 2" >> ${OUTPREFIX}.ini          #absorbing boundaries: migration beyond the border is lost
### dispersal_lattice_dims by defualt sets a square patches structure sqrt(nPop) x sqrt(nPop) and if sqrt(nPop) is not a whole number quantiNemo2 will return an error!
# n_rows_cols=$(echo "sqrt(($nPop))" | bc)
# echo -e "dispersal_lattice_dims ($n_rows_cols, $n_rows_cols)" >> ${OUTPREFIX}.ini        #structure of population in a sqrt(n) x sqrt(n) square grid
###### Setting the genetic map
Rscript ${SRC_DIR}/genomic_prediction/src/GPASim_1.1_build_genome.r \
  7 \
  97.7,151.5,63.3,119.2,89.1,115.2,113.7  \
  $nLoci
  ### Input:
    ### 7 chromosomes in Lolium
    ### Chromosome length in cM from Pfeifer, 2013 (The Perennial Ryegrass GenomeZipper: Targeted Use of Genome Resources for Comparative Grass Genomics)
    ### number of loci
  ### Output: genome_loci_*.temp
echo -e "quanti_genome { {1: $(cat genome_loci_1.temp)}" >> ${OUTPREFIX}.ini
echo -e "\t\t{2: $(cat genome_loci_2.temp)}" >> ${OUTPREFIX}.ini
echo -e "\t\t{3: $(cat genome_loci_3.temp)}" >> ${OUTPREFIX}.ini
echo -e "\t\t{4: $(cat genome_loci_4.temp)}" >> ${OUTPREFIX}.ini
echo -e "\t\t{5: $(cat genome_loci_5.temp)}" >> ${OUTPREFIX}.ini
echo -e "\t\t{6: $(cat genome_loci_6.temp)}" >> ${OUTPREFIX}.ini
echo -e "\t\t{7: $(cat genome_loci_7.temp)} }" >> ${OUTPREFIX}.ini
### Genotype configuration
echo -e "quanti_loci_1 $(echo $nLoci - $nBGS | bc)" >> ${OUTPREFIX}.ini                #[trait of interest: QTL-based] total number of loci: neutral and QTL
echo -e "quanti_loci_2 $nBGS" >> ${OUTPREFIX}.ini                #[background selection trait]total number of loci: neutral and QTL
echo -e "quanti_all $nAlleles" >> ${OUTPREFIX}.ini              #number of alleles per locus; and since we're simulating SNPs we want have a maximum of 5, i.e. A,T,C,G and DEL excluding N; but for simplicity biallelic is preferred
echo -e "quanti_nb_trait 2" >> ${OUTPREFIX}.ini                 #number of traits: 2 - one for the quantitative trait of interest and the other as background selection simulator
    ###### Build the quanti_allelic_file for the QTL and background selection (BGS) loci specs
    echo -e "# Quantitative Alleles Specifications File" > QTL.spec ### loci specifications for the trait of interest
    echo -e "[FILE_INFO]{" >> QTL.spec
    echo -e "  col_locus 1" >> QTL.spec
    echo -e "  col_allele 2" >> QTL.spec
    echo -e "  col_allelic_value 3" >> QTL.spec
    echo -e "  col_mut_freq 4" >> QTL.spec
    echo -e "  col_ini_freq 5" >> QTL.spec
    echo -e "}" >> QTL.spec
    echo -e "#locus\tallele\tvalue\tmut_freq\tini_freq" >> QTL.spec
    echo -e "# Quantitative Alleles Specifications File" > BGS.spec ### loci specifications for the background selection trait
    echo -e "[FILE_INFO]{" >> BGS.spec
    echo -e "  col_locus 1" >> BGS.spec
    echo -e "  col_allele 2" >> BGS.spec
    echo -e "  col_allelic_value 3" >> BGS.spec
    echo -e "  col_mut_freq 4" >> BGS.spec
    echo -e "  col_ini_freq 5" >> BGS.spec
    echo -e "}" >> BGS.spec
    echo -e "#locus\tallele\tvalue\tmut_freq\tini_freq" >> BGS.spec
    cp QTL.spec QTL_WT.spec
    cp BGS.spec BGS_WT.spec
    ###### Empirical allele count distribution in Lolium (Lolium2018 Illumina sequencing): script found in "misc/quantinemo2_geno_sim_00_QTL_build_testing.r" commented-out
    echo -e "ALLELE_COUNT_1,0.0" > Allele_counts_dist.spec
    # echo -e "ALLELE_COUNT_2,0.673749892530331" >> Allele_counts_dist.spec
    # echo -e "ALLELE_COUNT_3,0.257344297077753" >> Allele_counts_dist.spec
    # echo -e "ALLELE_COUNT_4,0.0660887235376553" >> Allele_counts_dist.spec
    # echo -e "ALLELE_COUNT_5,0.00281708685426013" >> Allele_counts_dist.spec
    echo -e "ALLELE_COUNT_2,1.0" >> Allele_counts_dist.spec
    echo -e "ALLELE_COUNT_3,0.0" >> Allele_counts_dist.spec
    echo -e "ALLELE_COUNT_4,0.0" >> Allele_counts_dist.spec
    echo -e "ALLELE_COUNT_5,0.0" >> Allele_counts_dist.spec ### input for the GPASim_01_build_loci_spec.r script
    # echo -e "ALLELE_COUNT_6,0.0" >> Allele_counts_dist.spec
    ### NOTE: Allele_counts_dist.spec should align with the "quanti_all" variable above!
    if [ $(cat Allele_counts_dist.spec | wc -l) -lt $nAlleles ] || [ $(head -n${nAlleles} Allele_counts_dist.spec | tail -n1 | cut -d, -f2) == 0.0 ] || [ $(head -n${nAlleles} Allele_counts_dist.spec | tail -n1 | cut -d, -f2) == 0 ]
    then
      echo "The hard-coded 'Allele_counts_dist.spec' file does not match the number of alleles set!"
      echo "Exiting now!"
      exit 1
    fi
    ###### Simulate the loci specifications for trait of interest (QTL + zero-effect) and background selection trait
    Rscript ${SRC_DIR}/genomic_prediction/src/GPASim_1.2_build_loci_spec.r \
      $nLoci \
      $nQTL \
      $nBGS \
      Allele_counts_dist.spec \
      $allele_eff_model \
      $nIndividuals
      ### Input:
        ### total number of loci in the genome
        ### number of QTL for the trait of interest
        ### number of QTL for the background selection trait
        ### frequency table of the number of alleles per locus
        ### probability distribution model to sample the allele effects from (the same for both trait of interest and the background selection trait)
        ### effective population size or in this case the number of individuals
      ### Output:
        ### (1) "BGS.spec"
        ### (2) "BGS_idx_map_to_genome.spec"
        ### (3) "QTL.spec" (NOTE: includes the QTL of the trait of interest and zero-effect loci and EXCLUDING the BGS loci)
        ### (4) "QTL_idx_map_to_genome.spec"
        ### (5) "QTL_min_max_GEBV.spec"
echo -e "quanti_allelic_file_1 QTL.spec" >> ${OUTPREFIX}.ini      #allelic file for the QTL (trait of interest) + zero-effect loci
echo -e "quanti_allelic_file_2 BGS.spec" >> ${OUTPREFIX}.ini      #alellic file for the background selection loci
echo -e "quanti_locus_index_1 {$(uniq QTL_idx_map_to_genome.spec)}" >> ${OUTPREFIX}.ini
echo -e "quanti_locus_index_2 {$(uniq BGS_idx_map_to_genome.spec)}" >> ${OUTPREFIX}.ini
echo -e "quanti_ini_allele_model 0" >> ${OUTPREFIX}.ini         #genotypes are set to be maximally polymorph while set to 1 for monomorph or fixed alleles per genotype
echo -e "quanti_mutation_rate 0" >> ${OUTPREFIX}.ini            #mutation rate (set to no mutation for simplicity now)
echo -e "quanti_mutation_model 0" >> ${OUTPREFIX}.ini           #random mutation model (see pages 56-60)
### Phenotype settings
echo -e "quanti_environmental_model 0" >> ${OUTPREFIX}.ini      #environmental variance constant set as quanti_heritability=50%
echo -e "quanti_heritability 0.50" >> ${OUTPREFIX}.ini          #narrow-sense heritability set to 50%
echo -e "quanti_selection_model 2" >> ${OUTPREFIX}.ini          #directional selection (Richards, 1959)
echo -e "quanti_dir_sel_min 0.0" >> ${OUTPREFIX}.ini          #minimum fitness
echo -e "quanti_dir_sel_max 1.0" >> ${OUTPREFIX}.ini          #maximum fitness
echo -e "quanti_dir_sel_growth_rate $selection" >> ${OUTPREFIX}.ini     #maximal slope fixed to 1
MIN_GEBV=$(tail -n1 QTL_min_max_GEBV.spec | cut -f1)                      #minimum genomic breeding value
MAX_GEBV=$(tail -n1 QTL_min_max_GEBV.spec | cut -f2)                      #maximum genomic breeding value
SEL_THRESH=$(echo "scale=7; $MIN_GEBV + (((($MAX_GEBV - $MIN_GEBV)) * ((1-$selection)) ))" | bc)                     #half-way between minimum and maximum GEBV
echo -e "quanti_dir_sel_max_growth $SEL_THRESH" >> ${OUTPREFIX}.ini         #maximal slope half-way between minimum and maximu possible GEBV
# echo -e "quanti_dir_sel_symmetry 1" >> ${OUTPREFIX}.ini                   #symmetric logistic selection curve
### summary Statistics Output
# echo -e "stat_log_time $nGen" >> ${OUTPREFIX}.ini               #output summary statistics for the last generation only
# echo -e "stat {q.adlt.fst" >> ${OUTPREFIX}.ini
# echo -e "q.adlt.fst.wc_pair" >> ${OUTPREFIX}.ini
# echo -e "q.adlt.R2}" >> ${OUTPREFIX}.ini
### Save genotype and phenotype data
echo -e "quanti_genot_logtime $nGen" >> ${OUTPREFIX}.ini        #write genotype output for the final generation only
# echo -e "quanti_save_genotype 1" >> ${OUTPREFIX}.ini            #output in FSTAT format: main genotype block after the loci ID: n x l --> number of individuals x number of loci: [1:5][1:5] --> genotype format, e.g. 25 for a heterozygote containing the 2nd allele and the 5th allele
echo -e "quanti_save_genotype 5" >> ${OUTPREFIX}.ini            #output in PLINK format: ped and map files
echo -e "quanti_genot_filename $OUTPREFIX" >> ${OUTPREFIX}.ini
echo -e "quanti_phenot_logtime $nGen" >> ${OUTPREFIX}.ini       #write phenotype output for the final generation only
# echo -e "quanti_phenot_logtime 10" >> ${OUTPREFIX}.ini        #test
# echo -e "quanti_save_phenotype 1" >> ${OUTPREFIX}.ini
echo -e "quanti_save_phenotype 0" >> ${OUTPREFIX}.ini           #do not output a separate phenotype file since it will be included in the plink output files
echo -e "quanti_phenot_filename $OUTPREFIX" >> ${OUTPREFIX}.ini

### WRITE OUT CHROMOSOME LENGTHS in bp (@Ansari2016) and cM (@Pfeifer2013)
echo -e "Chrom\tMbp\tcM" > Lperenne_genome.spec
echo -e "1\t470.30\t97.7" >> Lperenne_genome.spec
echo -e "2\t429.90\t151.5" >> Lperenne_genome.spec
echo -e "3\t406.56\t63.3" >> Lperenne_genome.spec
echo -e "4\t369.03\t119.2" >> Lperenne_genome.spec
echo -e "5\t339.41\t89.1" >> Lperenne_genome.spec
echo -e "6\t322.62\t115.2" >> Lperenne_genome.spec
echo -e "7\t284.07\t113.7" >> Lperenne_genome.spec

### GENERATE LOCI SPECIFICATIONS FLE: CHROM, r (cM), POS (bp) and the QTL spec file for the trait of interest (excluding the allele effect of the BGS QTL)
ls genome_loci_*.temp | sed s/.temp//g | rev | cut -d'_' -f1 > chrom_id.temp
cat genome_loci_*.temp > r_pos.temp
paste -d' ' chrom_id.temp r_pos.temp > chrom_r.temp ### locus positions in cM
Rscript ${SRC_DIR}/genomic_prediction/src/GPASim_1.3_genome_QTL_spec_parsing.r \
  chrom_r.temp \
  Lperenne_genome.spec \
  QTL.spec \
  QTL_idx_map_to_genome.spec
### Input:
### (1) locus positions per chromosome (NO HEADER;space-delimited; 1 chrom/line)
### (2) chromosome specifications (HEADER: Chrom, Mbp, cM; tab-delimited)
### (3) QTL specifications file (HEADER: metadata-see above; tab-delimited; excludes BGS QTL)
### (4) QTL indices that maps to the genome (NO HEADER; tab-delimited; 1 column)
### Output:
### (1) "GENOME_SPEC.csv" (HEADER: CHROM, r, POS)
### (2) "QTL_SPEC.csv" (HEADER: CHROM, POS, ALLELE, EFFECT)

### GENERATE INTITAL GENOTYPES
if [ $GRADIENT == 0 ]
then
  echo -e "Using NULL gradient of allele effects across the landscape."
elif [ $GRADIENT == 1 ]
then
  echo -e "Using a gradient of allele effects across the landscape."
  echo -e "Alleles diffusing from the first row of populations."
  ### WILDTYPE
  sed s/BGS.spec/BGS_WT.spec/g ${OUTPREFIX}.ini > ${OUTPREFIX}_WILDTYPE.ini
  sed -i s/QTL.spec/QTL_WT.spec/g ${OUTPREFIX}_WILDTYPE.ini
  sed -i s/"generations $nGen"/"generations 1"/g ${OUTPREFIX}_WILDTYPE.ini
  sed -i s/"quanti_genot_logtime $nGen"/"quanti_genot_logtime 1"/g ${OUTPREFIX}_WILDTYPE.ini
  sed -i s/"quanti_phenot_logtime $nGen"/"quanti_phenot_logtime 1"/g ${OUTPREFIX}_WILDTYPE.ini
  sed -i s/"quanti_save_genotype 5"/"quanti_save_genotype 1"/g ${OUTPREFIX}_WILDTYPE.ini ### output dat genotype instead of ped files
  ${SRC_DIR}/quantinemo_linux/quantinemo ${OUTPREFIX}_WILDTYPE.ini
  cp ${OUTPREFIX}_*/*.dat WILDTYPE.dat ### P(alleles with effects) = 0.00
  rm -R ${OUTPREFIX}_*/
  ### UNIFORM
  sed s/"generations $nGen"/"generations 1"/g ${OUTPREFIX}.ini > ${OUTPREFIX}_UNIFORM.ini
  sed -i s/"quanti_genot_logtime $nGen"/"quanti_genot_logtime 1"/g ${OUTPREFIX}_UNIFORM.ini
  sed -i s/"quanti_phenot_logtime $nGen"/"quanti_phenot_logtime 1"/g ${OUTPREFIX}_UNIFORM.ini
  sed -i s/"quanti_save_genotype 5"/"quanti_save_genotype 1"/g ${OUTPREFIX}_UNIFORM.ini ### output dat genotype instead of ped files
  ${SRC_DIR}/quantinemo_linux/quantinemo ${OUTPREFIX}_UNIFORM.ini
  cp ${OUTPREFIX}_*/*.dat UNIFORM.dat ### P(allele_max) = 1/2N
  rm -R ${OUTPREFIX}_*/
  ### MERGE
  NCOL=$(echo "sqrt ($nPop)" | bc)
  ROW1_NLINES=$(echo "(1 + ${nLoci}) + (${NCOL} * ${nIndividuals})" | bc)
  TAILROW1_NLINES=$(echo "${ROW1_NLINES} + 1" | bc)
  head -n${ROW1_NLINES} UNIFORM.dat > GENOTYPES_INI.dat           ### header +row 1 diffusers
  tail -n+${TAILROW1_NLINES} WILDTYPE.dat >> GENOTYPES_INI.dat   ### wildtype populations
  ### APPEND TO THE INI FILE
  echo -e "quanti_ini_genotypes GENOTYPES_INI.dat" >> ${OUTPREFIX}.ini
elif [ $GRADIENT == 2 ]
then
  echo -e "Using a gradient of allele effects across the landscape."
  echo -e "Alleles diffusing from the first and last rows of populations."
  ### WILDTYPE
  sed s/BGS.spec/BGS_WT.spec/g ${OUTPREFIX}.ini > ${OUTPREFIX}_WILDTYPE.ini
  sed -i s/QTL.spec/QTL_WT.spec/g ${OUTPREFIX}_WILDTYPE.ini
  sed -i s/"generations $nGen"/"generations 1"/g ${OUTPREFIX}_WILDTYPE.ini
  sed -i s/"quanti_genot_logtime $nGen"/"quanti_genot_logtime 1"/g ${OUTPREFIX}_WILDTYPE.ini
  sed -i s/"quanti_phenot_logtime $nGen"/"quanti_phenot_logtime 1"/g ${OUTPREFIX}_WILDTYPE.ini
  sed -i s/"quanti_save_genotype 5"/"quanti_save_genotype 1"/g ${OUTPREFIX}_WILDTYPE.ini ### output dat genotype instead of ped files
  ${SRC_DIR}/quantinemo_linux/quantinemo ${OUTPREFIX}_WILDTYPE.ini
  cp ${OUTPREFIX}_*/*.dat WILDTYPE.dat ### P(alleles with effects) = 0.00
  rm -R ${OUTPREFIX}_*/
  ### UNIFORM
  sed s/"generations $nGen"/"generations 1"/g ${OUTPREFIX}.ini > ${OUTPREFIX}_UNIFORM.ini
  sed -i s/"quanti_genot_logtime $nGen"/"quanti_genot_logtime 1"/g ${OUTPREFIX}_UNIFORM.ini
  sed -i s/"quanti_phenot_logtime $nGen"/"quanti_phenot_logtime 1"/g ${OUTPREFIX}_UNIFORM.ini
  sed -i s/"quanti_save_genotype 5"/"quanti_save_genotype 1"/g ${OUTPREFIX}_UNIFORM.ini ### output dat genotype instead of ped files
  ${SRC_DIR}/quantinemo_linux/quantinemo ${OUTPREFIX}_UNIFORM.ini
  cp ${OUTPREFIX}_*/*.dat UNIFORM.dat ### P(allele_max) = 1/2N
  rm -R ${OUTPREFIX}_*/
  ### MERGE
  NCOL=$(echo "sqrt ($nPop)" | bc)
  ROW1_NLINES=$(echo "(1 + ${nLoci}) + (${NCOL} * ${nIndividuals})" | bc)
  TAILROW1_NLINES=$(echo "${ROW1_NLINES} + 1" | bc)
  TOTAL_NLINES=$(echo "(1 + ${nLoci}) + (${nPop} * ${nIndividuals})" | bc)
  SANDWICH_NLINES=$(echo "$TOTAL_NLINES - $ROW1_NLINES - (${NCOL} * ${nIndividuals})" | bc)
  ROW2_NLINES=$(echo "${NCOL} * ${nIndividuals}" | bc)
  head -n${ROW1_NLINES} UNIFORM.dat > GENOTYPES_INI.dat                                       ### header +row1 diffusers
  tail -n+${TAILROW1_NLINES} WILDTYPE.dat | head -n${SANDWICH_NLINES}  >> GENOTYPES_INI.dat  ### wild type populations
  tail -n${ROW2_NLINES} UNIFORM.dat >> GENOTYPES_INI.dat                                     ### last row diffusers
  ### APPEND TO THE INI FILE
  echo -e "quanti_ini_genotypes GENOTYPES_INI.dat" >> ${OUTPREFIX}.ini
else
  echo -e "Sorry ${GRADIENT} is not a valid gradient input."
  echo -e "   - use 0 for none (i.e. alleles with effects randomly scattered across the landscape)"
  echo -e "   - use 1 for alleles with effects diffusing from the first row of populations"
fi
### clean-up
rm *.temp QTL_WT.spec BGS_WT.spec

### EXECUTE
${SRC_DIR}/quantinemo_linux/quantinemo ${OUTPREFIX}.ini
rm GENOTYPES_INI.dat UNIFORM.dat WILDTYPE.dat || echo -e "Using NULL gradient of allele effects across the landscape."
