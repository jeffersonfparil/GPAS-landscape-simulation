#!/bin/bash

### Input parameter description
# ( 1) ${SRC_DIR}           directory where qunatinemo2 folder, plink, emmax, gemma, gcta, and genomic_prediction folder are located
# ( 3) ${rep}               replicate number
# ( 4) ${nIndividuals}      number of individuals to simulate
# ( 5) ${nLoci}             number of loci to simulate (neutral + QTL)
# ( 6) ${nQTL}              number of QTL among the all the loci simulated (should range from 1 to nLoci-1)
# ( 7) ${nBGS}              number of background selection loci (should range from 1 to nLoci-1; conditional on nQTL where nQTL + nBGS <= nLoci)
# ( 8) ${nAlleles}          number of alleles per loci, e.g. 5 for A,T,C,G, and DEL (NOTE: We're now assuming multi-allelic SNP loci are just noise so we're simulating only biallelic loci!)
# ( 9) ${allele_eff_model}  probability distribution to sample the allele effects from (i.e. CHISQ, NORM, & UNIF)
# (10) ${nGen}              number of generations to simulate (approximating the number of generations since the introduction of Lolium in Australia in the 1880s minus the time when extensive herbicide application started)
# (11) ${nPop}              number of populations/subpopulations to simulate (NOTE: must have a natural number square-root)
# (12) ${migration}         migration rate across the populations (surrently using the 1D stepping stone model see line 117)
# (13) ${selection}         proportion of individual selected based on trait of interest (nQTL) defined as proportion of the most fit individuals selected, where fitness is defined using the generalized logistic curve (Richards, 1959): ranges from 0.00 to 1.00; lower s means more intense selection
# (14) ${bg_selection}      proportion of individual selected based on background selection trait (nBGS) defined as proportion of the most fit individuals selected, where fitness is defined using the generalized logistic curve (Richards, 1959): ranges from 0.00 to 1.00; lower s means more intense selection
# (15) ${GRADIENT}          uniformly distributed non-wildtype alleles
# (16) ${NTHREADS}          number of computing cores or threads to run in parallel for some of the modules
# (17) ${NPOOLS}            number of pools to genotype per population (NOTE: for across population sample 1 population:1 pool)
# (18) ${NLIB}              maximum number of individually barcoded libraries for sequencing (Indi-seq and Pool-seq) (NOTE: May need to have a minimum limit at 1,000 to keep GCTA from failing due to small population size leading to REML convergence failure)
# (19) ${NSEQ_MAX}          maximum number of multiplexed sequencing experiments (limiting factor in sampling strategy optimization)
# (20) SPARTAN_ACCOUNT      Spartan account name
# (21) SPARTAN_PARTITION    Sparton node partition
# (22) SPARTAN_CORES        Partition-specific number of cores
# (23) SPARTAN_RAM          Partition-specific RAM

### Hardcoded input parameters
DIR=/data/scratch/projects/punim0543/jparil
# DIR=/data/cephfs/punim0543/jparil
SRC_DIR=${DIR}/Softwares
nIndividuals=1000
nLoci=10000
nAlleles=2
allele_eff_model=CHISQ
nGen=200 ### boosted up from 150 to 200
nPop=100
NPOOLS=5
NLIB=384 ### reduced to an even 384 samples or 4 x 96-well plates
NSEQ_MAX=10
SPARTAN_PARTITION=snowy

### Spartan account name
SPARTAN_ACCOUNT=punim0543
### maximum computing time limit in days-hours:minutes:seconds
SPARTAN_TIMELIMIT="2-0:0:00"
### Spartan node partition to use and the corrensponding comfortably nearish to maximum available cores and RAM
if [ $SPARTAN_PARTITION == cloud ]
then
  SPARTAN_CORES=12
  SPARTAN_RAM=95
elif [ $SPARTAN_PARTITION == longcloud ]
then
  SPARTAN_CORES=12
  SPARTAN_RAM=95
elif [ $SPARTAN_PARTITION == physical ]
then
  SPARTAN_CORES=32
  SPARTAN_RAM=250
elif [ $SPARTAN_PARTITION == bigmem ]
then
  SPARTAN_CORES=36
  SPARTAN_RAM=1535
elif [ $SPARTAN_PARTITION == snowy ]
then
  SPARTAN_CORES=32
  SPARTAN_RAM=100 ### reduced to 110 or less according to Sean Crosby 20200806
elif [ $SPARTAN_PARTITION == mig ]
then
  SPARTAN_CORES=32
  SPARTAN_RAM=120
  SPARTAN_ACCOUNT=punim0594
else
  echo -e "Unrecognised Spartan partition."
  echo -e "Select from:"
  echo -e "\t - cloud"
  echo -e "\t - longcloud"
  echo -e "\t - physical"
  echo -e "\t - bigmem"
  echo -e "\t - snowy"
  echo -e "\t - mig"
  exit
fi

### Iterate across these simulation parameters ### TEST 1 SET FIRST! 2020-04-25
for rep in 2 3 4 5
do
for nQTL in 100 10 50
do
for nBGS in 100
do
for migration in 0.001 0.0001 0.01 
do
for selection in 0.05 0.10 0.50
do
for bg_selection in 0.10
do
for GRADIENT in 1 2 0
do
### build the slurm script
echo -e '#!/bin/bash' > GPASim_${rep}REP_${nQTL}NQTL_${nBGS}NBGS_${migration}MIG_${selection}SEL_${bg_selection}BGSEL_${GRADIENT}GRAD.slurm
echo -e "
# Partition for the job:
#SBATCH --partition=${SPARTAN_PARTITION}
# Multithreaded (SMP) job: must run on one node and the cloud partition
#SBATCH --nodes=1
# The name of the job:
#SBATCH --job-name=GPASim_${rep}REP_${nQTL}NQTL_${nBGS}NBGS_${migration}MIG_${selection}SEL_${bg_selection}BGSEL_${GRADIENT}GRAD
# The project ID which this job should run under:
#SBATCH --account=${SPARTAN_ACCOUNT}
# Maximum number of tasks/CPU cores used by the job:
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=${SPARTAN_CORES}
# The amount of memory in megabytes per process in the job:
#SBATCH --mem=${SPARTAN_RAM}GB
# The maximum running time of the job in days-hours:mins:sec
#SBATCH --time=${SPARTAN_TIMELIMIT}
# Send yourself an email when the job:
# aborts abnormally (fails)
#SBATCH --mail-type=FAIL
# begins
# #SBATCH --mail-type=BEGIN
# ends successfully
# #SBATCH --mail-type=END
# Use this email address:
#SBATCH --mail-user=jparil@student.unimelb.edu.au

### load modules
module load parallel/20181222-spartan_gcc-6.2.0.lua
module load Julia/1.1.1-spartan_gcc-6.2.0.lua
module load R/3.6.1-GCC-6.2.0

### Execute
${SRC_DIR}/genomic_prediction/src/GPASim_7.1_execute_per_node.sh \
    ${DIR} \
    ${SRC_DIR} \
    ${rep} \
    ${nIndividuals} \
    ${nLoci} \
    ${nQTL} \
    ${nBGS} \
    ${nAlleles} \
    ${allele_eff_model} \
    ${nGen} \
    ${nPop} \
    ${migration} \
    ${selection} \
    ${bg_selection} \
    ${GRADIENT} \
    ${SPARTAN_CORES} \
    ${NPOOLS} \
    ${NLIB} \
    ${NSEQ_MAX}
" >> GPASim_${rep}REP_${nQTL}NQTL_${nBGS}NBGS_${migration}MIG_${selection}SEL_${bg_selection}BGSEL_${GRADIENT}GRAD.slurm
### submit into queue
sbatch GPASim_${rep}REP_${nQTL}NQTL_${nBGS}NBGS_${migration}MIG_${selection}SEL_${bg_selection}BGSEL_${GRADIENT}GRAD.slurm
done
done
done
done
done
done
done

# ##################
# ### MONITORING ###
# ##################
# ### alias whatup='squeue -u jparil | tail; cd /data/scratch/projects/punim0543/jparil; echo RUNNING: $(squeue -u jparil | grep "jparil  R" | wc -l); echo PENDING: $(squeue -u jparil | grep "jparil PD" | wc -l); echo FINISHED: $(if [ $(tail -n1 slurm-* | grep -B 1 "FINALLY DONE" | sort | uniq | wc -l) -eq 2 ]; then echo $(tail -n1 slurm-* | grep -B 1 "FINALLY DONE" | sort | uniq | tail -n+2 | cut -d" " -f2 | wc -l) ; else echo $(tail -n1 slurm-* | grep -B 1 "FINALLY DONE" | sort | uniq | tail -n+3 | cut -d" " -f2 | wc -l); fi); echo ""'
# if [ $(tail -n1 slurm-* | grep -B 1 "FINALLY DONE" | sort | uniq | wc -l) -eq 2 ]
# then
# tail -n1 slurm-* | grep -B 1 "FINALLY DONE" | sort | uniq | tail -n+2 | cut -d' ' -f2 > FINISHED_slurm.list
# else
# tail -n1 slurm-* | grep -B 1 "FINALLY DONE" | sort | uniq | tail -n+3 | cut -d' ' -f2 > FINISHED_slurm.list
# fi
# while IFS= read -r f
# do
#   DIR=$(grep "GPASim_" $f | grep "GRAD.ini" | grep "^Reading settings file" | head -n1 | cut -d"'" -f2 | sed 's/.ini//g')
#   if [ $(echo $f) == $(head -n1 FINISHED_slurm.list) ]
#   then
#   echo $DIR > FINISHED_dir.list
#   else
#   echo $DIR >> FINISHED_dir.list
#   fi
# done < FINISHED_slurm.list


# ### check before moving slurm.out files into their respective folder
# cat FINISHED_*
# ### move slurm.out files into their respective folder
# for i in $(seq 1 $(cat FINISHED_slurm.list | wc -l))
# do
#   slurm=$(head -n${i} FINISHED_slurm.list | tail -n1)
#   dir=$(head -n${i} FINISHED_dir.list | tail -n1)
#   mv $slurm $dir
#   mv ${dir}.slurm $dir
# done


# ### scp finished simulation directories from and into ssh_uni
# ssh_uni
# cd /data/Lolium/Quantitative_Genetics/LOLSIM_2020/
# scp -o PubkeyAuthentication=no jparil@spartan.hpc.unimelb.edu.au:/data/scratch/projects/punim0543/jparil/FINISHED_dir.list .
# for dir in $(cat FINISHED_dir.list)
# do
#   echo $dir
#   ### delete me after scp-ing I have a password!!!
#   sshpass -p ####### scp -r -o PubkeyAuthentication=no jparil@spartan.hpc.unimelb.edu.au:/data/scratch/projects/punim0543/jparil/${dir} .
# done


# ### BACK ON SPARTAN
# ### clean-up
# for dir in $(cat FINISHED_dir.list)
# do
#   echo $dir
#   rm -R $dir
# done
