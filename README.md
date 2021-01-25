# GPAS Landscape Simulation

This repository contains the code for simulating the landscapes used in the manuscript: [**Optimising sampling design for the genomic analysis of quantitative traits in natural populations**](https://www.authorea.com/users/369690/articles/488504-optimising-sampling-design-for-the-genomic-analysis-of-quantitative-traits-in-natural-populations).

## Hardware requirement:
- a high-performance computing cluster (≥ 100 GB of RAM per node with 32 compute cores)

## Software requirement:
- Linux operating system running a Debian-based package manager (preferrably Ubuntu)

## Source code description:
- `GPASim_0.0_installation.sh` installs the required software packages with some notes at the bottom of the script for manual installation of R and Julia packages.
- `GPASim_1.0_simulate.sh` and associated scripts simulate the landscapes under various genetic and demographic scenarios.
- `GPASim_2.0_parse.sh` and associated scripts parse the output of the previous steps.
- `GPASim_3.0_across_population_sampling.sh` and associated scripts perform the stratified sampling strategy.
- `GPASim_4.0_filter_build_covariates.sh` and associated scripts build the genetic relatedness and Fst matrices for GPAS.
- `GPASim_5.0_GPAS.sh` and associated scripts perform GPAS using Indi-seq and Pool-seq genotype data, as well as assess the performance of GPAS, i.e. the accuracies of QTL detection and genomic prediction.
- `GPASim_6.0_merging_and_abc_optim.sh` and associated scripts merges and summarises the output of GPAS (the ABC optimisation component was dropped from the final methodology)
- `GPASim_7.0_slurmer.sh` is the [Spartan-specific](https://dashboard.hpc.unimelb.edu.au/) script to generate multiple [Slurm](https://slurm.schedmd.com/documentation.html) scripts to run on indivudal nodes of a high-performance computing cluster. The input parameters used in the paper are harcoded in lines 85 to 97:
  + `rep` - number of replications
  + `nQTL` - number of QTL
  + `nBGS` - number of background QTL (ignore this since this is constant and does not affect the trait of interest)
  + `migration` - migration rate
  + `selection` - selection intensity for the trait of interest
  + `bg_selection` - background selection intesity (ignore this since this is constant and does not affect the trait of interest)
  + `GRADIENT` - causal allele diffusion gradient (code: 0 for uniform distribution across the landscape; 1 for the top rows of populations as the origin of causal alleles; and 2 for the top and bottom rows of populations as the origin of the causal alleles)
