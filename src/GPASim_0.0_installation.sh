#!/bin/bash

###################################################################
### PREPARE SPARTAN FOR QUANTINEMO SIMULATION AND DATA ANALYSIS ###
###################################################################

### SAMPLE USAGE
# projectID=punim0543
# username=$USER
# DIR=/data/cephfs/${projectID}/${username}
# # DIR=/data/Lolium/Quantitative_Genetics/LOLSIM_2019
# ./GPWASim_00_installation.sh $DIR

### INPUT
DIR=$1

### create installation directory
mkdir ${DIR}
mkdir ${DIR}/Softwares
mkdir ${DIR}/Scripts
mkdir ${DIR}/Output
echo -e "Your installation directory is in $DIR"

### prepare quantiNemo2, genomic_prediction.git, popoolation2, and julia1 libraries
###### download quantiNemo2 if not yet installed
if [ $(ls ${DIR}/Softwares/quantinemo_linux | wc -l) -eq 0 ]
then
  echo "########################################"
  echo "Downloading QuantiNemo2"
  echo "########################################"
  cd ${DIR}/Softwares
  wget https://www2.unil.ch/popgen/softwares/quantinemo/files/quantinemo_linux.zip
  unzip quantinemo_linux.zip
  rm quantinemo_linux.zip
  echo "########################################"
fi
###### clone the genomic_prediction git repo not yet installed
if [ $(ls ${DIR}/Softwares/genomic_prediction | wc -l) -eq 0 ]
then
  echo "########################################"
  echo "Cloning genomic_prediction.git repository"
  echo "########################################"
  cd ${DIR}/Softwares
  git clone https://gitlab.com/jeffersonfparil/genomic_prediction.git
  echo "########################################"
fi
###### download plink if not yet installed
if [ $(ls ${DIR}/Softwares/plink | wc -l) -eq 0 ]
then
  echo "########################################"
  echo "Downloading PLINK"
  echo "########################################"
  cd ${DIR}/Softwares
  wget http://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20200121.zip
  unzip plink_linux_x86_64_20200121.zip
  rm plink_linux_x86_64_20200121.zip
  echo "########################################"
fi
###### download emmax if not yet installed
if [ $(ls ${DIR}/Softwares/emmax-intel64 | wc -l) -eq 0 ]
then
  echo "########################################"
  echo "Downloading EMMAX"
  echo "########################################"
  cd ${DIR}/Softwares
  wget http://csg.sph.umich.edu//kang/emmax/download/emmax-intel-binary-20120210.tar.gz
  tar -xvzf emmax-intel-binary-20120210.tar.gz
  rm emmax-intel-binary-20120210.tar.gz
  echo "########################################"
fi
###### download gemma if not yet installed
if [ $(ls ${DIR}/Softwares/gemma-0.98.1-linux-static | wc -l) -eq 0 ]
then
  echo "########################################"
  echo "Downloading GEMMA"
  echo "########################################"
  cd ${DIR}/Softwares
  wget https://github.com/genetics-statistics/GEMMA/releases/download/0.98.1/gemma-0.98.1-linux-static.gz
  gunzip gemma-0.98.1-linux-static.gz
  chmod +x gemma-0.98.1-linux-static
  echo "########################################"
fi
###### download gcta if not yet installed
if [ $(ls ${DIR}/Softwares/gcta_1.93.0beta/gcta64 | wc -l) -eq 0 ]
then
  echo "########################################"
  echo "Downloading GCTA"
  echo "########################################"
  cd ${DIR}/Softwares
  wget https://cnsgenomics.com/software/gcta/bin/gcta_1.93.0beta.zip
  unzip gcta_1.93.0beta.zip
  rm gcta_1.93.0beta.zip
  echo "########################################"
fi
###### download popoolation2 if not yet installed
if [ $(ls ${DIR}/Softwares/popoolation2_1201 | wc -l) -eq 0 ]
then
  echo "########################################"
  echo "Downloading Popoolation2"
  echo "########################################"
  cd ${DIR}/Softwares
  wget https://sourceforge.net/projects/popoolation2/files/popoolation2_1201.zip
  unzip popoolation2_1201.zip -d popoolation2_1201
  rm popoolation2_1201.zip
  echo "########################################"
fi
###### download NPSTAT if not yet installed
if [ $(ls ${DIR}/Softwares/npstat | wc -l) -eq 0 ]
then
  echo "########################################"
  echo "Downloading NPSTAT"
  echo "########################################"
  cd ${DIR}/Softwares
  git clone https://github.com/lucaferretti/npstat.git
  cd npstat/
  make || module load GSL/2.5-intel-2018.u4; make || sudo apt-get install libgsl-dev; make
  echo "########################################"
fi

#### NOTE: Spartan interactive
# sinteractive --partition=snowy -c 32 --mem=120000 --time=0-8:00:00
# module load parallel/20181222-spartan_gcc-6.2.0.lua
# module load Julia/1.1.1-spartan_gcc-6.2.0.lua
# module load R/3.5.2-GCC-6.2.0
# module load GSL/2.5-intel-2018.u4
#### NOTE: download julia 1.1 if not yet installed
# if [ $(ls ${DIR}/Softwares/julia-1.1.1 | wc -l) -eq 0 ]
# then
#   echo "########################################"
#   echo "Downloading Julia 1.1.1"
#   echo "########################################"
#   cd ${DIR}/Softwares
#   wget https://julialang-s3.julialang.org/bin/linux/x64/1.1/julia-1.1.1-linux-x86_64.tar.gz
#   tar -xvzf julia-1.1.1-linux-x86_64.tar.gz
#   rm julia-1.1.1-linux-x86_64.tar.gz
#   echo "########################################"
# fi
##### NOTE: INSTALL R PACKAGES MANUALLY!!! (FOR julia's RCall and ggmix)
# module load R/3.6.1-GCC-6.2.0
# R
# install.packages("glmnet")
# install.packages("RColorBrewer")
# install.packages("abc")
# q(save='no')
# # cd ~/R/
# # Rlib=$(ls | tail -n1)
# # cd $Rlib
# # Rver=$(ls | tail -n1)
# # cd $Rver
# # git clone https://github.com/cran/agricolae.git
# #### NOTE: INSTALL PYTHON'S MATPLOTLIB MANUALLY!!! (FOR julia's PyCall and Plots' plotting)
# # module load Python/2.7.13-GCC-6.2.0-bare
# # pip install matplotlib --user
# ### NOTE: INSTALL JULIA LIBRARIES MANUALLY
# module load R/3.6.1-GCC-6.2.0
# module load Julia/1.1.1-spartan_gcc-6.2.0.lua
# julia
# using Pkg
# Pkg.add([
# "CSV",
# "Distributions",
# "DataFrames",
# "DelimitedFiles",
# "LinearAlgebra",
# "Optim",
# "UnicodePlots",
# "ColorBrewer",
# "ProgressMeter",
# "Statistics",
# "StatsBase",
# "Distributed",
# "SharedArrays",
# "GeneticVariation",
# "RCall"
# ])
# using Pkg
# Pkg.add(PackageSpec(url="https://github.com/jeffersonfparil/GWAlpha.jl.git", rev="master"))
# exit()
# ### NOTE: may need to manualy download tar.gz dependencies from github and scp into spartan into the directories specified by the error messages at Pkg.build("GWAlpha")
# ### e.g.
# ### (1) https://github.com/JuliaBinaryWrappers/Arpack_jll.jl/releases/download/Arpack-v3.5.0+2/Arpack.v3.5.0.x86_64-linux-gnu-libgfortran3.tar.gz
# ### (2) https://github.com/staticfloat/RmathBuilder/releases/download/v0.2.0-1/libRmath.x86_64-linux-gnu.tar.gz
# ### (3) https://github.com/JuliaPackaging/Yggdrasil/releases/download/OpenSpecFun-v0.5.3+0/OpenSpecFun.v0.5.3.x86_64-linux-gnu-gcc4.tar.gz
# ### (4) https://github.com/bicycle1885/ZlibBuilder/releases/download/v1.0.4/Zlib.v1.2.11.x86_64-linux-gnu.tar.gz