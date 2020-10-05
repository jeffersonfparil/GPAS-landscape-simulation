#!/bin/bash
### WRITE POPULATION_COMBINATORICS.JL THAT ADAPTS TO $NSEQ
### WILL BE CALLED BY "GPASIM_6.2_ABC_INPUT.JL" FOR PARALLEL EXECUTION

NSEQ=$1     ### number of sequencing experiments within budget

### adaptive exhaustive combinatorics (grows exponentially with nPops! Will be very inefficient!)
echo "using Distributed" > population_combinatorics.jl
echo "using DelimitedFiles" >> population_combinatorics.jl
echo "using DataFrames" >> population_combinatorics.jl
echo "using CSV" >> population_combinatorics.jl
echo "function population_combinatorics(populations_list)" >> population_combinatorics.jl
echo "n = length(populations_list)" >> population_combinatorics.jl
# echo "ncombin = factorial(BigInt(n)) / ( factorial("${NSEQ}")*factorial(BigInt(n-"${NSEQ}")) )" >> population_combinatorics.jl
# echo 'println(string("Iterating across ", ncombin, " combinations..." ))' >> population_combinatorics.jl
echo "counter = [1]" >> population_combinatorics.jl
### prepare iterator names
### max NSEQ=3,844=(26+26+10^2)
touch iterators_list.txt
for c1 in {A..Z} {a..z} {0..9}
do
    for c2 in {A..Z} {a..z} {0..9}
    do
      printf "%s\n" "$c1$c2" >> iterators_list.txt
    done
done
# iterator_old=$(cut -d' ' -f1 <(echo {a..z}))
iterator_old=$(head -n1 iterators_list.txt | tail -n1)
### define sequence jumps which will restrict the number of combinations into a less insame amount, as ceil(n/40)
echo "@time @sync @distributed for "${iterator_old}" in 1:Int(ceil(n/40)):(n-1)" >> population_combinatorics.jl
# echo "@time @sync @distributed for "${iterator_old}" in 1:(n-1)" >> population_combinatorics.jl
for i in $(seq 2 $NSEQ)
do
  # iterator_new=$(cut -d' ' -f$i <(echo {a..z}))
  iterator_new=$(head -n$i iterators_list.txt | tail -n1)
  ### define sequence jumps which will restrict the number of combinations into a less insame amount, as ceil(n/40)
  echo "for "${iterator_new}" in ("${iterator_old}"+1):Int(ceil(n/40)):(n-("${NSEQ}"-"${i}"))" >> population_combinatorics.jl
  iterator_old=${iterator_new}
done
echo "all_pops = vcat([" >> population_combinatorics.jl
for i in $(seq 1 $(echo $NSEQ - 1 | bc))
do
  # iterator=$(cut -d' ' -f$i <(echo {a..z}))
  iterator=$(head -n$i iterators_list.txt | tail -n1)
  echo 'split(populations_list['${iterator}'], ":"),' >> population_combinatorics.jl
done
# iterator=$(cut -d' ' -f$NSEQ <(echo {a..z}))
iterator=$(head -n$NSEQ iterators_list.txt | tail -n1)
echo 'split(populations_list['${iterator}'], ":")' >> population_combinatorics.jl
echo "]...)" >> population_combinatorics.jl
echo "if length(all_pops) == length(unique(all_pops))" >> population_combinatorics.jl
### output filename
printf "fname_out=string(\"POPULATION_COMBINATIONS_\", " >> population_combinatorics.jl
for i in $(seq 1 $(echo $NSEQ - 1 | bc))
do
  # iterator=$(cut -d' ' -f$i <(echo {a..z}))
  iterator=$(head -n$i iterators_list.txt | tail -n1)
  printf "${iterator}, \"_\", " >> population_combinatorics.jl
done
# iterator=$(cut -d' ' -f$NSEQ <(echo {a..z}))
iterator=$(head -n$NSEQ iterators_list.txt | tail -n1)
printf "${iterator}, \".csv\")" >> population_combinatorics.jl
# echo "; println(fname_out)" >> population_combinatorics.jl
echo "" >> population_combinatorics.jl
### save
echo 'DelimitedFiles.writedlm(fname_out, hcat([' >> population_combinatorics.jl
for i in $(seq 1 $(echo $NSEQ - 1 | bc))
do
  # iterator=$(cut -d' ' -f$i <(echo {a..z}))
  iterator=$(head -n$i iterators_list.txt | tail -n1)
  echo "${iterator}," >> population_combinatorics.jl
done
# iterator=$(cut -d' ' -f$NSEQ <(echo {a..z}))
iterator=$(head -n$NSEQ iterators_list.txt | tail -n1)
echo "${iterator}" >> population_combinatorics.jl
echo "]...), ',')" >> population_combinatorics.jl
echo "end" >> population_combinatorics.jl
for i in $(seq 1 $NSEQ)
do
  echo "end" >> population_combinatorics.jl
done
echo "return(0)" >> population_combinatorics.jl
echo "end" >> population_combinatorics.jl

### Load GPAS_OUTPUT_STREAMLINED.csv and execute
echo 'DATA = CSV.read("GPAS_OUTPUT_STREAMLINED.csv"; missingstring="NA")' >> population_combinatorics.jl
echo '### listing all possible NSEQ conbinations...' >> population_combinatorics.jl
echo 'if sort(unique(DATA.TRAIN_POP)) != sort(unique(DATA.TEST_POP))' >> population_combinatorics.jl
echo '  println("The identities of the training populations are not the same as that of the validation populations!")' >> population_combinatorics.jl
echo '  println(string("Please check the input file: ", ARGS[1]))' >> population_combinatorics.jl
echo '  exit()' >> population_combinatorics.jl
echo 'end' >> population_combinatorics.jl
echo '### list all unique samples (popualtion name combinations)' >> population_combinatorics.jl
echo 'populations_list = sort(unique(DATA.TRAIN_POP))' >> population_combinatorics.jl
echo '### execute' >> population_combinatorics.jl
echo 'population_combinatorics(populations_list)' >> population_combinatorics.jl

### execute the julia script
julia population_combinatorics.jl