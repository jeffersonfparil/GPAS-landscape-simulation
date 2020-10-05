### Build the ABC optimization input
###   specifically the summary statistics of the GPAS performance metrics
###   across different combinations of the 4 sampling strategies:
###       (1) Across population sampling with individual genotyping (Indi-seq)
###       (2) Across population sampling with pool genotyping (Pool-seq)
###       (3) Within population sampling with Indi-seq
###       (4) Within population sampling with Pool-seq
### NOTE: Each instance of a sampling strategy combination has a fixed number of sampling-x-genotyping samples,
###       i.e. the number of sequencing experiments (the limiting factor or resource; NSEQ=10 is a reasonable number, where each of the 10 multiplexed libraries has a maximum of 1000 individually barcoded libraries)

### Inputs:
### (1) GPAS_OUTPUT_STREAMLINED.csv (comma-separated; header: |col01:SAMPLING_SCHEME|col02:GENOTYPING_SCHEME|col03:TRAIN_POP|col04:TEST_POP|col05:NTRAIN|\
###                                                           |col06:NTEST|col07:MODEL|col08:COVARIATE|col09:AUC|col10:AUC_CORRECTED|\
###                                                           |col11:BONFERRONI_5PERCENT_TRUE_POSITIVE|col12:BONFERRONI_5PERCENT_FALSE_POSITIVE|col13:BONFERRONI_5PERCENT_QTL_ID|col14:RMSE|col15:CORRELATION|)
### (3) maximum number of sequencing libraries
### (2) number of simulated QTL
### Outputs:
### (1) ABC_OPTIM_INPUT_SUMMSTATS.csv (sampling strategies across random combinatorial sensible samples and the corresponsing GPAS performance summaries;
###                                  comma-separated; header: |col01:ACROSS_INDI|col02:ACROSS_POOL|col03:WITHIN_INDI|col04:WITHIN_POOL|col05:AUC_MEAN|col06:AUC_CORRECTED_MEAN|col07:RMSE_MEAN|col08:CORRELATION_MEAN|col09:TRUE_POSITIVE_RATE|col10:FALSE_DISCOVERY_RATE|)
### (2) POPULATION_COMBINATIONS_*.csv (output of the adaptive function in "population_combinatorics.jl"; comma-separated; headerless)
### (3) SUMMSTATS_*.csv (GPAS performance summary across random combinatorial sensible samples; comma-separated; headerless)

### load libraries
using DelimitedFiles
using DataFrames
using CSV
using Random
using Statistics
function random_combinatorics_and_summary_statistics(fname_input::String, MAX_NSEQ::Int64, NQTL::Int64, SRC_DIR::String, exhaustive::Bool=false)
    # fname_input="GPAS_OUTPUT_STREAMLINED.csv"; NQTL=10; SRC_DIR="/data/Lolium/Softwares/genomic_prediction/src"; exhaustive=false
    ### load data
    DATA = CSV.read(fname_input; missingstring="NA")
    DATA = DATA[.!ismissing.(DATA.TRAIN_POP), :]
    ### listing all possible NSEQ conbinations...
    if sort(unique(DATA.TRAIN_POP)) != sort(unique(DATA.TEST_POP))
      println("The identities of the training populations are not the same as that of the validation populations!")
      println(string("Please check the input file: ", ARGS[1]))
      exit()
    end
    ### list all unique samples (popualtion name combinations)
    populations_list = sort(unique(DATA.TRAIN_POP))
    println("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
    println("LIST ALL SENSIBLE POPULATION SAMPLING COMBINATIONS (no duplicates within each sampling strategy)")
    println("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
    ### iterate across different number of sequencing experiments (maximum is the number of available samples)
    for NSEQ in 1:minimum([MAX_NSEQ, length(populations_list)])
      println(string("Number of sequencing experiments = ", NSEQ))
      ### NOTE: The populations in each sample are unique (no duplicates for each sample!)
      if exhaustive == true
        println("##################################")
        println("Performing exhausive combinatorics")
        println("##################################")
        ### combinatorics in parallel output csv files with thw prefix "POPULATION_COMBINATIONS_"
        ### run in shell: GPASim_6.2_adaptive_combinatorics.sh which generates: POPULATION_COMBINATIONS_*
        run(`$SRC_DIR/GPASim_6.2_adaptive_combinatorics.sh $NSEQ`)
        ### load and merge the output population combination indices
        ls_dir = readdir()
        idx_fnames = ls_dir[.!isnothing.(match.(r"POPULATION_COMBINATIONS_", ls_dir))]
        idx = []
        for f in idx_fnames
          append!(idx, convert(Array{Int64},DelimitedFiles.readdlm(f, ',')))
        end
        rm.(idx_fnames) ### clean-up
        ### list all population combinations
        IDX = reshape(idx, NSEQ, Int(length(idx)/NSEQ))'
        POPULATION_COMBINATIONS = populations_list[IDX]
      else
        println("###################################")
        println("Performing randomised combinatorics")
        println("###################################")
        ### random sampling in julia
        # #@@@ tests @@@@@
        # using CSV
        # populations_list = sort(unique(CSV.read("GPAS_OUTPUT_STREAMLINED.csv"; missingstring="NA").TRAIN_POP))
        # NSEQ=5
        # #@@@@@@@@@@@@@@@
        ### separate within and across population samples
        names_within = populations_list[length.(split.(populations_list, ":")) .== 1]
        names_across = populations_list[length.(split.(populations_list, ":")) .> 1]
        ### set the target number of elements from which we will be random sampling from
        n = maximum( [2*length(names_across), 2*NSEQ] )
        ### calcuate the maximum number of combinations we can achieve to limit computation
        ### as a function of the within poulation samples (the ones we are certain will all be disjoint samples)
        n_across = length(names_across)
        n_within = n - n_across
        iter_max = factorial(BigInt(n_within)) / ( factorial(BigInt(NSEQ)) * factorial(BigInt(n_within-NSEQ)) )
        ### set the target number of sensible samples (capped to at most 10,000 samples!)
        n_output = minimum([iter_max, 10000])
        ### random sensible sampling
        OUTPUT = []
        while length(OUTPUT) < n_output*NSEQ
          ### fill-up the set for random sampling with all the across population samples first (because population samples are always less that the within population samples)
          names_merged = copy(names_across)
          ### then supplement with the within population sample to reach the target set size, n, for each 
          append!(names_merged, names_within[Random.randperm!(collect(1:length(names_within)))][1:n_within])

          rand_sample = names_merged[Random.randperm!(collect(1:n))][1:NSEQ]
          if length(unique(vcat(split.(rand_sample, ":")...))) == length(vcat(split.(rand_sample, ":")...))
            append!(OUTPUT, rand_sample)
          end
        end
        POPULATION_COMBINATIONS = permutedims(reshape(OUTPUT, NSEQ, Int(n_output)))
      end
      println("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
      println("EXTRACT THE GPAS PERFORMANCES WITHIN EACH SAMPLING STRATEGY")
      println("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
      @time @sync @distributed for i in 1:size(POPULATION_COMBINATIONS,1)
        # i=100
        # println(i)
        ### subset the populations
        pop_list = POPULATION_COMBINATIONS[i,:]
        subdata = DATA[[sum(pop_list .== x)>0 for x in DATA.TRAIN_POP], :]
        ### missing data will pop-up here in there due to convergence errors in GCTA - we are just going to exclude them
        SUBDATA = subdata[.!ismissing.(subdata.SAMPLING_SCHEME), :]
        if size(SUBDATA,1) > 0
          ### identify the sampling strategies available in this population combination random sample
          sampling_strategies = string.(SUBDATA.SAMPLING_SCHEME, "-", SUBDATA.GENOTYPING_SCHEME)
          N_across_indi = sum(sampling_strategies .== "ACROSS-INDIVIDUAL")
          N_across_pool = sum(sampling_strategies .== "ACROSS-POOL")
          N_within_indi = sum(sampling_strategies .== "WITHIN-INDIVIDUAL")
          N_within_pool = sum(sampling_strategies .== "WITHIN-POOL")
          ### random sampling of populations allocated to the available sampling strategies
          ###### random numbers multiplied by the number of observed sampling srategies in SUBDATA
          RAND_across_indi = Random.rand(1)[1] * N_across_indi
          RAND_across_pool = Random.rand(1)[1] * N_across_pool
          RAND_within_indi = Random.rand(1)[1] * N_within_indi
          RAND_within_pool = Random.rand(1)[1] * N_within_pool
          ###### random proportions based on the above
          RAND_P_across_indi = RAND_across_indi / sum(RAND_across_indi + RAND_across_pool + RAND_within_indi + RAND_within_pool)
          RAND_P_across_pool = RAND_across_pool / sum(RAND_across_indi + RAND_across_pool + RAND_within_indi + RAND_within_pool)
          RAND_P_within_indi = RAND_within_indi / sum(RAND_across_indi + RAND_across_pool + RAND_within_indi + RAND_within_pool)
          RAND_P_within_pool = RAND_within_pool / sum(RAND_across_indi + RAND_across_pool + RAND_within_indi + RAND_within_pool)
          ###### number of populations allocated to each of the available sampling strategies
          RAND_NPOP_across_indi = Int64(round(RAND_P_across_indi * NSEQ))
          RAND_NPOP_across_pool = Int64(round(RAND_P_across_pool * NSEQ))
          RAND_NPOP_within_indi = Int64(round(RAND_P_within_indi * NSEQ))
          RAND_NPOP_within_pool = maximum([0, NSEQ - (RAND_NPOP_across_indi + RAND_NPOP_across_pool + RAND_NPOP_within_indi)]) ### making restrcting to positive whole numbers
          ###### fraction of populations allocated to the available sampling strategies
          P_across_indi = RAND_NPOP_across_indi / NSEQ
          P_across_pool = RAND_NPOP_across_pool / NSEQ
          P_within_indi = RAND_NPOP_within_indi / NSEQ
          P_within_pool = RAND_NPOP_within_pool / NSEQ
          ###### identify these randomly sample populations and their corresponding sample strategy allocation
          rand_pop_list = pop_list[Random.randperm(length(pop_list))]
          rand_sampling_scheme = vcat(
            repeat(["ACROSS"], RAND_NPOP_across_indi),
            repeat(["ACROSS"], RAND_NPOP_across_pool),
            repeat(["WITHIN"], RAND_NPOP_within_indi),
            repeat(["WITHIN"], RAND_NPOP_within_pool)
            )
          rand_genotyping_scheme = vcat(
            repeat(["INDIVIDUAL"], RAND_NPOP_across_indi),
            repeat(["POOL"], RAND_NPOP_across_pool),
            repeat(["INDIVIDUAL"], RAND_NPOP_within_indi),
            repeat(["POOL"], RAND_NPOP_within_pool)
            )
          ### merge these random samples for summarizing GPAS performance
          summarize_me = []
          for j in 1:length(rand_pop_list)
            # j = 1
            push!(summarize_me, SUBDATA[(SUBDATA.TRAIN_POP .== rand_pop_list[j]) .&
                                        (SUBDATA.SAMPLING_SCHEME .== rand_sampling_scheme[j]) .&
                                        (SUBDATA.GENOTYPING_SCHEME .== rand_genotyping_scheme[j]), :]
                  )
          end
          SUMMARIZE_ME = vcat(summarize_me...)
          ### extract summary statistics
          AUC_MEAN = Statistics.mean(SUMMARIZE_ME.AUC[.!ismissing.(SUMMARIZE_ME.AUC)])
          AUC_CORRECTED_MEAN = Statistics.mean(SUMMARIZE_ME.AUC_CORRECTED[.!ismissing.(SUMMARIZE_ME.AUC_CORRECTED)])
          RMSE_MEAN = Statistics.mean(SUMMARIZE_ME.RMSE[.!ismissing.(SUMMARIZE_ME.RMSE)])
          CORRELATION_MEAN = Statistics.mean(SUMMARIZE_ME.CORRELATION[.!ismissing.(SUMMARIZE_ME.CORRELATION)])
          bonferroni_5percent_true_positive = unique(vcat(split.(SUMMARIZE_ME.BONFERRONI_5PERCENT_TRUE_POSITIVE[.!ismissing.(SUMMARIZE_ME.BONFERRONI_5PERCENT_TRUE_POSITIVE)],":")...))
          bonferroni_5percent_true_positive = bonferroni_5percent_true_positive[.!ismissing.(bonferroni_5percent_true_positive)]
          bonferroni_5percent_false_positive = unique(vcat(split.(SUMMARIZE_ME.BONFERRONI_5PERCENT_FALSE_POSITIVE[.!ismissing.(SUMMARIZE_ME.BONFERRONI_5PERCENT_FALSE_POSITIVE)],":")...))
          bonferroni_5percent_false_positive = bonferroni_5percent_false_positive[.!ismissing.(bonferroni_5percent_false_positive)]
          bonferroni_5percent_qtl = unique(vcat(split.(SUMMARIZE_ME.BONFERRONI_5PERCENT_QTL_ID[.!ismissing.(SUMMARIZE_ME.BONFERRONI_5PERCENT_QTL_ID)],":")...))
          bonferroni_5percent_qtl = bonferroni_5percent_qtl[.!ismissing.(bonferroni_5percent_qtl)]
          TRUE_POSITIVE_RATE = length(bonferroni_5percent_qtl) / NQTL
          FALSE_DISCOVERY_RATE = length(bonferroni_5percent_false_positive) / (length(bonferroni_5percent_false_positive) + length(bonferroni_5percent_true_positive))
          out = (NSEQ=NSEQ,
                ACROSS_INDI=P_across_indi,
                ACROSS_POOL=P_across_pool,
                WITHIN_INDI=P_within_indi,
                WITHIN_POOL=P_within_pool,
                AUC_MEAN=AUC_MEAN,
                AUC_CORRECTED_MEAN=AUC_CORRECTED_MEAN,
                RMSE_MEAN=RMSE_MEAN,
                CORRELATION_MEAN=CORRELATION_MEAN,
                TRUE_POSITIVE_RATE=TRUE_POSITIVE_RATE,
                FALSE_DISCOVERY_RATE=FALSE_DISCOVERY_RATE
              )
          fname_out = string("SUMMSTATS_", string(hash(join(pop_list, "-"))), "-", Int(round(rand()*1e10)), ".csv")
          DelimitedFiles.writedlm(fname_out, hcat(out...), ',')
        end ### non-empty SUBDATA
      end ### parallel for-loop
    end ### iterate 1:NSEQ sequencing experiments
    ### load and merge the ouput summary statistics for each sampling strategy
    ls_dir = readdir()
    OUT = []
    summstats_fnames = ls_dir[.!isnothing.(match.(r"SUMMSTATS_", ls_dir))]
    for f in summstats_fnames
      push!(OUT, DelimitedFiles.readdlm(f, ','))
    end
    ### write-out
    FNAME_OUTPUT = "ABC_OPTIM_INPUT_SUMMSTATS.csv"
    DelimitedFiles.writedlm(FNAME_OUTPUT, hcat(["NSEQ",
                                                "ACROSS_INDI",
                                                "ACROSS_POOL",
                                                "WITHIN_INDI",
                                                "WITHIN_POOL",
                                                "AUC_MEAN",
                                                "AUC_CORRECTED_MEAN",
                                                "RMSE_MEAN",
                                                "CORRELATION_MEAN",
                                                "TRUE_POSITIVE_RATE",
                                                "FALSE_DISCOVERY_RATE"]...), ',')
    io=open(FNAME_OUTPUT,"a"); DelimitedFiles.writedlm(io, vcat(OUT...), ','); close(io)
    return(0)
end
