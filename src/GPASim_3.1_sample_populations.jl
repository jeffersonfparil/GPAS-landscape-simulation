##########################################################################
### SAMPLE MULTIPLE POPULATIONS CONSTRAINED BY THE NUMBER OF LIBRARIES ###
##########################################################################
### NOTE: divide the landscape into equally-sized rectangles and sample the central population in a jittered way
###       then list the sample sizes for each population constrained by the number of libraries (nLib)
### load libraries
using DelimitedFiles
### inputs: NOTE: length(pop_names) should have a whole number square-root! (quantinemo2 in GPASim_1.0_simulate.sh would have failed otherwise! See line GPASim_1.0_simulate.sh:106)
# ARGS = ["POP_PREFIXES.temp", "1000"]
pop_names = convert(Array{String,1}, DelimitedFiles.readdlm(ARGS[1])[:,1])
nLib = parse(Int64, ARGS[2])
n = Int(sqrt(length(pop_names)))
### prepare square matrix landscape (L where populations are arranged by row)
L = permutedims(reshape(pop_names, n, n))
### dividing the square
nrows = []
ncols = []
for i in 1:n
  for j in 1:n
    append!(nrows, Int(floor(n/i)))
    append!(ncols, Int(floor(n/j)))
  end
end
# F for factorization matrix: |col1:nrows|col2:ncols|col3:npops|col4:abs(nrows-ncols)|
F = hcat(nrows, ncols, nrows .* ncols, abs.(nrows .- ncols))
F = F[F[:,3] .>= 4, :]                ### at least 4 population samples (less than 4 and we won't have enough degress of freedom to build models)
F = F[F[:,1] .>= 2, :]                ### at least two divisions across rows
F = F[(mod.(n, F[:,1]) .== 0) .& (mod.(n, F[:,2]) .== 0), :] ### exclude imperfect divisions across rows and columns
F = F[reverse(sortperm(F[:,1])),:]    ### prioritize more divisions across rows since we simulated the QTL diffusion gradient to be across rows
F = F[sortperm(F[:,4]),:]             ### sort by increaing abs(nrow-ncol)
F = F[sortperm(F[:,3]),:]             ### sort by increasing number of population samples
# D for division matrix: |col1:nrows|col2:ncols|col3:npops|
D = zeros(Int64, length(unique(F[:,3])), 3)
for i in 1:length(unique(F[:,3]))
  # for each number of population samples select the factorization which minimises the difference between nrows and ncols;
  # which also maximizes nrows where we expect the QTL diffusion gradient to be parallel with
  D[i,:] = F[F[:,3] .== unique(F[:,3])[i],:][1, 1:3]
end
### sampling the central(-ish) population per bounding rectangle
POP_NAMES = []
for idx in 1:size(D,1)
  # idx = 4
  nRows = D[idx,1]
  nCols = D[idx,2]
  nPops = D[idx,3]
  # jitter to center(-ish) or to the boundaries (bottom for rows and right for columns)
  for jitterer_1 in [0, 1]
    for jitterer_2 in [0, 1]
      idx_rows = repeat(cumsum(repeat([Int(round(n/nRows))], nRows)) .- (jitterer_1*Int(round(n/nRows/2))), inner=nCols) ### upper limit of each row minus number of steps to central-ish + repeat to fit with the always as large or smaller number of columns
      idx_cols = repeat(cumsum(repeat([Int(round(n/nCols))], nCols)) .- (jitterer_2*Int(round(n/nCols/2))), outer=nRows) ### upper limit of each column minus number of steps to central-ish + repeat to fit with the always as large or larger number of rows
      names_array = []
      for i in 1:(nRows*nCols)
        push!(names_array, L[idx_rows[i], idx_cols[i]])
      end
      push!(POP_NAMES, join(names_array, ";"))
    end
  end
  # jitter to center(-ish) or to the boundaries (top for rows and left for columns)
  for jitterer_1 in [0, 1]
    for jitterer_2 in [0, 1]
      idx_rows = repeat(cumsum(repeat([Int(round(n/nRows))], nRows)) .- (Int(round(n/nRows))-1) .+ (jitterer_1*Int(round(n/nRows/2))), inner=nCols) ### upper limit of each row minus number of steps to central-ish + repeat to fit with the always as large or smaller number of columns
      idx_cols = repeat(cumsum(repeat([Int(round(n/nCols))], nCols)) .- (Int(round(n/nCols))-1) .+ (jitterer_2*Int(round(n/nCols/2))), outer=nRows) ### upper limit of each column minus number of steps to central-ish + repeat to fit with the always as large or larger number of rows
      names_array = []
      for i in 1:(nRows*nCols)
        push!(names_array, L[idx_rows[i], idx_cols[i]])
      end
      push!(POP_NAMES, join(names_array, ";"))
    end
  end
end
unique!(POP_NAMES)
### output population names per across population samples and the corresponding number of individuals per population
for sample in POP_NAMES
  # sample = POP_NAMES[end]
  populations = split(sample, ";")
  population_count = length(populations)
  population_sizes = repeat([Int(floor(nLib / population_count))], population_count)
  # add remainder samples (if mod(nLib, population_sizes[1])!=0) to the central(-ish) population
  population_sizes[Int(round(population_count/2))] = population_sizes[Int(round(population_count/2))] + (nLib - sum(population_sizes))
  # ouput: population ID and sample size per population: headerless, tab-delimited, |col1:pop_ID|col2:sample_sizes|
  fname_out = string("MULTIPOP_", population_count, "_", hash(sample), ".idx")
  isfile(fname_out) ? rm(fname_out) : nothing
  DelimitedFiles.writedlm(fname_out, hcat(populations, population_sizes))
end
