##################################
### Landscape Characterization ###
##################################
### Inputs:
# ARGS = ["LANDSCAPE.sync", "100"]
filename_sync = ARGS[1]            ### filename of the Pool-seq data where each population is a single pool
n_pop = parse(Int64, ARGS[2])   ### total number of populations in the landscape
### Ouputs:
### (1) Phenotypic summary statistics per population (comma-separated; header: |col1:POP|col2:SIZE|col3:MEAN|col4:VAR|; string(join(split(ARGS[1], ".")[1:end-1], "."), ".stat"))
### (2) Pairwise fixation indices (comma-separated; headerless; string(join(split(ARGS[1], ".")[1:end-1], "."), ".fst"))
### (3) Heatmaps of phenotype means, variances and Fst ("LANDSCAPE.png"; portable network graphics format)

using DelimitedFiles
using Statistics
using RCall
using Distributed
Distributed.addprocs(length(Sys.cpu_info()))
@everywhere using GWAlpha

### Phenotype variation across the landscape
POP=[]; SIZE=[]; MEAN=[]; VAR=[];
for f in readdir()[.!isnothing.(match.(r"POP_", readdir())) .& .!isnothing.(match.(r".fam", readdir()))]
  # f = readdir()[.!isnothing.(match.(r"POP_", readdir())) .& .!isnothing.(match.(r".fam", readdir()))][1]
  y = DelimitedFiles.readdlm(f)[:,end]
  push!(POP, join(split(f, ".")[1:end-1], "."))
  push!(SIZE, length(y))
  push!(MEAN, Statistics.mean(y))
  push!(VAR, Statistics.var(y))
end
OUT = (POP=POP,SIZE=SIZE,MEAN=MEAN,VAR=VAR)
out_fname = string(join(split(ARGS[1], ".")[1:end-1], "."), ".stat")
DelimitedFiles.writedlm(out_fname, hcat(string.(keys(OUT))...), ',')
io = open(out_fname, "a"); DelimitedFiles.writedlm(io, hcat(OUT...), ','); close(io)
### Fst across the landscape
GWAlpha.relatedness_module.Fst_pairwise(filename_sync=filename_sync, window_size=100000, pool_sizes=repeat([100], n_pop))
mv(string(join(split(ARGS[1], ".")[1:end-1], "."), "_COVARIATE_FST.csv"), string(join(split(ARGS[1], ".")[1:end-1], "."), ".fst"), force=true)
FST = DelimitedFiles.readdlm(string(join(split(ARGS[1], ".")[1:end-1], "."), ".fst"), ',')
### Plotting
@rput n_pop
@rput POP
@rput MEAN
@rput VAR
@rput FST
R"POP=unlist(POP); MEAN=unlist(MEAN); VAR=unlist(VAR); FST=unlist(FST)"
R"VAR_min = min(VAR)"
R"VAR_max = max(VAR)"
R"FST_min = min(FST)"
R"FST_max = max(FST)"
R"VAR = (VAR-VAR_min) / (VAR_max-VAR_min)"
R"FST = (FST-FST_min) / (FST_max-FST_min)"
R"png('LANDSCAPE.png', width=4000, height=2000)"
R"par(mfrow=c(2,2))"
R"m = sqrt(n_pop)"
R"x=rep(1:m, times=m)"
R"y=rep(1:m, each=m)"
R"n_colors = 1e3"
R"par(cex=2)"
R"colors = colorRampPalette(c('#f7fcf0', '#e0f3db', '#ccebc5', '#a8ddb5', '#7bccc4', '#4eb3d3', '#2b8cbe', '#0868ac', '#084081'))(n_colors)"
### population names
R"plot(x=c(0,m+1), y=c(0,m+1), type='n', xaxt='n', yaxt='n', xlab=NA, ylab=NA, main='Landscape of Populations'); text(x, y, lab=POP)"
### means
R"image(matrix(MEAN, nrow=m, byrow=FALSE), col=colors[ceiling(c(min(MEAN*n_colors)):ceiling(max(MEAN*n_colors)))], xaxt='n', yaxt='n', main=paste0('Phenotype Means\n', '(min=', round(min(MEAN),4), '; max=', round(max(MEAN),4), ')'))"
R"par(new=TRUE);plot(x=c(1,m), y=c(1,m), type='n', xaxt='n', yaxt='n', xlab=NA, ylab=NA); text(x, y, lab=round(MEAN,2))"
### scaled variances
R"image(matrix(VAR, nrow=m, byrow=FALSE), col=colors[ceiling(c(min(VAR*n_colors)):ceiling(max(VAR*n_colors)))], xaxt='n', yaxt='n', main=paste0('Phenotype Variances (Mapped into 0.0-1.0 scale)\n', '(min=', round(VAR_min,4), '; max=', round(VAR_max,4), ')'))"
R"par(new=TRUE);plot(x=c(1,m), y=c(1,m), type='n', xaxt='n', yaxt='n', xlab=NA, ylab=NA); text(x, y, lab=round(VAR,2))"
### Fst
R"par(cex=1)"
R"image(FST, col=colors[ceiling(c(min(FST*n_colors)):ceiling(max(FST*n_colors)))], xaxt='n', yaxt='n', main='Pairwise Fixation Indices (Mapped into 0.0-1.0 scale)')"
R"dev.off()"
