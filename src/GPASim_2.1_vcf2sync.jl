#####################################
### VCF TO SYNC CONVERSION SCRIPT ###
#####################################
### load libraries
using GeneticVariation
using DelimitedFiles
using Distributions
### arguments
# ARGS = ["POP_01.vcf", "POP_01.fam", "5"]
# ARGS = ["POP_01.vcf", "POP_01.fam", "1"]
fname = ARGS[1] ### vcf file name
fam = DelimitedFiles.readdlm(ARGS[2])[:,end] ### fam file
npools = parse(Int64, ARGS[3])
out_fname = string(split(fname, ".vcf")[1], ".sync") ### output sync filename
out_phen_py_fname = string(split(fname, ".vcf")[1], ".py") ### output phenotype py-format filename
out_phen_csv_fname = string(split(fname, ".vcf")[1], ".csv") ### output phenotype csv-format filename (COL1: pool_sizes, COL2: means)
### sort individuals and dividie into pools
pheno = (hcat(collect(1:length(fam)), fam)[sortperm(fam), :])
n = size(pheno,1)
pool_sizes = repeat([Int(floor(n/npools))], npools)
pool_sizes[Int(ceil(npools/2))] = pool_sizes[Int(ceil(npools/2))] + (n - sum(pool_sizes)) ### add remainder if any to the central(ish) pool
cum_pool_sizes = cumsum(pool_sizes)
npools > 1 ? (perc = (cum_pool_sizes ./ sum(pool_sizes))[1:end-1]) : (perc = (cum_pool_sizes ./ sum(pool_sizes)))
q = Distributions.percentile(convert(Array{Float64,1},pheno[:,2]), perc)
pool_means = []
idx = []
for p in 1:npools
    ini = append!([1], cum_pool_sizes)[p]
    fin = append!([1], cum_pool_sizes)[p+1]
    push!(idx, convert(Array{Int64,1},pheno[ini:fin, 1]))
    append!(pool_means, mean(pheno[ini:fin, 2]))
end
### define parsing function
function parse_vcf_2_sync_per_locus(record::GeneticVariation.VCF.Record, individual_ids::Array{Int64,1})::String
    sort!(individual_ids)
    gen = vcat(GeneticVariation.VCF.genotype(record)...)[individual_ids]
    ref = parse(Int, GeneticVariation.VCF.ref(record))
    alt = 3 - ref ### ref and alt are either 1 or 2
    n_alt = sum(parse.(Int, vcat(split.(gen, "/")...)))
    n_ref = 2*length(gen) - n_alt
    sync_col4 = zeros(Int, 6)
    sync_col4[[ref,alt]] = [n_ref, n_alt]
    sync_col4 = join(string.(sync_col4), ":")
    return(sync_col4)
end
### open vcf stream
vcf = GeneticVariation.VCF.Reader(open(fname, "r"));
# RECORD = [VCF.Record("20\t14370\trs6054257\tG\tA\t29\tPASS\tNS=3;DP=14;AF=0.5;DB;H2\tGT:GQ:DP:HQ\t0|0:48:1:51,51\t1|0:48:8:51,51")]
### open sync stream
if isfile(out_fname)
    rm(out_fname) ### remove if it exists
end
io = open(out_fname, "a")
### iterate across loci
for record in vcf
  # RECORD[1] = record
  # record = RECORD[1]
  chr = GeneticVariation.VCF.chrom(record)
  pos = GeneticVariation.VCF.pos(record)
  ### ref and alt alleles are simulated as "1" and "2" and are not canonical vcf encoding
  ref = parse(Int, GeneticVariation.VCF.ref(record))
  alt = 3 - ref ### ref and alt are either 1 or 2
  ### parse allele counts per pool
  sync_row_out = [chr, pos, ref]
  for pool in 1:npools
      individual_ids = idx[pool]
      push!(sync_row_out, parse_vcf_2_sync_per_locus(record, individual_ids))
  end
  ### modify ref allele to either "A" for 1 and "T" for 2
  sync_row_out[3] = ["A", "T"][sync_row_out[3]]
  ### append into the sync stream
  DelimitedFiles.writedlm(io, hcat(sync_row_out...))
end
close(io)
close(vcf)
### output pooled phenotype files in py and csv formats
line1 = "Pheno_name='pheno';"
line2 = string("sig=", std(pheno[:,2]), ";")
line3 = string("MIN=", minimum(pheno[:,2]), ";")
line4 = string("MAX=", maximum(pheno[:,2]), ";")
line5 = string("perc=[", join(perc, ","), "];")
line6 = string("q=[", join(q, ","), "];")
PY_PHEN = vcat(line1, line2, line3, line4, line5, line6)
CSV_PHEN = hcat(pool_sizes, pool_means)
DelimitedFiles.writedlm(out_phen_py_fname, PY_PHEN)
DelimitedFiles.writedlm(out_phen_csv_fname, CSV_PHEN, ',')
