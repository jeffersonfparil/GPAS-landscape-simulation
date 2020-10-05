#############################################
###         GENOMIC PREDICTION (GP)       ###
### using Pool sequencing (Pool-seq) data ###
#############################################
### 2-step genomic prediction using "GWAlpha_SNPwise" output
### Step 1: Calculate the polygenic "risk" score using test population's genotype data (.bed,.bim,&.fam) and training populations's SNP effects (.gwas)
### Step 2: Predict the test population's phenotypes using the training population's linear model connecting the actual phenotypes with plink-derived polygenic "risk" scores, i.e.:
###             - y_train ~ po + p1*polygenic_score_train
###             - y_test_predicted = p0 + p1*polygenic_score_test
### Straight-forward predition for the non-itertive models (i.e. ["ML", "REML"]_["LS", "GLMNET"])
#############################################
### Inputs:
### (1) filename of the training population's output of genome-wide association study using pool sequencing data (Pool-GWAS)
###     - string(prefix, "_GWAlpha_SNPwise.gwas")
###     - string(prefix, "_GWAlpha_REML_LS_", fst_id, ".gwas")
###     - string(prefix, "_GWAlpha_REML_LS_", fst_id, ".ranef")
###     - string(prefix, "_GWAlpha_REML_RR_", fst_id, ".gwas")
###     - string(prefix, "_GWAlpha_REML_RR_", fst_id, ".ranef")
###     - string(prefix, "_GWAlpha_REML_GELMNET_", fst_id, ".gwas")
###     - string(prefix, "_GWAlpha_REML_GELMNET_", fst_id, ".ranef")
###     - string(prefix, "_GWAlpha_REML_LASSO_", fst_id, ".gwas")
###     - string(prefix, "_GWAlpha_REML_LASSO_", fst_id, ".ranef")
###     - *.gwas: HEADER: |col1:CHROM|col2:POS|col3:ALLELES|col4:FREQ|col5:ALPHA|col6:PVALUES|col7:LOD|
###     - *.ranef: HEADERLESS: |col1:random effects|
### (2) filename of the test population's full SNP data in sync format (FULLSNPSET_*.sync)
### Output:
### (1) Genomic prection cross-validation statistics (tab-delimited; |col01:TRAIN_name|col02:TRAIN_model|col03:TRAIN_rancovar|col04:TEST_name|col05:n_train|col06:n_test|col07:p0|col08:p1|col09:r2adj|co10l:y_true|col11:y_pred|)
###     - For GWAlpha_SNPwise: string(TRAIN_name, "_", TRAIN_model, "-", TEST_name, ".gp")
###     - For GWAlpha_REML_[LS, RR, GLMNET, LASSO]_[HIVERT, WEIRCOCK]: string(TRAIN_name, "_", TRAIN_model, "_", TRAIN_rancovar, "-", TEST_name, ".gp")

using DelimitedFiles
using Statistics
using LinearAlgebra
using GWAlpha

function poolgp_parallelization(TRAIN_fname_gwas, TEST_fname_sync)
    # # ARGS = ["POP_01_GWAlpha_SNPwise.gwas", "FULLSNPSET_POP_01.sync"]
    # # ARGS = ["POP_01_GWAlpha_MIXED_REML_WEIRCOCK.gwas", "FULLSNPSET_POP_02.sync"]
    # # ARGS = ["POP_01_GWAlpha_RIDGE.gwas", "FULLSNPSET_POP_02.sync"]
    # # ARGS = ["MULTIPOP_8_9961816979469537538_GWAlpha_SNPwise.gwas", "FULLSNPSET_MULTIPOP_8_9961816979469537538.sync"]
    # # ARGS = ["MULTIPOP_8_9961816979469537538_GWAlpha_REML_GLMNET_WEIRCOCK.gwas", "FULLSNPSET_POP_02.sync"]
    # TRAIN_fname_gwas = ARGS[1]
    # TEST_fname_sync = ARGS[2]

    ### start parsing training and test population names
    TRAIN_id = split(split(basename(TRAIN_fname_gwas), ".")[1], "_")
    TEST_name =  split(split(basename(TEST_fname_sync), "FULLSNPSET_")[2], ".")[1]

    ### parse test population sync into csv and load
    fname_geno_csv = string(join(split(TEST_fname_sync, ".")[1:end-1], "."), "_ALLELEFREQ.csv")
    if isfile(fname_geno_csv) == false
        GWAlpha.sync_processing_module.sync_parse(TEST_fname_sync)
    end
    GENO = DelimitedFiles.readdlm(fname_geno_csv, ',')

    ### load allelic effects from the training population
    gwas = DelimitedFiles.readdlm(TRAIN_fname_gwas, ',', header=true)[1]
    idx = [ sum(sum(reshape(GENO[i, 1:3],1,3) .== gwas[:,1:3], dims=2) .== 3) > 0 for i in 1:size(GENO,1) ]

    ### extract the testing genotype and the training allelic effects
    X = convert(Array{Float64,2}, GENO[idx,4:end])'
    b_hat = convert(Array{Float64,1}, gwas[:,5])

    ### prediction for iterative (GWAlpha_SNPwise) and non-iterative (GWAlpha_MIXED_REML_[WEIRCOCK, HIVERT], and GWAlpha_[RR, LASSO, GLMNET]) models
    if match(r"SNPwise", TRAIN_fname_gwas) != nothing
        println("Performing 2-step genomic prediction: y_pred = p0 + p1 * polygenic_score; where polygenic_score = X * b_hat_iterative")
        ### prepare names
        TRAIN_name = join(TRAIN_id[1:(end-2)], "_")
        TRAIN_model = join(TRAIN_id[(end-1):end], "_")
        TRAIN_rancovar = NaN
        TRAIN_fname_pheno = string(TRAIN_name, ".csv")
        TRAIN_fname_p0p1 = string(TRAIN_name, "_", TRAIN_model, ".p0p1")
        OUT_fname = string(TRAIN_name, "_", TRAIN_model, "-", TEST_name, ".gp")
        ### compute the polygenic scores and predict the phenotypes
        polygenic_score = X * b_hat
        ### solve y ~ p0 + p1 * polygenic_score
        ### NOTE: run auto-cross-validation first i.e. training == test or validation population
        if TRAIN_name == TEST_name
            y_train = DelimitedFiles.readdlm(TRAIN_fname_pheno, ',')[:,2]
            x_polygenic = hcat(repeat([1], inner=length(polygenic_score)), polygenic_score)
            p0, p1 = LinearAlgebra.inv(x_polygenic' * x_polygenic) * (x_polygenic' * y_train)
            y_train_centered = y_train .- mean(y_train); err = y_train - (x_polygenic * [p0,p1]); n = length(y_train)
            r2 = 1 - ((err' * err) / (y_train_centered' * y_train_centered))
            r2adj = 1 - ((1-r2) * ((n-1)/(n-1-1)))
            p0_p1_r2adj_n_train = (p0=p0,p1=p1,r2adj=r2adj, n_train=length(y_train))
            DelimitedFiles.writedlm(TRAIN_fname_p0p1, hcat(string.(keys(p0_p1_r2adj_n_train))...), '\t')
            io=open(TRAIN_fname_p0p1, "a"); DelimitedFiles.writedlm(io, hcat(p0_p1_r2adj_n_train...), '\t'); close(io)
            DelimitedFiles.readdlm(TRAIN_fname_p0p1, header=true)[1]
        end
        ### extract the solution to y ~ p0 + p1 * polygenic_score
        p0, p1, r2adj, n_train = DelimitedFiles.readdlm(TRAIN_fname_p0p1, header=true)[1]
        ### predict phenotypes
        y_pred = hcat(repeat([1], inner=length(polygenic_score)), polygenic_score) * [p0, p1]
        ### prepare output file
        y_test_true = DelimitedFiles.readdlm(string(TEST_name,".csv"), ',')[:,2]
        OUT = (TRAIN_name=repeat([TRAIN_name],length(y_pred)),
                TRAIN_model=repeat([TRAIN_model],length(y_pred)),
                TRAIN_rancovar=repeat([TRAIN_rancovar],length(y_pred)),
                TEST_name=repeat([TEST_name],length(y_pred)),
                n_train=repeat([n_train],length(y_pred)),
                n_test=repeat([length(y_pred)],length(y_pred)),
                p0=repeat([p0],length(y_pred)),
                p1=repeat([p1],length(y_pred)),
                r2adj=repeat([r2adj],length(y_pred)),
                y_true=y_test_true,
                y_pred=y_pred)
        DelimitedFiles.writedlm(OUT_fname, hcat(string.(keys(OUT))...))
        io=open(OUT_fname, "a"); DelimitedFiles.writedlm(io, hcat(OUT...)); close(io)
    else
        println("Performing genomic prediction: y_pred = X*b_hat")
        if match(r"MIXED", TRAIN_fname_gwas) != nothing
            ### prepare names
            TRAIN_name = join(TRAIN_id[1:(end-4)], "_")
            TRAIN_model = join(TRAIN_id[(end-3):(end-1)], "_")
            TRAIN_rancovar = TRAIN_id[end]
            OUT_fname = string(TRAIN_name, "_", TRAIN_model, "_", TRAIN_rancovar, "-", TEST_name, ".gp")
        else
            ### prepare names
            TRAIN_name = join(TRAIN_id[1:(end-2)], "_")
            TRAIN_model = join(TRAIN_id[(end-1):end], "_")
            TRAIN_rancovar = NaN
            OUT_fname = string(TRAIN_name, "_", TRAIN_model, "-", TEST_name, ".gp")
        end
        ### predict phenotype
        y_pred = (hcat(repeat([1.00], size(X,1)), X) * b_hat)
        ### prepare output file
        y_test_true = DelimitedFiles.readdlm(string(TEST_name,".csv"), ',')[:,2]
        OUT = (TRAIN_name=repeat([TRAIN_name],length(y_pred)),
                TRAIN_model=repeat([TRAIN_model],length(y_pred)),
                TRAIN_rancovar=repeat([TRAIN_rancovar],length(y_pred)),
                TEST_name=repeat([TEST_name],length(y_pred)),
                n_train=repeat([size(DelimitedFiles.readdlm(string(TRAIN_name, ".csv")),1)],length(y_pred)),
                n_test=repeat([length(y_pred)],length(y_pred)),
                p0=repeat([NaN],length(y_pred)),
                p1=repeat([NaN],length(y_pred)),
                r2adj=repeat([NaN],length(y_pred)),
                y_true=y_test_true,
                y_pred=y_pred)
        DelimitedFiles.writedlm(OUT_fname, hcat(string.(keys(OUT))...))
        io=open(OUT_fname, "a"); DelimitedFiles.writedlm(io, hcat(OUT...)); close(io)
    end
    return 0
end