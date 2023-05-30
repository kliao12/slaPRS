# slaPRS
Here is the repository for our method for constructing polygenic risk scores in admixed individuals: stacking local ancestry PRS (slaPRS). See https://zenodo.org/record/7987848#.ZHZS0C1h30o for link to download example plink files (too large to host on github)

## Installation Instructions
devtools::install_github('kliao12/slaPRS') <br>
library(bigsnpr) #method uses bigsnpr package to handle genotype/phenotype data and perform C+T in R <br>
library(data.table) <br>
library(caret); library(glmnet) #method uses native R packages for elastic net fitting <br>

## Example workflow
library(slaPRS)

### 1) Import example plink files containing genotype and phenotype (in fam file) information using bigsnpr package
bigSNP_file <- load_plinkFile("./path/to/example_plink_files.bed") #Use this is importing plink file 

str(bigSNP_file) #Look at inputted format

#Get alias for commonly used objects
G   <- bigSNP_file$genotypes; CHR <- bigSNP_file$map$chromosome
orig_POS <- bigSNP_file$map$physical.pos; 

### 2) Import summary statistics for each population with columns: chr, pos, a0, a1, beta, p val
### Note: example_sumstats provided has both populations in one file so need to separate
all_sumstats <- slaPRS::example_sumstats

pop1_stats <- all_sumstats[,c(2,3,4,5,6,7)]
pop2_stats <- all_sumstats[,c(2,3,4,5,8,9)]

### 3) Keep only sites in common between genotype and summary statistics
N <- length(bigSNP_file$fam$sample.ID)
used.bigSNP <- snp_attach(snp_subset(bigSNP_file, ind.row = rows_along(bigSNP_file$fam), ind.col = which(orig_POS %in% all_sumstats$POS)))
G <- used.bigSNP$genotypes; CHR <- used.bigSNP$map$chromosome; POS <- used.bigSNP$map$physical.pos

### 2) Match alleles between summary statistics and plink/rds file. Note not actually needed b/c simulations ensured matched alleles but used to illustrate
matched_betas_pvals <- match_stats(used.bigSNP, pop1_stats, pop2_stats)

### 3) Get local ancestry in each window across chromosomes. Note this can take a few minutes
all_la <- slaPRS::example_la

start <- POS[1]; end <- POS[length(POS)]
local_la_df <- get_la_example(all_la, 5000000) #Right now have to run this by hand, github won't load. Maybe need to run document() first

### 4) Split into training and testing samples 
ind.train <- sample(nrow(G), 7000); ind.test <- setdiff(rows_along(G), ind.train)

### 5) Compute stacked PRS in testing samples
stacked_PRS <- run_slaPRS(used.bigSNP, matched_betas_pvals, local_la_df, full_model = TRUE, ind.train, ind.test)

### 6) Compare to C+T in testing samples 
CT_PRS_list <- run_CT(used.bigSNP, matched_betas_pvals, rsquared = 0.10, ind.train, ind.test)

### 7) Compute correlation for both
test_y <- as.numeric(bigSNP.obj$fam$affection)[ind.test]
stacked_cor <- cor.test(stacked_PRS, test_y)$estimate
pop1_cor <- cor.test(CT_PRS_list$pop1_CT, test_y)$estimate; pop2_cor <- cor.test(CT_PRS_list$pop2_CT, test_y)$estimate

print(paste0("cor(stacked PRS, y) = ", stacked_cor))
print(paste0("cor(pop1 C+T PRS, y) = ", pop1_cor)); print(paste0("cor(pop2 C+T PRS, y) = ", pop2_cor))
