# ---------------------------------------------------------------
# Test CellDMC and TCA false positive rate
# ---------------------------------------------------------------

## Aim: To test the false positive rate of CellDMC and TCA using ALSPAC DNAm data and random phenotypes

## Date: 2021-06-07

## NOTE: need to use R version > 4 to install TCA (not gcc 9.10) - module add languages/R-4.0.3-bioconductor-gcc7.1.0

## pkgs
# BiocManager::install("EpiDISH")
# install.packages("TCA", dependencies = TRUE)
library(tidyverse) # tidy code and data
library(EpiDISH) # CellDMC 
library(TCA) # TCA
# library(usefunc) # own package of useful functions

## 
args <- commandArgs(trailingOnly = TRUE) 
split <- as.numeric(args[1])
# split <- 1
split1 <- (split - 1) * 10 + 1
split2 <- split * 10
message("split = ", split1, " to ", split2)

# ---------------------------------------------------------------
# Setup
# ---------------------------------------------------------------

## Samplesheet
samplesheet_file <- "/panfs/panasas01/sscm/ms13525/aries-release-v4/data/samplesheet/data.Robj"
load(samplesheet_file)
# just mums
samplesheet <- samplesheet[samplesheet$time_point == "FOM", ]

## DNAm data
aries_meth_file <- "/panfs/panasas01/sscm/ms13525/aries-release-v4/data/betas/data.Robj"
load(aries_meth_file)
meth <- beta[, samplesheet$Sample_Name]
rm(beta)

## Cell counts
cc_file="/panfs/panasas01/sscm/ms13525/aries-release-v4/data/derived/cellcounts/blood gse35069/data.txt"
cell_counts <- read_tsv(cc_file) %>%
	dplyr::filter(IID %in% colnames(meth))

phenotypes <- map_dfc(1:1000, function(x) {set.seed(x); sample(c(0,1), ncol(meth), replace=T)})
colnames(phenotypes) <- paste0("p", 1:1000)
p_to_keep <- paste0("p", split1:split2)
phenotypes <- dplyr::select(phenotypes, one_of(p_to_keep))

phenotypes <- phenotypes %>%
	mutate(Sample_Name = colnames(meth)) %>%
	# left_join(cell_counts, by = c("Sample_Name" = "IID")) %>%
	dplyr::select(Sample_Name, everything())

stopifnot(all(cell_counts$IID == phenotypes$Sample_Name))

cell_counts2 <- as.matrix(dplyr::select(cell_counts, -IID))
rownames(cell_counts2) <- cell_counts$IID
cell_counts2[sign(cell_counts2) == -1] <- 0

# ---------------------------------------------------------------
# tests
# ---------------------------------------------------------------

## Need to re-estimate cell props so they all add up to 1 for TCA to work...
reest_cell_props <- function(cc_mat)
{
	out_mat <- map_dfr(1:nrow(cc_mat), function(x) {
		cc_mat[x,]/sum(cc_mat[x,])
	})
	return(as.matrix(out_mat))
}

cell_counts_test <- cell_counts2[1:50, ]
meth_test <- meth[,1:50]
phenotypes_test <- phenotypes[1:50, ]

cell_counts_test <- reest_cell_props(cell_counts_test)
rownames(cell_counts_test) <- colnames(meth_test)

## CellDMC
start_time <- proc.time()
celldmc_res <- lapply(p_to_keep, function(phen) {
	print(phen)
	CellDMC(meth_test, phenotypes_test[[phen]], cell_counts_test)
})
names(celldmc_res) <- p_to_keep
time_taken <- proc.time() - start_time
time_taken # ~30 mins

## TCA
start_time <- proc.time()
tca_res <- lapply(p_to_keep, function(phen) {
	print(phen)
	tca_phen <- matrix(phenotypes_test[[phen]])
	colnames(tca_phen) <- phen
	rownames(tca_phen) <- phenotypes_test$Sample_Name
	out <- tca(X = meth_test, C1 = tca_phen, W = cell_counts_test)
	return(out)
})
names(tca_res) <- p_to_keep
time_taken <- proc.time() - start_time
time_taken # ~100 mins

# ---------------------------------------------------------------
# output results
# ---------------------------------------------------------------



## lambda
get_lambda <- function(pvals) {
	pvals <- na.omit(pvals)
	obs <- qchisq(pvals, df=1, lower.tail = FALSE)
	lambda <- median(obs)/qchisq(0.5, df = 1)
	boot.medians <- sapply(1:100, function(i) median(sample(obs, replace = T)))
	se <- sd(boot.medians / qchisq(0.5, df = 1))
	out <- tibble(lambda_est = lambda, lambda_se = se)
	return(out)
	# median(qchisq(pvals, df = 1, lower.tail = F), na.rm = T) / qchisq(0.5, 1)
}

cell_types <- names(celldmc_res[[p_to_keep[1]]]$coe)

celldmc_lamb <- lapply(p_to_keep, function(phen) {
	res <- celldmc_res[[phen]]
	out <- map_dfr(cell_types, function(ct) {
		p <- res$coe[[ct]]$p
		out <- get_lambda(p)
		out$n_sig <- sum(p < 1e-7)
		return(out)
	}) %>% 
		mutate(cell = cell_types)
	return(out)
})
names(celldmc_lamb) <- p_to_keep

tca_lamb <- lapply(p_to_keep, function(phen) {
	res <- tca_res[[phen]]
	out <- map_dfr(cell_types, function(ct) {
		column <- paste0(ct, ".", phen)
		p <- res$gammas_hat_pvals[, column]
		out <- get_lambda(p)
		out$n_sig <- sum(p < 1e-7)
		return(out)
	}) %>% 
		mutate(cell = cell_types)
	return(out)
})
names(tca_lamb) <- p_to_keep


## save it all! 
out_dir <- "results/temp"
out_file <- paste0("celldmc_res_split", split, ".RData")
save(celldmc_lamb, file = file.path(out_dir, out_file))

out_file <- paste0("tca_res_split", split, ".RData")
save(tca_lamb, file = file.path(out_dir, out_file))

print("FIN!")




















