# ---------------------------------------------------------------
# Test CellDMC and TCA false positive rate
# ---------------------------------------------------------------

## Aim: To test the false positive rate of CellDMC and TCA using ALSPAC DNAm data and random phenotypes

## Date: 2021-06-07

## NOTE: need to use R version > 4 to install TCA (not gcc 9.10) - module add languages/R-4.0.3-bioconductor-gcc7.1.0

## pkgs
# BiocManager::install("EpiDISH")
# install.packages("TCA", dependencies = TRUE)
# install.packages("omicwas")
# BiocManager::install("TOAST") 
library(tidyverse) # tidy code and data
library(EpiDISH) # CellDMC 
library(TCA) # TCA
library(omicwas) # omicwas
library(TOAST) # TOAST
library(aries) # getting aries cell counts
library(usefunc) # own package of useful functions - devtools::install_github("thomasbattram/usefunc")

## args
args <- commandArgs(trailingOnly = TRUE) 
meth_file <- args[1]
cc_file <- args[2]
outfile <- args[3]

# setwd("SCRATCH_SPACE")
# meth_file <- "data/aries_fom.RData" 
# cc_file <- "ARIES_DATA_DIR/aries-blood-cell-counts.txt" 
# aries_dir <- "/user/work/ms13525/aries"
# outfile <- "/home/tb13101/projects/cell-specific-effects/results/temp/tcareg_res_split70.RData"

## data
methods <- c("celldmc", "tca", "tcareg", "omicwas", "toast")
get_method <- function(outfile) {
	gsub(".*results/temp/", "", outfile) %>%
		gsub("_res_split.*", "", .)
}
method <- get_method(outfile)

get_split <- function(outfile) as.numeric(str_extract(outfile, "[0-9]"))
split <- get_split(outfile)
split1 <- (split - 1) * 10 + 1
split2 <- split * 10
message("split = ", split1, " to ", split2)

aries <- aries.select(aries_dir, time.point = "FOM")
cell_counts <- aries$cell.counts[["blood-gse35069"]]

# ---------------------------------------------------------------
# Setup
# ---------------------------------------------------------------

## DNAm data
meth <- new_load(meth_file)

## Cell counts
cell_counts <- cell_counts[rownames(cell_counts) %in% colnames(meth), ]
# cell_counts <- read_tsv(cc_file) %>%
# 	dplyr::filter(IID %in% colnames(meth))

phenotypes <- map_dfc(1:1000, function(x) {set.seed(x); sample(c(0,1), ncol(meth), replace=T)})
p_names <- paste0("p", 1:1000)
colnames(phenotypes) <- p_names
p_to_keep <- paste0("p", split1:split2)
phenotypes <- dplyr::select(phenotypes, one_of(p_to_keep))

phenotypes <- phenotypes %>%
	mutate(Sample_Name = colnames(meth)) %>%
	# left_join(cell_counts, by = c("Sample_Name" = "IID")) %>%
	dplyr::select(Sample_Name, everything())

stopifnot(all(rownames(cell_counts) == phenotypes$Sample_Name))

cell_counts[sign(cell_counts) == -1] <- 0

## Need to re-estimate cell props so they all add up to 1 for TCA to work...
reest_cell_props <- function(cc_mat)
{
	out_mat <- map_dfr(1:nrow(cc_mat), function(x) {
		cc_mat[x,]/sum(cc_mat[x,])
	})
	return(as.matrix(out_mat))
}

cell_counts_test <- cell_counts[1:50, ]
meth_test <- meth[,1:50]
phenotypes_test <- phenotypes[1:50, ]

cell_counts_test <- reest_cell_props(cell_counts_test)
rownames(cell_counts_test) <- colnames(meth_test)

# meth_test <- meth_test[1:50, 1:50]

# ---------------------------------------------------------------
# tests
# ---------------------------------------------------------------

## CellDMC function
run_celldmc <- function(p_to_keep)
{
	out <- lapply(p_to_keep, function(phen) {
		print(phen)
		CellDMC(meth_test, phenotypes_test[[phen]], cell_counts_test)
	})
	names(out) <- p_to_keep
	return(out)
}

## TCA function
run_tca <- function(p_to_keep)
{
	out <- lapply(p_to_keep, function(phen) {
		print(phen)
		tca_phen <- matrix(phenotypes_test[[phen]])
		colnames(tca_phen) <- phen
		rownames(tca_phen) <- phenotypes_test$Sample_Name
		out <- tca(X = meth_test, C1 = tca_phen, W = cell_counts_test)
		return(out)
	})
	names(out) <- p_to_keep
	return(out)
}

## TCA-reg function
run_tcareg <- function(p_to_keep)
{
	tca_mdl <- tca(X = meth_test, W = cell_counts_test, constrain_mu = TRUE)
	out <- lapply(p_to_keep, function(phen) {
		print(phen)
		tca_phen <- matrix(phenotypes_test[[phen]])
		colnames(tca_phen) <- phen
		rownames(tca_phen) <- phenotypes_test$Sample_Name
		out <- tcareg(X = meth_test, tca.mdl = tca_mdl, y = tca_phen, test = "marginal_conditional")
		return(out)
	})
	names(out) <- p_to_keep
	return(out)
}

## omicwas function
run_omicwas <- function(p_to_keep)
{
	out <- lapply(p_to_keep, function(phen) {
		seed <- which(p_names %in% phen)
		set.seed(seed)
		Y <- meth_test[sample(1:nrow(meth_test), 1000), ]
		omicwas_phen <- matrix(phenotypes_test[[phen]])
		colnames(omicwas_phen) <- phen
		rownames(omicwas_phen) <- phenotypes_test$Sample_Name
		res <- ctassoc(X = omicwas_phen, W = cell_counts_test, Y = Y, 
					   test = "nls.logit", regularize = TRUE)
		out <- res$coefficients %>% dplyr::filter(term == phen)
		return(out)
	})
	names(out) <- p_to_keep
	return(out)	
}

## toast function
run_toast <- function(p_to_keep)
{
	out <- lapply(p_to_keep, function(phen) {
		print(phen)
		toast_phen <- data.frame(pheno = phenotypes_test[[phen]])
		rownames(toast_phen) <- phenotypes_test$Sample_Name
		Design_out <- makeDesign(toast_phen, cell_counts_test)
		fitted_model <- fitModel(Design_out, meth_test)
		res <- csTest(fitted_model, coef = "pheno", 
		                    cell_type = NULL, contrast_matrix = NULL)
		return(res)
	})
	names(out) <- p_to_keep
	return(out)
}

function_name <- paste0("run_", method)
sim_func <- match.fun(function_name)
# p_to_keep <- p_to_keep[1]
start_time <- proc.time()
res <- sim_func(p_to_keep)
time_taken <- proc.time() - start_time
time_taken

# ---------------------------------------------------------------
# output results
# ---------------------------------------------------------------

## extract results
extract_results <- function(p_to_keep, method, res, sig_val = 1e-7)
{
	cell_types <- colnames(cell_counts_test)
	lamb_out <- lapply(p_to_keep, function(phen) {
		p_res <- res[[phen]]
		out <- map_dfr(cell_types, function(ct) {
			if (method == "celldmc") {
				p <- p_res$coe[[ct]]$p
			} else if (method == "tca") {
				column <- paste0(ct, ".", phen)
				p <- p_res$gammas_hat_pvals[, column]
			} else if (method == "tcareg") {
				p <- p_res$pvals[, ct]
			} else if (method == "omicwas") {
				p <- p_res[p_res$celltype == ct, "p.value", drop = TRUE]
				sig_val <- 0.05/1000
			} else if (method == "toast") {
				p <- p_res[[ct]]$p_value
			}
			out <- get_lambda(p)
			out$n_sig <- sum(p < sig_val)
			return(out)
		})	%>% 
			mutate(cell = cell_types)
		return(out)	
	})
	names(lamb_out) <- p_to_keep
	return(lamb_out)
}


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

out_res <- extract_results(p_to_keep, method, res)

## save it all! 
save(out_res, file = outfile)

print("FIN!")















