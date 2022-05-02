# ----------------------------------------------------------------------
# Testing whether can generate SEs from TCA/TCAREG/TOAST
# ----------------------------------------------------------------------





library(tidyverse) # tidy code and data
library(EpiDISH) # CellDMC 
library(TCA) # TCA
library(omicwas) # omicwas
library(TOAST) # TOAST
library(aries) # getting aries cell counts
library(boot) # bootstrapping
library(usefunc) # own package of useful functions - devtools::install_github("thomasbattram/usefunc")

# setwd("SCRATCH_SPACE")
# meth_file <- "data/aries_fom.RData" 
# aries_dir <- "/user/work/ms13525/aries"
# all_res_output <- "results/bootstrap-tests-res.RData"

## data
methods <- c("tca", "tcareg", "toast")

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

phenotypes <- map_dfc(1:10, function(x) {set.seed(x); sample(c(0,1), ncol(meth), replace=T)})
p_names <- paste0("p", 1:10)
colnames(phenotypes) <- p_names
p_to_keep  <- p_names

phenotypes <- phenotypes %>%
	mutate(Sample_Name = colnames(meth)) %>%
	# left_join(cell_counts, by = c("Sample_Name" = "IID")) %>%
	dplyr::select(Sample_Name, everything())

stopifnot(all(rownames(cell_counts) == phenotypes$Sample_Name))

cell_counts[sign(cell_counts) == -1] <- 0

## Need to re-estimate cell props so they all add up to 1 for TCA to work...
#' @param cc_mat matrix. cell proportions matrix with rows as people and columns as cell types
reest_cell_props <- function(cc_mat)
{
	out_mat <- map_dfr(1:nrow(cc_mat), function(x) {
		cc_mat[x, ] / sum(cc_mat[x, ])
	})
	return(as.matrix(out_mat))
}

cell_counts_test <- cell_counts
meth_test <- meth[1:50, ]
phenotypes_test <- phenotypes

cell_counts_test <- reest_cell_props(cell_counts_test)
rownames(cell_counts_test) <- colnames(meth_test)
cells <- colnames(cell_counts_test)

## OUT TABLE THANG
## Cell type | CpG | beta | P | boot mean | boot median | boot SE |


## METHODS

run_tcareg <- function(pt, mt, cct, phen)
{
	tca_mdl <- tca(X = mt, W = cct, constrain_mu = TRUE)
	tca_phen <- matrix(pt[[phen]])
	colnames(tca_phen) <- phen
	rownames(tca_phen) <- pt$Sample_Name
	out <- tcareg(X = mt, tca.mdl = tca_mdl, y = tca_phen, test = "marginal_conditional")
	return(out)
}

start_time <- proc.time()
tcareg_res <- run_tcareg(phenotypes_test, meth_test, cell_counts_test, "p1")
time_taken <- proc.time() - start_time

sort_tcareg <- function(tcareg_res, cells)
{
    tcareg_betas <- tcareg_res$beta
    tcareg_p <- tcareg_res$pvals
    tcareg_out <- lapply(cells, function(cell) {
    	tibble(cpg = rownames(tcareg_betas), 
    		   beta = tcareg_betas[, grep(cell, colnames(tcareg_betas))], 
    		   p = tcareg_p[, grep(cell, colnames(tcareg_p))])   
    })
    names(tcareg_out) <- cells
    return(tcareg_out)
}

s_tcareg_res <- sort_tcareg(tcareg_res, cells)

## TCA-reg function
boot_tcareg <- function(data, indices, phen, mt, cct)
{
	pt <- data[indices, ]
	mt <- mt[, indices]
	cct <- cct[indices, ]
	tca_mdl <- tca(X = mt, W = cct, constrain_mu = TRUE)
	tca_phen <- matrix(pt[[phen]])
	colnames(tca_phen) <- phen
	rownames(tca_phen) <- pt$Sample_Name
	out <- tcareg(X = mt, tca.mdl = tca_mdl, y = tca_phen, test = "marginal_conditional")
	return(out$beta)
}

output_tcareg <- boot(data = phenotypes_test, statistic = boot_tcareg, R=1000, 
			   phen = "p1", mt = meth_test, cct = cell_counts_test)

## output$t0 = matrix with cell-types as columns and cpgs as rows
## output$t = matrix with rows as replicates and columns as individual values of the output$t0 matrix
##		-- need to figure out how the columns map to the t0 matrix cells 

c1 <- output_tcareg$t0[1,1] ## top left cell of t0
summary(output_tcareg$t[,1]) ## first column
sd(output_tcareg$t[,1])

c2 <- output_tcareg$t0[2,1] ## one row below top left cell of t0
summary(output_tcareg$t[,2]) ## second column of t
sd(output_tcareg$t[,2])

## randomly select 50 and check the dist
# set.seed(2)
# to_select <- sample(1:ncol(output$t), 49)
# plots <- lapply(to_select, function(x) {
# 	plot_dat <- tibble(beta = output$t[,x])
# 	ggplot(plot_dat, aes(x = beta)) + 
# 		geom_histogram(fill = "blue", colour = "black") + 
# 		theme_void()
# })
# cow_out <- cowplot::plot_grid(plotlist = plots, nrow = 7, ncol = 7)
# ggsave("tcareg-beta-dist-test.png", plot = cow_out)


## TCA
run_tca <- function(pt, phen, mt, cct)
{
	tca_phen <- matrix(pt[[phen]])
	colnames(tca_phen) <- phen
	rownames(tca_phen) <- pt$Sample_Name
	out <- tca(X = mt, C1 = tca_phen, W = cct)
	return(out)
}

start_time <- proc.time()
tca_res <- run_tca(phenotypes_test, "p1", meth_test, cell_counts_test)
time_taken <- proc.time() - start_time

sort_tca <- function(tca_res, cells) 
{
    tca_beta <- tca_res$gammas_hat
    tca_p <- tca_res$gammas_hat_pvals
    tca_out <- lapply(cells, function(cell) {
    	tibble(cpg = rownames(tca_beta), 
    		   beta = tca_beta[, grep(cell, colnames(tca_beta))], 
    		   p = tca_p[, grep(cell, colnames(tca_p))])
    }) 
    names(tca_out) <- cells
    return(tca_out)
}

## FOR ALL METHODS OUTPUT TABLE LIKE:
# | CpG | Beta | P |
# (for each cell type)

s_tca_res <- sort_tca(tca_res, cells)

boot_tca <- function(data, indices, phen, mt, cct)
{
	pt <- data[indices, ]
	mt <- mt[, indices]
	cct <- cct[indices, ]
	tca_phen <- matrix(pt[[phen]])
	colnames(tca_phen) <- phen
	rownames(tca_phen) <- pt$Sample_Name
	out <- tca(X = mt, C1 = tca_phen, W = cct)
	return(out$gammas_hat)
}

output_tca <- boot(data = phenotypes_test, statistic = boot_tca, R=1000, 
			   phen = "p1", mt = meth_test, cct = cell_counts_test)

## TOAST
run_toast <- function(pt, phen, mt, cct)
{
	toast_phen <- data.frame(pheno = pt[[phen]])
	rownames(toast_phen) <- pt$Sample_Name
	Design_out <- makeDesign(toast_phen, cct)
	fitted_model <- fitModel(Design_out, mt)
	res <- csTest(fitted_model, coef = "pheno", 
	                    cell_type = NULL, contrast_matrix = NULL)
	return(res)
}

start_time <- proc.time()
toast_res <- run_toast(phenotypes_test, "p1", meth_test, cell_counts_test)
time_taken <- proc.time() - start_time

toast_res <- toast_res[names(toast_res) != "joint"]
cols_to_get <- c("beta", "beta_var", "effect_size", "p_value")
s_toast_res <- lapply(cells, function(x) {
    out <- toast_res[[x]] %>%
        dplyr::select(one_of(cols_to_get)) %>%
        rename(p = p_value) %>%
        rownames_to_column(var = "cpg")
    return(out)
})
names(s_toast_res) <- cells

boot_toast <- function(data, indices, phen, mt, cct, cols_to_get)
{
	pt <- data[indices, ]
	mt <- mt[, indices]
	cct <- cct[indices, ]
	toast_phen <- data.frame(pheno = pt[[phen]])
	# rownames(toast_phen) <- pt$Sample_Name
	Design_out <- makeDesign(toast_phen, cct)
	fitted_model <- fitModel(Design_out, mt)
	res <- csTest(fitted_model, coef = "pheno", 
	                    cell_type = NULL, contrast_matrix = NULL)
	toast_res <- res[names(res) != "joint"]
	# cols_to_get <- "beta"
	out <- map_dfr(toast_res, cols_to_get) %>%
        as.data.frame %>%
        as.matrix
    rownames(out) <- rownames(toast_res[[1]])
	return(out)
}

output_toast <- boot(data = phenotypes_test, statistic = boot_toast, R=1000, 
			   phen = "p1", mt = meth_test, cct = cell_counts_test, cols_to_get = "beta")

output_toast_eff <- boot(data = phenotypes_test, statistic = boot_toast, R=1000, 
			   phen = "p1", mt = meth_test, cct = cell_counts_test, cols_to_get = "effect_size")

# TAB OUT (one per method-cell type combo [18 tables])

# | CpG | B ori | SE ori | P ori | Mean B | Median B | SE | P |

get_p <- function(beta, se) 
{
	z <- beta / se
	p <- 2 * pnorm(-abs(z))
	return(p)
}

get_se <- function(p, beta) 
{
	z <- qnorm(p/2)
	se <- abs(beta / z)
	return(se)
}

cell_res_mapping <- tibble(ct = cells, 
						   max = c(50, 100, 150, 200, 250, 300), 
						   min = c(1, 51, 101, 151, 201, 251))
methods <- c("tca", "tcareg", "toast")
all_out <- lapply(methods, function(method) {
	bres <- get(paste0("output_", method))
	res <- get(paste0("s_", method, "_res"))
	method_out <- lapply(cells, function(cell) {
		crm_temp <- cell_res_mapping[cell_res_mapping$ct == cell, ]
		cols <- seq(crm_temp$min, crm_temp$max)
		cpg_out <- map_dfr(cols, function(x) {
			B <- bres$t[, x]
			tibble(Mean_B = mean(B), 
				   Median_B = median(B), 
				   SE = sd(B), 
				   P = get_p(Mean_B, SE))
		})
		cpg_out$cpg <- rownames(bres$t0)
		comb_out <- res[[cell]] %>%
			left_join(cpg_out) %>%
			mutate(SE_ori = get_se(p, beta)) %>%
			dplyr::select(CpG = cpg, B_ori = beta, SE_ori, P_ori = p, 
						  Mean_B, Median_B, SE, P)
		return(comb_out)
	})
	names(method_out) <- cells
	return(method_out)
})
names(all_out) <- methods

all_res_out <- list(tca = output_tca, tcareg = output_tcareg, toast = output_toast, summary = all_out)
save(all_res_out, file = all_res_output)

## Checking one tab
all_out$tca$Bcell

## Testing comp of SEs
se_dat <- all_out$tca$Bcell %>%
	dplyr::select(CpG, SE_ori, SE) %>%
	pivot_longer(-CpG, names_to = "calculation", values_to = "SE") %>%
	mutate(calculation = ifelse(calculation == "SE", "boot", "no-boot"))

outplot <- ggplot(se_dat, aes(x = calculation, y = SE, group = factor(CpG))) + 
	geom_point() + 
	geom_path(position = "identity")

ggsave("test.png", outplot)

## Need to get outliers

#' Get outliers based on the Tukey method
#' 
#' @param vals numeric vector for which to find outliers
#' @param largest_diff logical. TRUE = extract largest value (differences should be absolute)
#' @return all outliers in the input vector
get_outliers <- function(vals, largest_diff = FALSE)
{
	q <- quantile(vals)
	iqr <- IQR(vals)
	too_hi <- vals > q[4] + 3 * iqr
	if (largest_diff) {
		out <- vals == max(vals)
	} else {
		out <- too_hi
	}
	return(out)
} 

all_outliers <- lapply(methods, function(method) {
	method_out <- lapply(cells, function(cell) {
		temp_out <- all_out[[method]][[cell]]
		temp_out <- temp_out %>%
			mutate(se_diff = abs(SE - SE_ori), 
				   p_diff = abs(P - P_ori)) %>%
			mutate(se_outlier = get_outliers(se_diff, largest_diff = TRUE), 
				   p_outlier = get_outliers(p_diff, largest_diff = TRUE)) %>%
			dplyr::filter(se_outlier|p_outlier)
		return(temp_out)
	})
	names(method_out) <- cells
	return(method_out)
})
names(all_outliers) <- methods

#' Get CpGs from list of outliers from each method and/or cell type
#' 
#' @param outlier_dat list (methods) of lists (cell types) containing data from "all_outliers"
#' @param all logical. TRUE = get all outlier CpGs for all methods and cell types. FALSE = extract CpGs for each method separately. Default = FALSE.
#' @param stat character vector. stats to get outlier CpGs for. Must be "se", "p", or c("se", "p")
#' @return character vector of CpGs that have been identified as outliers 
get_cpg_outliers <- function(outlier_dat, all = FALSE, stat = c("se", "p"))
{
	out <- lapply(outlier_dat, function(d) {
		unlist(map(d, "CpG"))
	})
	if (all) out <- unlist(out)
	return(out)
}

cpg_outliers <- get_cpg_outliers(all_outliers)

# list (cells)
# 

lapply(methods, function(method) {
	ol <- cpg_outliers[[method]]
	bres <- get(paste0("output_", method))
	res <- get(paste0("s_", method, "_res"))
	hists <- lapply(names(ol), function(x) {
		print(x)
		cpg <- ol[x]
		cell <- substr(x, 1, nchar(x) - 1)
		crm_temp <- cell_res_mapping[cell_res_mapping$ct == cell, ]
		cols <- seq(crm_temp$min, crm_temp$max)
		bres_temp <- bres$t[, cols]
		colnames(bres_temp) <- rownames(bres$t0)
		plot_dat <- data.frame(Beta = bres_temp[, cpg])
		hist_plot <- ggplot(plot_dat, aes(x = Beta)) + 
			geom_histogram(fill = "blue", colour = "black") + 
			theme_bw()
		return(hist_plot)
	})		
	test_plots <- cowplot::plot_grid(plotlist = hists, nrow = 6)
	ggsave(paste0("test-hists-", method, ".png"), plot = test_plots)
	return(NULL)
})

lapply(methods, function(method) {
	ol <- cpg_outliers[[method]]
	hists <- lapply(names(ol), function(x) {
		print(x)
		cpg <- ol[x]
		plot_dat <- data.frame(Meth_beta = meth_test[cpg, ])
		hist_plot <- ggplot(plot_dat, aes(x = Meth_beta)) + 
			geom_histogram(fill = "blue", colour = "black") + 
			theme_bw()
		return(hist_plot)
	})		
	test_plots <- cowplot::plot_grid(plotlist = hists, nrow = 6)
	ggsave(paste0("test-hists-meth-", method, ".png"), plot = test_plots)
	return(NULL)
})

## NICE TABLES (ONE PER METHOD)

# | Cell type | CpG | SE (Z) | SE (boot) | SE diff | P | P (boot) | P diff |

nice_tab <- lapply(methods, function(method) {
	outlier_tabs <- all_outliers[[method]]
	out <- map_dfr(outlier_tabs, function(x) {
		x %>%
			dplyr::select(CpG, 
						  `SE (Z)` = SE_ori, `SE (boot)` = SE, `SE diff` = se_diff,
						  `P` = P_ori, `P (boot)` = P, `P diff` = p_diff)
	}, .id = "Cell type")
	out <- out %>%
		tidy_nums
	return(out)
})
names(nice_tab) <- methods

openxlsx::write.xlsx(nice_tab, file = "results/bootstrap-outlier-tabs.xlsx")

## Plotting beta distributions

set.seed(2)
to_select <- sample(1:50, 10)
to_select_ct <- lapply(cells, function(x) {
	crm <- cell_res_mapping[cell_res_mapping$ct == x, ]
	ts <- seq(crm$min, crm$max)[to_select]
	return(ts)
})
names(to_select_ct) <- cells

methods <- c("tca", "tcareg", "toast", "toast_eff")
lapply(methods, function(method) {
	bres <- get(paste0("output_", method))
	plots <- lapply(unlist(to_select_ct), function(x) {
		plot_dat <- tibble(beta = bres$t[,x])
		ggplot(plot_dat, aes(x = beta)) + 
			geom_histogram(fill = "blue", colour = "black") + 
			theme_void()
	})
	cow_out <- cowplot::plot_grid(plotlist = plots, nrow = 10, ncol = 6)
	ggsave(paste0("results/", method, "-beta-dist-test.png"), plot = cow_out)
	return(NULL)
})


### BOOTSTRAP TESTS

test_boot_output <- function(data, indices) {
	out_vals <- sapply(1:4, function(x) rnorm(mean = x, n=1))
	if (all(indices == 1:nrow(data))) out_vals <- 1:4
	out <- matrix(out_vals, nrow=2)
	return(out)
} 
# Performing 1500 replications with boot 
boot_test_res <- boot(data=mtcars, statistic=test_boot_output, R=1500)
boot_test_res$t0
mean(boot_test_res$t[,1])
mean(boot_test_res$t[,2])
mean(boot_test_res$t[,3])
mean(boot_test_res$t[,4])










