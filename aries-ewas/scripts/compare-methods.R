# ----------------------------------------
# Compare EWAS results for methods for estimating cell-specific effects
# ----------------------------------------

## Aim: Compare methods for estimating cell specific-effects

## Date: 2022-01-23

## pkgs
library(tidyverse) # tidy code and data
library(ewaff) # for EWAS functions
library(aries) # for extracting aries samplesheet
library(usefunc) # own package of useful functions

args <- commandArgs(trailingOnly = TRUE)
meta_file <- args[1]
ewas_resdir <- args[2]
traits_file <- args[3]
out_file <- args[4]

# setwd("SCRATCH_SPACE")
# meta_file <- "data/metadata.tsv"
# ewas_resdir <- "results/ewas-res"
# traits_file <- "data/traits.txt"
# ## traits_file <- "data/test-traits.txt"
# out_file <- "results/summary.tsv"

## data
meta_file <- read_tsv(meta_file)
traits <- readLines(traits_file)
names(traits) <- c("BMI", "Smoking")
original_ewas_files <- c("data/DV__BMI__FOM1.txt", "results/ewas-res/r6010.tsv")
names(original_ewas_files) <- c("BMI", "Smoking")
original_ewas <- lapply(original_ewas_files, read_tsv)
names(original_ewas) <- names(original_ewas_files)
## function for reading in data for each trait

methods <- c("celldmc", "tca", "tcareg", "omicwas", "toast")
# methods <- c("celldmc", "tca", "omicwas", "toast")

sort_res <- function(res, cpg = NULL, p_threshold = 1e-5) 
{
	sig_res <- res$p %>%
		as.data.frame() %>%
		rownames_to_column(var = "CpG") %>%
		pivot_longer(!CpG, names_to = "cell_type", values_to = "p") %>%
		dplyr::filter(p < p_threshold)

	if (!is.null(cpg)) {
		sig_res <- sig_res %>%
			dplyr::filter(CpG %in% cpg)
	}

	beta_df <- res$beta %>%
		as.data.frame() %>%
		rownames_to_column(var = "CpG") %>%
		pivot_longer(!CpG, names_to = "cell_type", values_to = "beta")
	
	out_res <- sig_res %>%
		left_join(beta_df) 

	if ("se" %in% names(res)) {
		se_df <- res$se %>%
			as.data.frame() %>%
			rownames_to_column(var = "CpG") %>%
			pivot_longer(!CpG, names_to = "cell_type", values_to = "se")
		out_res <- out_res %>%
			left_join(se_df)
	}

	out_res <- out_res %>%
		dplyr::select(CpG, beta, one_of("se"), p, cell_type) %>%
		as_tibble
	return(out_res)
}

all_res <- lapply(traits, function(trait) {
	out <- lapply(methods, function(method) {
		print(method)
		res_file <- paste0(method, "/", trait, ".RData")
		res_path <- file.path(ewas_resdir, res_file)
		ewas_res <- new_load(res_path)
		out <- sort_res(ewas_res, p_threshold = 1e-7)
		return(out)
	})
	names(out) <- methods
	return(out)
})
names(all_res) <- names(traits)

## To do:
# 1. Get normal EWAS results for BMI from EWAS Catalog res
# 2. Make plots
# 	a. Normal EWAS "hits" vs. their results in specific cells
# 	b. Cell-specific "hits" that aren't in the normal EWAS
# 	c. Complete comparison of effect sizes? - maybe a heatmap with all methods + normal EWAS

## 2a.

get_hits <- function(res, p_threshold) dplyr::filter(res, p.value < p_threshold) %>% arrange(p.value)

## BMI
ori_res <- get_hits(original_ewas[["BMI"]], 1e-7) %>%
	dplyr::select(CpG = probeID, beta = estimate, se, p = p.value) %>%
	mutate(cell_type = "bulk")
cs_res <- lapply(methods, function(method) {
		res_file <- paste0(method, "/", "fm1ms111.RData")
		res_path <- file.path(ewas_resdir, res_file)
		ewas_res <- new_load(res_path)
		out <- sort_res(ewas_res, cpg = ori_res$CpG, p_threshold = 1)
		return(out)
})
names(cs_res) <- methods
cs_res <- cs_res[names(cs_res) != "omicwas"]
cs_res$celldmc <- dplyr::select(cs_res$celldmc, -se)

cs_res_df <- bind_rows(cs_res, .id = "method")

plot_dat <- ori_res %>%
	dplyr::select(-se) %>%
	mutate(method = "bulk") %>%
	bind_rows(cs_res_df) %>%
	mutate(effect_dir = ifelse(sign(beta) == -1, "neg", "pos"), logp = -log10(p))

cb_pal <- get_cb_palette()
cols <- cb_pal[1:2]
names(cols) <- c("neg", "pos")
abcg1_bar <- ggplot(plot_dat[plot_dat$method != "bulk",], aes(x = method, y = logp, fill = as.factor(effect_dir))) + 
	geom_bar(stat="identity") +
	geom_hline(yintercept = plot_dat[plot_dat$method == "bulk", "logp", drop=T], colour = cols[2]) +  
	scale_fill_manual(name = "DoE", values = cols) + 
	labs(y = "-log10(P)") + 
	theme_bw() + 
	theme(axis.text.x = element_text(angle=90)) + 
	facet_grid(~ cell_type)

ggsave("abcg1-comp.png", abcg1_bar)

## Smoking
smok_ori_res <- get_hits(original_ewas[["Smoking"]], 1e-7) %>%
	dplyr::select(CpG = probeID, beta = estimate, se, p = p.value) %>%
	mutate(cell_type = "bulk")
smok_cs_res <- lapply(methods, function(method) {
		res_file <- paste0(method, "/", "r6010.RData")
		res_path <- file.path(ewas_resdir, res_file)
		ewas_res <- new_load(res_path)
		out <- sort_res(ewas_res, cpg = smok_ori_res$CpG, p_threshold = 1)
		return(out)
})
names(smok_cs_res) <- methods
smok_cs_res <- smok_cs_res[names(smok_cs_res) != "omicwas"]

ahrr_res <- lapply(smok_cs_res, function(x) {dplyr::filter(x, CpG == "cg05951221")})

ahrr_res$celldmc <- dplyr::select(ahrr_res$celldmc, -se)

ahrr_res_df <- bind_rows(ahrr_res, .id = "method")

plot_dat <- smok_ori_res %>%
	dplyr::filter(CpG == "cg05951221") %>%
	dplyr::select(-se) %>%
	mutate(method = "bulk") %>%
	bind_rows(ahrr_res_df) %>%
	mutate(effect_dir = ifelse(sign(beta) == -1, "neg", "pos"), logp = -log10(p))


ahhr_bar <- ggplot(plot_dat[plot_dat$method != "bulk",], aes(x = method, y = logp, fill = as.factor(effect_dir))) + 
	geom_bar(stat="identity") +
	geom_hline(yintercept = plot_dat[plot_dat$method == "bulk", "logp", drop=T], colour = cols[1]) +  
	scale_fill_manual(name = "DoE", values = cols) + 
	labs(y = "-log10(P)") + 
	theme_bw() + 
	theme(axis.text.x = element_text(angle=90)) + 
	facet_grid(~ cell_type)

ggsave("ahrr-comp.png", ahhr_bar)


### tabs of both
smok_res_list <- lapply(all_res$Smoking, function(x){ 
	if(nrow(x) == 0) return(NULL)
	return(x)
})

smok_ori_all_res <- original_ewas$Smoking %>%

smok_res_df <- smok_res_list %>%
	bind_rows(.id = "method") %>%
	left_join(original_ewas$Smoking, by = c("CpG" = "probeID")) %>%
	dplyr::select(-Details)

write.table(smok_res_df, file = "smok_res_df_out.tsv", sep = "\t", quote = F, row.names = F, col.names = T)

cpg_of_int <- smok_res_df$

res_file <- paste0(method, "/", trait, ".RData")
res_path <- file.path(ewas_resdir, "tca/r6010.RData")
ewas_res <- new_load(res_path)

temp_out <- sort_res(ewas_res, cpg = cpg_of_int, p_threshold = 1)
write.table(temp_out, file = "smok-tca-cg23248424.tsv", sep = "\t", quote = F, row.names = F, col.names = T)

all_res3 <- lapply(methods, function(method) {
	print(method)
	res_file <- paste0(method, "/r6010.RData")
	res_path <- file.path(ewas_resdir, res_file)
	ewas_res <- new_load(res_path)
	out <- sort_res(ewas_res, p_threshold = 1, cpg = cpg_of_int)
	return(out)
})
names(all_res3) <- methods
all_res3 <- all_res3[names(all_res3) != "omicwas"]
all_res3$celldmc <- dplyr::select(all_res3$celldmc, -se)


all_res3_df <- bind_rows(all_res3, .id = "method")
write.table(all_res3_df, file = "smok-allmeth-cg23248424.tsv", sep = "\t", quote = F, row.names = F, col.names = T)










