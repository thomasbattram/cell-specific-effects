# ------------------------------------------------------------
# Examining the results of the test methods
# ------------------------------------------------------------

## Aim: examine results of simulations to assess false positives for CellDMC and TCA

## Date: 2021-06-11

## pkgs
library(tidyverse) # tidy code and data
library(cowplot) # arrange multiple plots neatly on same page
library(usefunc) # own package of useful functions

## data
res_dir <- "results/temp"
res_files <- list.files(res_dir, full.names = T)

tca_files <- grep("tca", res_files, value = T)

tca_res <- map_dfr(tca_files, function(file) {
	res <- new_load(file)
	return(bind_rows(res, .id = "phenotype"))
})

celldmc_files <- grep("cell", res_files, value = T)

celldmc_res <- map_dfr(celldmc_files, function(file) {
	res <- new_load(file)
	return(bind_rows(res, .id = "phenotype"))
})

# TEST 2!!
# tca_res <- new_load("results/tca_res_50cpgs.RData")
# tca_res <- bind_rows(tca_res, .id = "phenotype")
# celldmc_res <- new_load("results/celldmc_res_50cpgs.RData")
# celldmc_res <- bind_rows(celldmc_res, .id = "phenotype")

summary(tca_res)
summary(celldmc_res)

tca_res$method <- "TCA"
celldmc_res$method <- "celldmc"
all_res <- bind_rows(tca_res, celldmc_res)

## plots
celltypes <- unique(all_res$cell)

pal <- get_cb_palette()[1:2]
fdr_plot_list <- lapply(celltypes, function(ct) {
	ct_res <- all_res[all_res$cell == ct, ]
	ggplot(ct_res, aes(x = phenotype, y = n_sig, colour = method)) +
		geom_point() +
		scale_colour_manual(values = pal) + 
		labs(y = "False positives") +
		theme_bw() + 
		theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks = element_blank()) 
})
names(fdr_plot_list) <- celltypes

leg <- get_legend(fdr_plot_list[[1]])

fdr_plot_list <- lapply(fdr_plot_list, function(p) p + theme(legend.position = "none"))

prow <- plot_grid(plotlist = fdr_plot_list, nrow = 3, 
				  labels = celltypes)

out_plots <- plot_grid(prow, leg, ncol = 2, rel_widths = c(3, .4))
ggsave("results/celldmc-tca-false-positives.pdf", plot = out_plots)
# ggsave("results/celldmc-tca-false-positives-50cpgs.pdf", plot = out_plots)
