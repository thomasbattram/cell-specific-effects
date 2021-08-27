# ------------------------------------------------------------
# Examining the results of the test methods
# ------------------------------------------------------------

## Aim: examine results of simulations to assess false positives for CellDMC and TCA

## Date: 2021-06-11

## pkgs
library(tidyverse) # tidy code and data
library(cowplot) # arrange multiple plots neatly on same page
library(usefunc) # own package of useful functions

## args
args <- commandArgs(trailingOnly = TRUE)
temp_results <- args[1]
fp_outfile <- args[2]
lamb_outfile <- args[3]
fp_lamb_outfile <- args[4]

# temp_results <- c("results/temp/toast_res_split10.RData results/temp/toast_res_split11.RData results/temp/celldmc_res_split10.RData results/temp/celldmc_res_split11.RData results/temp/tca_res_split10.RData results/temp/tca_res_split11.RData results/temp/tcareg_res_split10.RData results/temp/tcareg_res_split11.RData results/temp/omicwas_res_split10.RData results/temp/omicwas_res_split11.RData")
# fp_outfile <- "results/false-positives-plot.pdf"
# lamb_outfile <- "results/lambdas-plot.pdf"
# fp_lamb_outfile <- "results/false-positives-and-lambdas-plot.pdf"

# ------------------------------------------------------------
# Put data together from temporary files
# ------------------------------------------------------------

methods <- c("tca", "celldmc", "tcareg", "toast", "omicwas") 
bind_data <- function(res_files, method)
{
	filenames <- grep(paste0(method, "_res"), res_files, value = T)
	message("There are ", length(filenames), " files for the ", method, " method.")
	out_res <- map_dfr(filenames, function(file) {
		res <- new_load(file)
		return(bind_rows(res, .id = "phenotype"))
	})
	return(out_res)
}

files <- unlist(str_split(temp_results, " "))
list_res <- lapply(methods, bind_data, res_files = files)
names(list_res) <- methods
all_res <- bind_rows(list_res, .id = "method")

# ------------------------------------------------------------
# Plot results
# ------------------------------------------------------------
celltypes <- unique(all_res$cell)

## False positives

pal <- get_cb_palette()[1:length(unique(all_res$method))]
fdr_plot_list <- lapply(celltypes, function(ct) {
	ct_res <- all_res[all_res$cell == ct, ]
	ggplot(ct_res, aes(x = phenotype, y = n_sig, colour = method)) +
		geom_point() +
		scale_colour_manual(values = pal) + 
		labs(y = "False positives", title = ct) +
		theme_classic() + 
		theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks = element_blank()) 
})
names(fdr_plot_list) <- celltypes

leg <- get_legend(fdr_plot_list[[1]])

fdr_plot_list <- lapply(fdr_plot_list, function(p) p + theme(legend.position = "none"))

# ggsave("results/testplot.pdf", plot = fdr_plot_list[[1]])

prow <- plot_grid(plotlist = fdr_plot_list, nrow = length(celltypes)/2)

out_plots <- plot_grid(prow, leg, ncol = 2, rel_widths = c(3, .4))
ggsave(fp_outfile, plot = out_plots)
# ggsave("results/celldmc-tca-false-positives-50cpgs.pdf", plot = out_plots)

## lambda plots

pal <- get_cb_palette()[1:length(unique(all_res$method))]
p_lamb <- ggplot(all_res, aes(x = cell, y = lambda_est, fill = method)) + 
	geom_violin() + 
	geom_boxplot(position=position_dodge(0.9), width = 0.2) + 
	# stat_summary(fun=median, geom="point", size=2) + 
	labs(y = "Lambda") + 
	scale_fill_manual(values = pal) + 
	theme_bw()
ggsave(lamb_outfile, plot = p_lamb)

## Do large lambdas correspond to high false positive rate?

arrange(all_res, desc(n_sig))

pal <- get_cb_palette()[1:length(unique(all_res$method))]
falselamb_plot_list <- lapply(celltypes, function(ct) {
	ct_res <- all_res[all_res$cell == ct, ]
	ggplot(ct_res, aes(x = lambda_est, y = n_sig, colour = method)) +
		geom_point() +
		scale_colour_manual(values = pal) + 
		labs(y = "False positives", title = ct) +
		theme_classic() + 
		theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks = element_blank()) 
})
names(falselamb_plot_list) <- celltypes

leg <- get_legend(falselamb_plot_list[[1]])

falselamb_plot_list <- lapply(falselamb_plot_list, function(p) p + theme(legend.position = "none"))

prow <- plot_grid(plotlist = falselamb_plot_list, nrow = length(celltypes)/2)

out_plots <- plot_grid(prow, leg, ncol = 2, rel_widths = c(3, .4))
ggsave(fp_lamb_outfile, plot = out_plots)















