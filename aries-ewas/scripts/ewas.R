# ----------------------------------------
# ewas script
# ----------------------------------------

## Aim: Run EWAS using cell-specific effects methods across all traits.

## Date: 2022-01-18

## pkgs
library(tidyverse) # tidy code and data
library(ewaff) # for EWAS functions
library(EpiDISH) # celldmc
library(TCA) # TCA
library(omicwas) # omicwas
library(TOAST) # TOAST
library(aries) # for extracting aries samplesheet
library(usefunc) # own package of useful functions

args <- commandArgs(trailingOnly = TRUE)
phen_file <- args[1]
meth_file <- args[2]
meta_file <- args[3]
svs_file <- args[4]
pcs_file <- args[5]
aries_dir <- args[6]
out_files <- args[7]
failed_outfile <- args[8]

# setwd("SCRATCH_SPACE")
# phen_file <- "data/aries-fom-phenotype-data.tsv"
# meth_file <- "../sims/data/aries_fom.RData"
# meta_file <- "data/metadata.tsv"
# svs_file <- "data/svs/r6010.tsv"
# ## svs_file <- "data/svs/LDLD_FOM1.tsv"
# pcs_file <- "data/FOM_pcs.eigenvec"
# aries_dir <- "/user/work/ms13525/aries"
# out_files <- "results/ewas-res/celldmc/r6010.RData results/ewas-res/tca/r6010.RData results/ewas-res/tcareg/r6010.RData results/ewas-res/toast/r6010.RData results/ewas-res/omicwas/r6010.RData"
# ## out_files <- "results/ewas-res/celldmc/LDLD_FOM1.RData results/ewas-res/tca/LDLD_FOM1.RData results/ewas-res/tcareg/LDLD_FOM1.RData results/ewas-res/toast/LDLD_FOM1.RData results/ewas-res/ewaff/LDLD_FOM1.RData"
# ## out_files <- "results/ewas-res/r6010.tsv"
# failed_outfile <- "results/ewas-res/failed-ewas.tsv"

out_files <- unlist(str_split(out_files, " "))
# if (length(out_files) != 5) {
#     message(out_files)
#     stop("Length of out_files is not 5!")
# }

## Get trait name from outfile
get_trait <- function(out_file)
{
    file <- basename(out_file)
    trait <- gsub("\\..*", "", file)
    return(trait)
}

trait <- get_trait(out_files[1])

## read in data
pheno_dat <- read_tsv(phen_file) %>%
    dplyr::select(aln, alnqlet, qlet, one_of(trait))
meta_dat <- read_tsv(meta_file) %>%
    dplyr::filter(alspac_name == trait)
meth <- new_load(meth_file)
aries <- aries.select(aries_dir, time.point = "FOM")
aries$samples$ALN <- as.character(aries$samples$ALN)

pcs <- read_delim(pcs_file, delim = " ", col_names = F) 
head(pcs)
colnames(pcs) <- c("FID", "IID", paste0(rep("PC", times = 20), 1:20))
pcs$ALN <- gsub("[A-Z]", "", pcs[["FID"]])
pcs <- dplyr::select(pcs, -IID, -FID)
pcs_to_keep <- paste0(rep("PC", times = 10), 1:10)
pc_cols <- pcs_to_keep

svs <- read_tsv(svs_file)

# -------------------------------------------------------
# sort data for running EWAS
# -------------------------------------------------------

# methylation data
# mdata <- impute_matrix(meth)

covs <- c("age", pc_cols, grep("sv", colnames(svs), value = T))
phen_dat <- pheno_dat %>%
    mutate(aln = as.character(aln)) %>%
    left_join(pcs, by = c("aln" = "ALN")) %>%
    left_join(aries$samples, by = c("aln" = "ALN")) %>%
    left_join(svs) %>%
    dplyr::select(Sample_Name, aln, all_of(c(trait, covs))) %>%
    dplyr::filter(!is.na(Sample_Name)) %>%
    na.omit(.)

cell_counts <- aries$cell.counts[["blood-gse35069"]]
cell_counts <- cell_counts[rownames(cell_counts) %in% phen_dat$Sample_Name, ]

idx <- match(phen_dat$Sample_Name, rownames(cell_counts))
cell_counts <- cell_counts[idx, ]
stopifnot(all(rownames(cell_counts) == phen_dat$Sample_Name))

## Negative values are soooo close to zero, and values need to be 0-1 for some methods
cell_counts[sign(cell_counts) == -1] <- 0

## Need to re-estimate cell props so they all add up to 1 for TCA to work...
reest_cell_props <- function(cc_mat)
{
    out_mat <- map_dfr(1:nrow(cc_mat), function(x) {
        cc_mat[x,]/sum(cc_mat[x,])
    })
    return(as.matrix(out_mat))
}

rownam <- rownames(cell_counts)
cell_counts <- reest_cell_props(cell_counts)
rownames(cell_counts) <- rownam

# ----------------------------------------
# EWAS functions
# ----------------------------------------

run_celldmc <- function(temp_phen, temp_meth, phen, cc, covs, IID)
{
    cov_form <- as.formula(paste0("~ ", paste(covs, collapse = " + ")))
    cov_mat <- model.matrix(cov_form, data = temp_phen)
    out <- CellDMC(beta.m = temp_meth, pheno.v = temp_phen[[phen]], frac.m = cc, cov.mod = cov_mat)
    return(out)
}

sort_celldmc <- function(celldmc_res) 
{
    celldmc_coefs <- celldmc_res$coe
    cols_to_get <- c("Estimate", "SE", "p")
    celldmc_out <- lapply(cols_to_get, function(x) {
        out <- map_dfr(celldmc_coefs, x) %>%
            as.data.frame
        rownames(out) <- rownames(celldmc_coefs[[1]])
        return(out)
    })
    names(celldmc_out) <- c("beta", "se", "p")
    return(celldmc_out)
}

run_tca <- function(temp_phen, temp_meth, phen, cc, covs, IID)
{
    phen_vals <- temp_phen[[phen]]
    if (is.binary(phen_vals)) phen_vals <- ifelse(phen_vals == "yes", 1, 0)
    tca_phen <- matrix(phen_vals)
    colnames(tca_phen) <- phen
    rownames(tca_phen) <- temp_phen[[IID]]
    tca_covs <- as.matrix(temp_phen[, !colnames(temp_phen) %in% c(phen, IID, "aln")])
    rownames(tca_covs) <- temp_phen[[IID]]
    out <- tca(X = temp_meth, C1 = tca_phen, W = cc, C2 = tca_covs)
    return(out)
}

sort_tca <- function(tca_res) 
{
    tca_beta <- tca_res$gammas_hat
    tca_p <- tca_res$gammas_hat_pvals
    tca_out <- list(beta = tca_beta[, grep(trait, colnames(tca_beta))], 
                      p = tca_p[, grep(trait, colnames(tca_p))])
    tca_out <- lapply(tca_out, function(x) {
        colnames(x) <- gsub("\\..*", "", colnames(x))
        return(x)
    })
    return(tca_out)
}

 
run_tcareg <- function(temp_phen, temp_meth, phen, cc, covs, IID)
{
    phen_vals <- temp_phen[[phen]]
    if (is.binary(phen_vals)) phen_vals <- ifelse(phen_vals == "yes", 1, 0)
    tca_phen <- matrix(phen_vals)
    colnames(tca_phen) <- phen
    rownames(tca_phen) <- temp_phen[[IID]]
    tca_covs <- as.matrix(temp_phen[, !colnames(temp_phen) %in% c(IID, phen, "aln")])
    rownames(tca_covs) <- temp_phen[[IID]]    
    tca_mdl <- tca(X = temp_meth, W = cc, C2 = tca_covs, constrain_mu = TRUE)
    out <- tcareg(X = temp_meth, tca.mdl = tca_mdl, y = tca_phen, test = "marginal_conditional")
    return(out)
}

sort_tcareg <- function(tcareg_res)
{
    tcareg_betas <- tcareg_res$beta
    tcareg_p <- tcareg_res$pvals
    tcareg_out <- list(beta = tcareg_res$beta, 
                         p = tcareg_res$pvals)
    return(tcareg_out)
}

# run_omicwas <- function(temp_phen, temp_meth, phen, cc, covs, IID, seed = 2)
# {
#     set.seed(seed)
#     Y <- temp_meth[sample(1:nrow(temp_meth), 1000), ]
#     # Y <- temp_meth
#     phen_vals <- temp_phen[[phen]]
#     if (is.binary(phen_vals)) phen_vals <- ifelse(phen_vals == "yes", 1, 0)
#     omicwas_phen <- matrix(phen_vals)
#     colnames(omicwas_phen) <- phen
#     rownames(omicwas_phen) <- temp_phen[[IID]]
#     omicwas_covs <- as.matrix(temp_phen[, !colnames(temp_phen) %in% c(phen, IID, "aln")])
#     rownames(omicwas_covs) <- temp_phen[[IID]]
#     res <- ctassoc(X = omicwas_phen, W = cc, Y = Y, C = omicwas_covs, 
#                    test = "nls.logit", regularize = TRUE)
#     ## Matrix (or vector) of covariates; samples x covariates. X, W, Y, C should benumeric.
#     out <- res$coefficients %>% dplyr::filter(term == phen)
#     return(out) 
# }

# sort_omicwas <- function(omicwas_res)
# {
#     cols_to_get <- c("estimate", "p.value")
#     omicwas_out <- lapply(cols_to_get, function(x) {
#         out <- omicwas_res %>%
#             dplyr::select(response, celltype, one_of(x)) %>%
#             pivot_wider(names_from = celltype, values_from = one_of(x)) %>%
#             as.data.frame
#         rownames(out) <- out$response
#         out <- out[, !colnames(out) == "response"]
#         return(out)
#     })
#     names(omicwas_out) <- c("beta", "p")
#     return(omicwas_out)
# }

run_toast <- function(temp_phen, temp_meth, phen, cc, covs, IID)
{
    toast_phen <- temp_phen[, !colnames(temp_phen) %in% c(IID, "aln")] %>%
        as.data.frame
    if (is.binary(toast_phen[[phen]])) {
        toast_phen[[phen]] <- ifelse(toast_phen[[phen]] == "yes", 1, 0)
        toast_phen[[phen]] <- as.factor(toast_phen[[phen]])
    }
    rownames(toast_phen) <- temp_phen[[IID]]
    Design_out <- makeDesign(toast_phen, cc)
    fitted_model <- fitModel(Design_out, temp_meth)
    res <- csTest(fitted_model, coef = phen, 
                  cell_type = NULL, contrast_matrix = NULL)
    return(res)
}

sort_toast <- function(toast_res)
{
    toast_res <- toast_res[names(toast_res) != "joint"]
    cols_to_get <- c("beta", "effect_size", "p_value")
    toast_out <- lapply(cols_to_get, function(x) {
        out <- map_dfr(toast_res, x) %>%
            as.data.frame
        rownames(out) <- rownames(toast_res[[1]])
        return(out)
    })
    names(toast_out) <- c("beta", "effect_size", "p")
    return(toast_out)
}

run_ewaff <- function(temp_phen, temp_meth, phen, covs, IID)
{
    model <- as.formula(paste0("methylation ~ ", paste(c(phen, covs), collapse = " + ")))
    out <- ewaff.sites(model, variable.of.interest = phen,
                              methylation = temp_meth, data = temp_phen, method = "glm",
                              generate.confounders = NULL, family = "gaussian")    
    return(out)
}

sort_ewaff <- function(ewaff_res)
{
    list(beta = ewaff_res$table$estimate, 
         se = ewaff_res$table$se, 
         p = ewaff_res$table$p)
}

run_ewas <- function(phen, p_dat, cc, meth_dat, IID, method, covs) 
{
    # Match meth to Pheno
    temp_meth <- meth_dat[, colnames(meth_dat) %in% p_dat[[IID]]]
    temp_meth <- meth_dat[, na.omit(match(p_dat[[IID]], colnames(meth_dat)))]
    temp_phen <- p_dat[match(colnames(temp_meth), p_dat[[IID]]), ]
    temp_cc <- cc[rownames(cc) %in% temp_phen$Sample_Name, ]
    
    if (!all(temp_phen[[IID]] == colnames(temp_meth))) stop("phenotype and DNAm data not matched.")
    ## FOR TESTS
    # temp_meth <- temp_meth[1:50, 1:50]
    # temp_phen <- temp_phen[1:50, ]
    # temp_cc <- temp_cc[1:50, ]

    function_name <- paste0("run_", method)
    ewas_func <- match.fun(function_name)
    # p_to_keep <- p_to_keep[1]
    start_time <- proc.time()
    message("FUNCTION TIME")
    res <- tryCatch({
        ewas_func(temp_phen, temp_meth, phen, temp_cc, covs, IID)
    }, error = function(e) {
        usr_m <- paste0("Error in EWAS of ", phen, " using ", method, ".")
        err_msg(e, r_msg = TRUE, user_msg = usr_m, to_return = phen)        
    })
    time_taken <- proc.time() - start_time
    print(time_taken) # 68 mins (4105.556 seconds) for CellDMC
    return(res)
}

# ----------------------------------------
# Run the EWAS
# ----------------------------------------

get_method <- function(out_file) basename(dirname(out_file))

cell_types <- colnames(cell_counts)

cell_counts_copy <- cell_counts
cell_counts <- cell_counts_copy
# cell_counts <- cell_counts[1:50, ]
meth_copy <- meth
meth <- meth_copy
# meth <- meth[1:50, 1:50]
phen_dat_copy <- phen_dat
# phen_dat <- phen_dat[1:50, ]
phen_dat <- phen_dat_copy

lapply(out_files, function(out_file) {
    method <- get_method(out_file)
    message("Running analyses using ", method)
    ewas_res <- run_ewas(phen = trait, 
                         p_dat = phen_dat, 
                         cc = cell_counts, 
                         meth_dat = meth,
                         IID = "Sample_Name", 
                         method = method, 
                         covs = covs)
    if (is.character(ewas_res)) {
        if (ewas_res == trait) {
            failed_out <- data.frame(phenotype = trait, ewas_method = method)
            write.table(failed_out, file = failed_outfile, 
                        sep = "\t", col.names = F, row.names = F, quote = F, append = T)
            message("Writing out failed EWAS to ", failed_outfile)
            return(NULL)
        }
    }
    message("Finished analyses using ", method)
    sort_function_name <- paste0("sort_", method)
    sort_func <- match.fun(sort_function_name)
    out_res <- sort_func(ewas_res)
    message("Saving sorted results to ", out_file)
    save(out_res, file = out_file)
    rm(list = c("out_res", "ewas_res"))
})

print("FIN")
