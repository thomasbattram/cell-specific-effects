# -------------------------------------------------------
# Script to generate surrogate variables
# -------------------------------------------------------

## Aim: generate surrogate variables for each trait in the EWAS

## Date: 2022-01-17

## pkgs
library(tidyverse) # tidy data and code
library(sva) # calculating surrogate variables
library(SmartSVA) # calculating SVs
library(matrixStats) # imputing DNAm data
library(aries) # to access the samplesheet 
library(usefunc) # personal package of useful functions

args <- commandArgs(trailingOnly = TRUE)
phen_file <- args[1]
meth_file <- args[2]
meta_file <- args[3]
pcs_file <- args[4]
aries_dir <- args[5]
out_file <- args[6]
removed_out_file <- args[7]

# phen_file <- "data/aries-fom-phenotype-data.tsv"
# meth_file <- "../sims/data/aries_fom.RData"
# meta_file <- "data/metadata.tsv"
# pcs_file <- "data/FOM_pcs.eigenvec"
# aries_dir <- ""
# out_file <- "data/svs/fm1ms100.tsv"
# removed_out_file <- "data/svs/removed/fm1ms100.RData"

## Get trait name from outfile
get_trait <- function(out_file)
{
	file <- basename(out_file)
	trait <- gsub("\\..*", "", file)
	return(trait)
}

trait <- get_trait(out_file)

## read in data
pheno_dat <- read_tsv(phen_file) %>%
	dplyr::select(aln, alnqlet, qlet, one_of(trait))
meta_dat <- read_tsv(meta_file) %>%
	dplyr::filter(alspac_name == trait)
meth <- new_load(meth_file)
aries <- aries.select(aries_dir, time.point = "FOM")
# samplesheet 
samplesheet <- new_load("/user/work/tb13101/cell-specific-effects/aries-ewas/data/aries-samplesheet.Robj")
samplesheet <- samplesheet %>%
	dplyr::filter(Sample_Name %in% colnames(meth)) %>%
	as_tibble

pcs <- read_delim(pcs_file, delim = " ", col_names = F) 
head(pcs)
colnames(pcs) <- c("FID", "IID", paste0(rep("PC", times = 20), 1:20))
pcs$ALN <- gsub("[A-Z]", "", pcs[["FID"]])
pcs <- dplyr::select(pcs, -IID, -FID)
pcs_to_keep <- paste0(rep("PC", times = 10), 1:10)
pc_cols <- pcs_to_keep


# -------------------------------------------------------
# functions
# -------------------------------------------------------

# impute function from Matt
impute_matrix <- function(x, FUN = function(x) rowMedians(x, na.rm = T)) {
    idx <- which(is.na(x), arr.ind = T)
    if (length(idx) > 0) {
        v <- FUN(x)
        v[which(is.na(v))] <- FUN(matrix(v, nrow = 1))
        x[idx] <- v[idx[, "row"]]
    }
    return(x)
}

# function to add quotes for weird trait names
addq <- function(x) paste0("`", x, "`")

generate_svs <- function(trait, phen_data, meth_data, covariates = "", nsv, 
						 IID = "Sample_Name") {
	print("Starting SV generation")
	covs <- covariates[!grepl("^sv", covariates)]
	phen <- phen_data %>%
		dplyr::select(one_of(IID), one_of(trait, covs)) %>%
		.[complete.cases(.), ]
	
	mdat <- meth_data[, colnames(meth_data) %in% phen[[IID]]]
	phen <- phen %>%
		dplyr::filter(!!as.symbol(IID) %in% colnames(mdat))
	
	# models 
	trait_mod <- paste0("~ ", addq(trait))
	cov_mod <- paste(covs, collapse = " + ")
	if (covs != "") {
		full_mod <- paste(trait_mod, cov_mod, sep = " + ")
		fom <- as.formula(full_mod)
		# null model
		fom0 <- as.formula(paste0("~ ", cov_mod))
		mod0 <- model.matrix(fom0, data = phen)
	} else {
		fom <- as.formula(trait_mod)
		mod0 <- NULL
	}

	# full model - with variables of interest 
	mod <- model.matrix(fom, data = phen)

	# Estimate the surrogate variables
	tryCatch({
		svobj <- smartsva.cpp(mdat, mod, mod0, n.sv = nsv, VERBOSE = T)
		svs <- as.data.frame(svobj$sv, stringsAsFactors = F)
		svs[[IID]] <- phen[[IID]]
		# head(svs)
		colnames(svs)[1:nsv] <- paste0("sv", 1:nsv)
		return(svs)
	}, error = function(e) {err_msg(e, r_msg = TRUE, user_msg = paste("SV fail"))})
}

# -------------------------------------------------------
# sort data for generating SVs
# -------------------------------------------------------

# methylation data
mdata <- impute_matrix(meth)

covs <- c("age", pc_cols)
phen_dat <- pheno_dat %>%
	mutate(aln = as.character(aln)) %>%
	left_join(pcs, by = c("aln" = "ALN")) %>%
	left_join(samplesheet, by = c("aln" = "ALN")) %>%
	dplyr::select(Sample_Name, aln, alnqlet, qlet, all_of(c(trait, covs))) %>%
	dplyr::filter(!is.na(Sample_Name))

# -------------------------------------------------------
# Generate SVs
# -------------------------------------------------------

svs <- generate_svs(trait = trait, 
					phen_data = phen_dat, 
					meth_data = mdata, 
					covariates = covs, 
					nsv = 20, 
					IID = "Sample_Name")

# -------------------------------------------------------
# Check association between SVs and phenotype of interest
# -------------------------------------------------------

sv_nam <- grep("sv", colnames(svs), value=T)

sv_check_dat <- phen_dat %>%
	left_join(svs) %>%
	dplyr::select(one_of(trait, sv_nam)) %>%
	na.omit

sv_assoc <- map_dfr(sv_nam, function(x) {
	print(x)
	form <- as.formula(paste(x, trait, sep = " ~ "))
	fit <- lm(form, data = sv_check_dat)
	out_nums <- summary(fit)$coef[2, ]
	out <- as_tibble(t(as.matrix(out_nums))) %>%
		mutate(sv = x) %>%
		dplyr::select(sv, beta = Estimate, SE = `Std. Error`, t_val = `t value`, P = `Pr(>|t|)`)
	return(out)
})

## remove associations at P<0.01 (changing from P<0.01 as doing 20 tests here...)
sv_to_rm <- sv_assoc %>% 
	dplyr::filter(P < 0.01) %>%
	pull(sv)

svs_out <- svs %>% 
	dplyr::select(-one_of(sv_to_rm))

# -------------------------------------------------------
# Write out results and info about removed SVs
# -------------------------------------------------------

write.table(svs_out, file = out_file,
			sep = "\t", quote = F, col.names = T, row.names = F)

sv_rem_info <- lapply(sv_to_rm, function(x) {svs[[x]]})
names(sv_rem_info) <- sv_to_rm
sv_rem_info$sv_assoc <- sv_assoc

save(sv_rem_info, file = removed_out_file)
