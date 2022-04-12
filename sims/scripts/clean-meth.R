# -------------------------------------------------------
# Filter ARIES betas
# -------------------------------------------------------
# Version = v4

# -------------------------------------------------------
# Setup
# -------------------------------------------------------

## pkgs
library(tidyverse) # tidy code and data
library(meffil) # contains 450k features
library(usefunc) # own package of useful functions - devtools::install_github("thomasbattram/usefunc")

## args
args <- commandArgs(trailingOnly = TRUE)
meth_file <- args[1] 
samplesheet_file <- args[2]
detp_file <- args[3]
outfile <- args[4]

# ARIES_DATA_DIR = ""
# meth_file = paste0(ARIES_DATA_DIR, "/methylation/aries-v4/aries-meth-v4.Robj")
# samplesheet_file = paste0(ARIES_DATA_DIR, "/methylation/aries-v4/aries-samplesheet.Robj")
# detp_file = paste0(ARIES_DATA_DIR, "/methylation/aries-v4/aries-detp-v4.Robj")
# outfile = "data/aries_fom.RData"

# -------------------------------------------------------
# Load ARIES data and remove bad samples
# -------------------------------------------------------

## load aries data
samplesheet <- new_load(samplesheet_file)
samplesheet <- samplesheet %>%
	dplyr::filter(time_point == "FOM")

# samples to remove
sample_rm <- which(samplesheet$duplicate.rm == "Remove" | samplesheet$genotypeQCkids == "ETHNICITY" | 
				   samplesheet$genotypeQCkids == "HZT;ETHNICITY" | samplesheet$genotypeQCmums == "/strat")
samplesheet <- samplesheet[-sample_rm, ]

# methylation data 
beta <- new_load(meth_file)
meth <- beta[, samplesheet$Sample_Name]
rm(beta)

# detection p values
detp <- new_load(detp_file)
pvals <- detp[, samplesheet$Sample_Name]
rm(detp)

print("finished reading in stuff")

# -------------------------------------------------------
# Remove bad probes
# -------------------------------------------------------

## load annotation data
annotation <- meffil.get.features("450k")

## Filter meth data (remove sex chromosomes and SNPs and probes with high detection P-values)
pvalue_over_0.05 <- pvals > 0.05
count_over_0.05 <- rowSums(sign(pvalue_over_0.05))
Probes_to_exclude_Pvalue <- rownames(pvals)[which(count_over_0.05 > ncol(pvals) * 0.05)]
XY <- as.character(annotation$name[which(annotation$chromosome %in% c("chrX", "chrY"))])
SNPs.and.controls <- as.character(annotation$name[-grep("cg|ch", annotation$name)])
annotation<- annotation[-which(annotation$name %in% c(XY, SNPs.and.controls, Probes_to_exclude_Pvalue)), ]
print(length(annotation))
print(dim(meth))
meth <- base::subset(meth, row.names(meth) %in% annotation$name)
paste("There are now ", nrow(meth), " probes")
paste(length(XY), "were removed because they were XY")
paste(length(SNPs.and.controls), "were removed because they were SNPs/controls")
paste(length(Probes_to_exclude_Pvalue), "were removed because they had a high detection P-value")
rm(XY, SNPs.and.controls, pvals, count_over_0.05, pvalue_over_0.05, Probes_to_exclude_Pvalue)

filtered_vars <- c("detection_p_values", "on_XY", "SNPs/controls")

# COULD ALSO ADD ZHOU LIST HERE! 

# -------------------------------------------------------
# Filter ARIES betas
# -------------------------------------------------------

q <- rowQuantiles(meth, probs = c(0.25, 0.75), na.rm = T)
iqr <- q[, 2] - q[, 1]
too.hi <- which(meth > q[,2] + 3 * iqr, arr.ind=T)
too.lo <- which(meth < q[,1] - 3 * iqr, arr.ind=T)
if (nrow(too.hi) > 0) meth[too.hi] <- NA
if (nrow(too.lo) > 0) meth[too.lo] <- NA

dim(meth)
num_na <- apply(meth, 2, function(x){sum(is.na(x))})
rem_samp <- which(num_na > (0.05 * nrow(meth)))
meth <- meth[, -rem_samp]
dim(meth)

print(paste0("Number of samples removed = ", length(rem_samp)))

## impute matrix -- REQUIRED FOR CELL-SPECIFIC METHODS
impute_matrix <- function(x, FUN = function(x) matrixStats::rowMedians(x, na.rm = T)) {
    idx <- which(is.na(x), arr.ind = T)
    if (length(idx) > 0) {
        v <- FUN(x)
        v[which(is.na(v))] <- FUN(matrix(v, nrow = 1))
        x[idx] <- v[idx[, "row"]]
    }
    return(x)
}

meth <- impute_matrix(meth)

save(meth, file = outfile)

