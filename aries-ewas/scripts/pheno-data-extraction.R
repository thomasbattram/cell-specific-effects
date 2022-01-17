# -------------------------------------------------------
# Extracting ALSPAC data
# -------------------------------------------------------

## need IDs for people in ARIES

## Note on ALSPAC dates: FOM1 clinics started in 2008 and were completed by 2011. 
## T questionaires are 18 years post-birth (returned in 2012). S questionnaires are done 12 years 1 month 
## post birth (returned in 2007)

## Date: 2022-01-13

## pkgs
library(tidyverse) # tidy code and data
library(alspac) # to extract ALSPAC data
library(usefunc) # own package of useful functions

## data
ec_meta <- read_tsv("data/ec-aries-metadata.tsv")
IDs <- read_tsv("data/ARIES_sample_ids.txt") # THIS WAS TAKEN FROM THE EWAS CATALOG RDSF SPACE

alspac_data_dir <- "/Volumes/Data/"
setDataDir(alspac_data_dir)

data(current)
# data(useful)

# -------------------------------------------------------
# RUN THIS TO UPDATE THE DICTIONARIES
# -------------------------------------------------------
# current <- createDictionary("Current", name="current")
# useful <- createDictionary("Useful_data", name="useful")

# -------------------------------------------------------
# Filter out data not present enough in ARIES
# -------------------------------------------------------

# Read in the ARIES IDs and extract ones from timepoint of interest
IDs <- dplyr::filter(IDs, time_point == "FOM")
str(IDs)

# ------------------------------------------------------------------------------------
# Extract data 
# ------------------------------------------------------------------------------------

## Check ALSPAC variables in EC have the same name as they did before
which(!ec_meta$alspac_name %in% current$name) # 0!!

## ECZEMA VARS
current[grep("eczema", current$lab),]
# s1026 = A2q: Mother has had eczema in last 2 years (asked near to FOM1)
# s4140 = D1o: Mother has taken medication for eczema in the past 12 months (asked near to FOM1)
## SMOKING VARS
current[grep("smok", current$lab),]
# t5560 = E71: Respondent ever smoked in the past (4175 counts)
# r6010 = G2a: Respondent has ever been a smoker (7595 counts) -- CHECK DIFFERENCE IN NUMBERS!

## Smoking test
smoking_current <- current[current$name %in% c("t5560", "r6010"), ]
smoking_result <- extractVars(smoking_current)
smoking_result %>% 
	dplyr::filter(aln %in% IDs$ALN) %>% 
	pull(t5560) %>% 
	table
smoking_result %>% 
	dplyr::filter(aln %in% IDs$ALN) %>% 
	pull(r6010) %>% 
	table
# r0610 has many more people, so going with that phenotype.
vars_of_interest <- c(ec_meta$alspac_name, "s1026", "r6010")

new_current <- current %>%
	dplyr::filter(name %in% vars_of_interest)

# paths of interest
# PoI <- c()

# extraction
result <- extractVars(new_current)

# ------------------------------------------------------------------------------------
# Initial look at data 
# ------------------------------------------------------------------------------------
## finding the age of participants at each questionnaire

res <- result %>%
	dplyr::filter(aln %in% IDs$ALN) 

# check for aln and qlet columns.
grep("aln|qlet", colnames(res), value = T)
# Change if more than just aln, alnqlet, qlet
dim(res)
dim(new_current)
# delete extra columns minus the aln, qlet and alnqlet columns
col_rm <- colnames(res)[!colnames(res) %in% new_current$name]
col_rm <- col_rm[!col_rm %in% c("aln", "qlet", "alnqlet")]
res <- res[, !colnames(res) %in% col_rm]
qlet_cols <- grep("qlet", colnames(res), value = T)

meta_dat <- new_current %>%
	dplyr::select(obj, alspac_name = name, unedited_label = lab) %>%
	as_tibble

# ------------------------------------------------------------------------------------
# Start cleaning data
# ------------------------------------------------------------------------------------

extract_alspac_labels <- function(dat, alsp_dir)
{
	if (!all(c("path", "obj", "name") %in% colnames(dat))) stop("'dat' requires the columns 'path', 'obj', 'name' in columns")
	if (!file.exists(alsp_dir)) stop("ALSPAC directory not found.")
	
	uniq_obj <- unique(dat$obj)
	message("Extracting labels in ", length(uniq_obj), " ALSPAC files.")
	list_out <- lapply(uniq_obj, function(ob) {
		df <- dat[dat$obj == ob, ]
		full_path <- file.path(alsp_dir, unique(df$path), unique(df$obj))
		message("Reading in file: ", full_path)
		alsp_res <- haven::read_dta(full_path) %>%
			dplyr::select(all_of(df$name))
		val_labs <- lapply(df$name, function(nam) {labelled::val_labels(alsp_res[[nam]])})
		names(val_labs) <- df$name
		return(val_labs)
	})
	out <- unlist(list_out, recursive = FALSE)

	return(out)
}

## get labels for categorical variables
cat_vars <- c("r6010", "s1026")

lab_dat_in <- new_current
val_labs <- extract_alspac_labels(lab_dat_in, alsp_dir = alspac_data_dir)

unique_labs <- unique(unlist(val_labs))
which(sapply(val_labs, function(x) -4 %in% x))
which(sapply(val_labs, function(x) -11 %in% x))
which(sapply(val_labs, function(x) -110 %in% x))
which(sapply(val_labs, function(x) -2 %in% x))
which(sapply(val_labs, function(x) -3 %in% x))
which(sapply(val_labs, function(x) 0 %in% x))
which(sapply(val_labs, function(x) 99 %in% x))

## sort out values - e.g. remove missing values
val_labs
num_vars <- names(val_labs)[!names(val_labs) %in% cat_vars]
clean_res <- res %>%
	mutate(r6010 = case_when(r6010 == 1 ~ "yes", r6010 == 2 ~ "no"),
		   s1026 = case_when(s1026 %in% c(1, 2) ~ "yes", s1026 == 3 ~ "no"))

missing_vals <- c(-9999, -11, -10, -8, -1, -111, -110, -108, -101, -2, -3, -4, 0)
test_df <- data.frame(x = c(1, missing_vals), y = c(2, missing_vals))
for (i in missing_vals) {
	test_df[test_df == i] <- NA
}
for (i in missing_vals) {
	clean_res[clean_res == i] <- NA
}

## sanity check
table(clean_res$s1026)
table(clean_res$r6010)

summary(clean_res) # doesn't seem to be any negative values
## Check to see if there are negative values
neg_vals <- sapply(clean_res, function(x) {
	vals <- x[!is.na(x)]
	out <- try(sign(vals))
	if (inherits(out, "try-error")) return(NULL)
	any(out == -1)
})
neg_val_vars <- names(which(unlist(neg_vals)))
## do they make sense? 
meta_dat %>%
	dplyr::filter(alspac_name %in% neg_val_vars)

## proinsulin there...
clean_res$proinsulin_FOM1
pro_in_dat <- haven::read_dta(file.path(alspac_data_dir, "Current/Other/Samples/Mother", "Mother_samples_5b.dta"))
pro_in_dat <- pro_in_dat %>%
	dplyr::select(proinsulin_FOM1)

str(pro_in_dat) ## Seems as though it's just not labelled properly, almost certainly it's a missing value!
clean_res$proinsulin_FOM1 <- ifelse(clean_res$proinsulin_FOM1 == -99, NA, clean_res$proinsulin_FOM1)

# --------------------------------------------------------------
# remove phenotypes with too much missing data
# --------------------------------------------------------------
missing_dat <- map_df(seq_along(clean_res), function(x) {
	out <- data.frame(phen = colnames(clean_res[x]), na_count = sum(is.na(clean_res[[x]])))
	return(out)
})
sum(missing_dat$na_count > nrow(clean_res)/2) # 0 phenotypes have over 50% missing data

# ------------------------------------------------------------------------------------
# Finishing tidying data + saving it all
# ------------------------------------------------------------------------------------

meta_dat$n <- map_dbl(meta_dat$alspac_name, function(x) sum(!is.na(clean_res[[x]])))
meta_dat$binary <- map_lgl(meta_dat$alspac_name, function(x) is.binary(clean_res[[x]]))

## Write out meta-data
write.table(meta_dat, "data/metadata.tsv", 
			row.names = F, col.names = T, quote = F, sep = "\t")

## Write out phenotype data
save_alsp_data <- function(outdat, filename, outpath, password)
{
	ori_wd <- getwd()
	setwd(outpath)
	write.table(outdat, file = filename, 
				quote = F, col.names = T, row.names = F, sep = "\t")
	zip(gsub(".tsv", ".zip", filename), 
		files = filename, 
		flags = paste("--password", password))
	system(paste("rm", filename))
	setwd(ori_wd)
}

# Set new password each time
PASSWORD <- "" ## REMEMBER THIS!
save_alsp_data(clean_res, "aries-fom-phenotype-data.tsv", "data", PASSWORD)
