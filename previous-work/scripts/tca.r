# installing tca
library(devtools)
install_github("cozygene/TCA")

# load tca
library(TCA)	

samples <- data.matrix(read.csv('samplesheet_colremoved.csv', row.names = 1))
meth_data <- data.matrix(read.csv('RAmeth_matrix_linesremoved_10Kmostvar.csv', row.names = 1))
pheno <- t(samples)
k <- 6

# MUST convert to matrix
ra <- matrix(pheno[,1, drop=FALSE])
covars <- matrix(pheno[,-1, drop=FALSE])


# houseman cell proportions from GLINT 
cell_proportions <- data.matrix(read.csv('cell_proportions_houseman.csv', row.names=1))

# use of tca
tca_model <- tca(meth_data, cell_proportions, C1=covars, C2=ra, refit_W=TRUE)

# write tca output to file
saveRDS(tca_model, file = "RA_TCAmodels.rds")

# re-load tca output
tca_model <-  readRDS(file="RA_TCAmodels.rds")

# fitting a tca regression model

tca_reg <- tcareg(meth_data, tca_model, ra, C3 = NULL, test = "marginal",
      null_model = NULL, alternative_model = NULL, save_results = TRUE,
      output = "TCA", sort_results = TRUE, parallel = FALSE,
      num_cores = NULL, log_file = "TCA.log", features_metadata = NULL,
      debug = TRUE)

# save tca reg as r object
saveRDS(tca_reg, file="RA_TCAreg.rds")


tca_tensor <- tensor(meth_data, tca_model, parallel = FALSE,
                     log_file = "TCA.log", debug = FALSE)

saveRDS(tca_tensor, file="RA_TCAtensor.rds")

# write TCA tensor output to file to be parsed using python script
tca_output <- readRDS(file="RA_TCAtensor.rds")
tca_reg <- readRDS(file="RA_TCAreg.rds")

ct1 <- tca_output[[1]]
ct2 <- tca_output[[2]]
ct3 <- tca_output[[3]]
ct4 <- tca_output[[4]]
ct5 <- tca_output[[5]]
ct6 <- tca_output[[6]]
ct7 <- tca_output[[7]]
write.csv(ct1, file = "ct1.csv")
write.csv(ct2, file = "ct2.csv")
write.csv(ct3, file = "ct3.csv")
write.csv(ct4, file = "ct4.csv")
write.csv(ct5, file = "ct5.csv")
write.csv(ct6, file = "ct6.csv")

pvalct1 <- tca_reg[[1]][[9]]
pvalct2 <- tca_reg[[2]][[9]]
pvalct3 <- tca_reg[[3]][[9]]
pvalct4 <- tca_reg[[4]][[9]]
pvalct5 <- tca_reg[[5]][[9]]
pvalct6 <- tca_reg[[6]][[9]]
write.csv(pvalct1, file = "pvalct1.csv")
write.csv(pvalct2, file = "pvalct2.csv")
write.csv(pvalct3, file = "pvalct3.csv")
write.csv(pvalct4, file = "pvalct4.csv")
write.csv(pvalct5, file = "pvalct5.csv")
write.csv(pvalct6, file = "pvalct6.csv")

