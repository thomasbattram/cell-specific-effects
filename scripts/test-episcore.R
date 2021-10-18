# ----------------------------------------------------------------------
# Testing EPISCORE
# ----------------------------------------------------------------------

## Aim: to test EPISCORE to see if it's possible to derive a DNAm reference matrix with just single-cell 
## 		RNA-seq data

## Date: 2021-07-19


## pkgs
# remotes::install_github("aet21/EpiSCORE")
library(tidyverse) # tidy code and data
library(EpiSCORE) # generate DNAm ref 
library(usefunc) # own package of useful functions

## Way that episcore works:
## 	1. Construction of a tissue-specific mRNA expression reference matrix
##  2. Identification of imputable genes (just using roadmap or Stem-Cell-Matrix Compendium-2)
##  3. Imputation of DNAm reference matrix from mRNA expression matrix

# ----------------------------------------------------------------------
# Test following EpiSCORE vignette - https://github.com/aet21/EpiSCORE
# ----------------------------------------------------------------------

data(lungSS2mca1) ##loads in scRNA-Seq SmartSeq2 lung atlas
data(lung10Xmca1) ##loads in a subset of the scRNA-Seq 10X lung atlas for validation purpose
ls()

ncpct.v <- summary(factor(celltypeSS2.idx));
names(ncpct.v) <- celltypeSS2.v;
print(ncpct.v);

## Constructing the expression reference matrix
# Requires expression matrix, index vector that indicates cell type of each column of matrix, a vector of 
# the cell names in the index vector 
expref.o <- ConstExpRef(lungSS2mca1.m, celltypeSS2.idx, celltypeSS2.v, markspecTH = rep(3, 4));
print(dim(expref.o$ref$med));
head(expref.o$ref$med);
str(expref.o)

## Constructing DNAm reference matrix

## Both reference datasets (RMAP and SCM2) are used to create separate references and then merged to 
## increase marker gene coverage THIS IS ONLY BECAUSE for the overlapping genes in the two references
## there is very strong correlation in DNAm patterns
refMscm2.m <- ImputeDNAmRef(expref.o$ref$med,db="SCM2",geneID="SYMBOL");
refMrmap.m <- ImputeDNAmRef(expref.o$ref$med,db="RMAP",geneID="SYMBOL");

refMmg.m <- ConstMergedDNAmRef(refMscm2.m,refMrmap.m);
print(dim(refMmg.m));
head(refMmg.m);

