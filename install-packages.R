## Installing R packages needed for analyses

args <- commandArgs(trailingOnly = TRUE)
folders <- args[1]
folders <- unlist(strsplit(folders, " "))

pkgs <- list(cran = c("devtools",
                      "tidyverse", 
                      "BiocManager",   
                      "bookdown", 
                      "knitr", 
                      "kableExtra",
                      "cowplot", 
                      "RColorBrewer",  
                      "SmartSVA", 
                      "matrixStats",
                      "haven", 
                      "TCA", 
                      "omicwas"),
             bioc = c("sva", "EpiDISH", "TOAST"), 
             git = c("https://github.com/perishky/meffil", 
                     "https://github.com/perishky/ewaff", 
                     "https://github.com/explodecomputer/alspac",
                     "https://github.com/MRCIEU/aries", 
                     "https://github.com/thomasbattram/usefunc"))

for (pkg in pkgs$cran) {
  cat("R package:", pkg, "\n")
  installed <- installed.packages()[,"Package"]
  if (!pkg %in% installed)
     install.packages(pkg, repos = "https://www.stats.bris.ac.uk/R/")
}

for (pkg in pkgs$bioc) {
  cat("R package:", pkg, "\n")
  installed <- installed.packages()[,"Package"]
  if (!pkg %in% installed)
    BiocManager::install(pkg)
}

for (url in pkgs$git) {
  installed <- installed.packages()[,"Package"]
  pkg <- basename(url)
  cat("R package:", pkg, "\n")
  if (!pkg %in% installed)
    devtools::install_github(url)
}

## Checking dependencies
message("\n\n ***Checking whether packages used in scripts are installed.*** \n\n")

if (!require(renv)) {
  install.packages("renv")
}
library(renv)

dep <- dependencies(folders)
pkgs <- unique(dep$Package)
installed <- installed.packages()[,"Package"]
for (x in pkgs) {
  if (!x %in% installed) {
    warning(x, " NOT INSTALLED.")
  }
}

message("\n\nFinished checks. \n\n")