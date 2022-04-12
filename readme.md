# cell-specific effects

This repo is to store notes and potentially work on the effectiveness of methods that look to ascertain cell-specific DNAm effects in EWAS.

[`methods-list.xlsx`](docs/methods-list.xlsx) contains a list of the different methods with their PMID, paper title, first author, corresponding author + email, and how the method can be implemented.

Notes on the different methods from their respective papers can be found in [`running-notes.md`](docs/running-notes.md).

[`previous-work`](previous-work) contains work done by a summer student of Matt's that assessed how well TCA and HIRE tend to estimate actual cell-types. 

## Requirements

### Modules

Analyses requires R (version 4.0.3), Python (version 3.8.5) and tmux (version 3.2)

Example code below on how to add the modules needed on Bristol HPC:

`
module add lang/r/4.0.3-bioconductor-gcc
module add lang/python/anaconda/3.8.5-2021-Jupyter
module add tools/tmux/3.2
`

### Conda/snakemake

To check if snakemake is already installed then run 
`
module add lang/python/anaconda/3.8.5-2021-Jupyter
conda env list
`
and if a snakemake environment (ie. "path_to_env/snakemake") is present then it is installed otherwise continue with the instructions below.

To initiate conda on startup (required for activating environments as specified below) then run `module add lang/python/anaconda/3.8.5-2021-Jupyter; conda init bash` (assuming bash is shell of choice), exit bluepebble and then re-enter.

To insall snakemake, run the following code (change out "/home/tb13101" for your own directory):
`
srun --job-name "InteractiveJob" --nodes=1 --ntasks-per-node=1 --cpus-per-task=1 --time=1:00:00 --mem=5GB --pty bash
module add lang/python/anaconda/3.8.5-2021-Jupyter
cd
conda create -p /user/home/tb13101/conda-envs
conda install -p /user/home/tb13101/conda-envs -c conda-forge mamba
conda activate /user/home/tb13101/conda-envs
mamba create -c conda-forge -c bioconda -n snakemake snakemake
`

To check it's worked run:
`
conda activate /user/home/tb13101/conda-envs/envs/snakemake
snakemake --help
` 

### R packages

To install R packages run the [`install-packages.R`](install-packages.R) script and input the folders that need to be checked for R packages, e.g. `Rscript install-packages.R "sims aries-ewas"`. This uses the Bristol mirror for CRAN packages. Please manually change this if Bristol is not the right mirror.