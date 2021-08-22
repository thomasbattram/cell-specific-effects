# cell-specific effects

This repo is to store notes and potentially work on the effectiveness of methods that look to ascertain cell-specific DNAm effects in EWAS.

[`methods-list.xlsx`](methods-list.xlsx) contains a list of the different methods with their PMID, paper title, first author, corresponding author + email, and how the method can be implemented.

Notes on the different methods from their respective papers can be found in the [`notes`](notes) folder.

[`previous-work`](previous-work) contains work done by a summer student of Matt's that assessed how well TCA and HIRE tend to estimate actual cell-types. 

The false positive rate of the methods can be tested by running [`test-methods.R`](scripts/test-methods.R) - this requires submission to bc3. To plot the results run through [`test-methods-results.R`](scripts/test-methods-results.R).

## Requirements

### Modules

Please run the below code to add the appropriate modules before running the analyses:

`
module add lang/python/anaconda/3.8.5-2021-Jupyter
# module add apps/pandoc/2.13
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
qsub -I -l select=1:ncpus=1:mem=5G,walltime=1:00:00
module add lang/python/anaconda/3.8.5-2021-Jupyter
cd
conda create -p /home/tb13101/conda-envs
conda install -p /home/tb13101/conda-envs -c conda-forge mamba
conda activate /home/tb13101/conda-envs
mamba create -c conda-forge -c bioconda -n snakemake snakemake
`

To check it's worked run:
`
conda activate /home/tb13101/conda-envs/envs/snakemake
snakemake --help
` 

## Workflow

Snakemake workflow is presented below:

0. Start a tmux session - `tmux new -s cell-spec`, `cd` into correct folder and make the required folders if you need to (INCLUDING "job-errors-and-outputs")
1. Activate conda env - `conda activate /home/tb13101/conda-envs/envs/snakemake`
2. Edit the Snakefile template and name it "Snakefile"
3. Do a dry run with `snakemake -nrp`
4. Submit pipeline as a job:

`
snakemake -pr \
-j 100 \
--cluster "qsub \
	-N cell-spec \
	-l select=1:ncpus=6:mem=32G,walltime=24:00:00 \
	-o job-errors-and-outputs/cell-spec-{rule}-error-{jobid} \
	-e job-errors-and-outputs/cell-spec-{rule}-output-{jobid}"
`

5. Deactivate the tmux session (CTRL+b + d)
6. Check the completion by going back to the tmux session - `tmux a -t cell-spec`
7. __AFTER COMPLETION__ kill the tmux session when jobs are complete - `tmux kill-session -t cell-spec`


## To-do

0. Streamline processes using snakemake (maybe)
1. Test the methods on real phenotypes to see if they tend to output the same results
2. Generate a report - include the previous work in this!