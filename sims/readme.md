# sims

This section of the project looks to run cell-specific EWAS for 1000 simulated phenotypes using ARIES DNAm data and assess whether there is likely a high false positive rate for these EWAS. 

Broad steps:

1. Setup DNAm data and generate random phenotypes
2. Run analyses using each method
3. Check false positive rate

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

