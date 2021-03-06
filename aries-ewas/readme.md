# ARIES EWAS

This section of the project looks to run cell-specific EWAS for real phenotypes within ARIES and then compare the results for these EWAS with those that were conducted as part of The EWAS Catalog.

Broad steps:

1. Extract ALSPAC data
2. Run cell-specific EWAS
3. Compare cell-specific results with normal EWAS

## Extract ALSPAC data

1. Move phenotype meta data from The EWAS Catalog RDSF to `data/ec-aries-metadata.tsv`
2. Extract ALSPAC phenotype data using [`pheno-data-extraction.R`](scripts/pheno-data-extraction.R)
3. Move that data to the RDSF and to BluePebble

__NOTE:__ Make sure no ALSPAC data is saved locally!

## Run cell-specific EWAS

Snakemake workflow:

0. Start a tmux session - `tmux new -s cell-spec`
1. Activate conda env - `conda activate /home/user/tb13101/conda-envs/envs/snakemake`
2. Edit the Snakefile template and name it "Snakefile"
3. Do a dry run with `snakemake -nrp`
4. Submit pipeline as a job:
`
snakemake -prk \
-j 100 \
--cluster-config bp1-cluster.json \
--cluster "sbatch \
  --job-name={cluster.name} \
  --nodes={cluster.nodes} \
  --ntasks-per-node={cluster.ntask} \
  --cpus-per-task={cluster.ncpu} \
  --time={cluster.time} \
  --mem={cluster.mem} \
  --output={cluster.output} \
  --error={cluster.error}"
`
5. Deactivate the tmux session (CTRL+b + d)
6. Kill the tmux session when jobs are complete - `tmux kill-session -t cell-spec`