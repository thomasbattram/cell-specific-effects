#!/bin/bash

#SBATCH --job-name=gen-svs
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1:00:00
#SBATCH --mem=16GB
#SBATCH --array=1-389

## Scratch space
SCRATCH_WD="/user/work/tb13101/cell-specific-effects/aries-ewas"
WD="/user/home/tb13101/projects/cell-specific-effects/aries-ewas"
cd ${WD}

## input files for analyses
meta_file="${SCRATCH_WD}""/data/metadata.tsv"
# script="${WD}""/scripts/gen-svs.R"
pheno="${SCRATCH_WD}""/data/aries-fom-phenotype-data.tsv"
meth="${SCRATCH_WD}""/../sims/data/aries_fom.RData"
pcs="${SCRATCH_WD}""/data/FOM_pcs.eigenvec"


## Test time taken to run analyses (just run all in R)
## Add taking a trait for each array number! -- should be able to take from a text file of traits??
	# Test with code below
trait_n="${SLURM_ARRAY_TASK_ID}"
# trait_n="1"
trait=`sed "${trait_n}q;d" ${SCRATCH_WD}/data/traits.txt`
echo $trait

## output files
results="${SCRATCH_WD}""/data/svs/${trait}.tsv"
removed="${SCRATCH_WD}""/data/svs/removed/${trait}.RData"

## Echo input and output files
echo "Here are the input files:"
echo $meta_file
echo $pheno
echo $meth
echo $pcs
echo "Here are the output files:"
echo $results
echo $removed

## Run the script
Rscript scripts/gen-svs.R "${pheno}" "${meth}" "${meta_file}" "${pcs}" "" "${results}" "${removed}"