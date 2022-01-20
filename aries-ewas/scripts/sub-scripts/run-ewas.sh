#!/bin/bash

#SBATCH --job-name=run-ewas
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --mem=16GB
#SBATCH --array=1-100

## Scratch space
SCRATCH_WD="/user/work/tb13101/cell-specific-effects/aries-ewas"
WD="/user/home/tb13101/projects/cell-specific-effects/aries-ewas"
cd ${WD}

trait_n="${SLURM_ARRAY_TASK_ID}"
# trait_n="2"
trait=`sed "${trait_n}q;d" ${SCRATCH_WD}/data/traits.txt`
echo $trait

## input files for analyses
meta_file="${SCRATCH_WD}""/data/metadata.tsv"
# script="${WD}""/scripts/gen-svs.R"
pheno="${SCRATCH_WD}""/data/aries-fom-phenotype-data.tsv"
meth="${SCRATCH_WD}""/../sims/data/aries_fom.RData"
pcs="${SCRATCH_WD}""/data/FOM_pcs.eigenvec"
svs="${SCRATCH_WD}""/data/svs/${trait}.tsv"
aries_dir=""

## output files
RES_DIR="${SCRATCH_WD}""/results/ewas-res/"
results=""${RES_DIR}"omicwas/${trait}.RData "${RES_DIR}"toast/${trait}.RData "${RES_DIR}"tcareg/${trait}.RData "${RES_DIR}"tca/${trait}.RData ${RES_DIR}""celldmc/${trait}.RData"
# got to put methods at end if want to use them first apparently...


## Echo input and output files
echo "Here are the input files:"
echo $meta_file
echo $pheno
echo $meth
echo $svs
echo $pcs
echo $aries_dir
echo "Here are the output files:"
echo $results

## Run the script
time Rscript scripts/ewas.R "${pheno}" "${meth}" "${meta_file}" "${svs}" "${pcs}" "${aries_dir}" "${results}"