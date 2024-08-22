#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=110gb
#SBATCH -J create_hybrid_genotypes
#SBATCH -o /home/hirschc1/burns756/empirical_sim/eo_files/%x_%j.out
#SBATCH -e /home/hirschc1/burns756/empirical_sim/eo_files/%x_%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=burns756@umn.edu
#SBATCH --no-requeue

module load R/3.6.0

# go to project folder
cd ~/empirical_sim

Rscript scripts/create_hybrid_genotypes.R ${HYBINFO} ${RILGENO} ${OUT} --het-parents=${HET}
