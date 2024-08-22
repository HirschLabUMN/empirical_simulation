#!/bin/bash
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --mem=20gb
#SBATCH -J anova_real_traits_hybrids
#SBATCH -o /home/hirschc1/burns756/empirical_sim/eo_files/%x_%j.out
#SBATCH -e /home/hirschc1/burns756/empirical_sim/eo_files/%x_%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=burns756@umn.edu
#SBATCH --no-requeue

module load R/4.2.2-openblas

# go to project folder
cd ~/empirical_sim

# run anova
Rscript scripts/anova_real_traits_hybrids.R ${PHENO}
