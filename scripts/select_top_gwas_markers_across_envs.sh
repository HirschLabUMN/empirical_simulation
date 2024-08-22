#!/bin/bash
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=50gb
#SBATCH -J select_top_gwas_markers_across_envs
#SBATCH -o /home/hirschc1/burns756/empirical_sim/eo_files/%x_%j.out
#SBATCH -e /home/hirschc1/burns756/empirical_sim/eo_files/%x_%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=burns756@umn.edu
#SBATCH --no-requeue


module load R/3.6.0

# go to project folder
cd ~/empirical_sim

for nmarkers in 1 3 5 10 20 30 40 50 60 70 80 90 100 200 300 400 500 600 700 800 900 1000; do
  # select markers
  Rscript scripts/select_top_gwas_markers_across_envs.R \
          ${GWASFILE} \
          ${LDFILE} \
          ${nmarkers} \
          ${FOLDER} \
          --ld-threshold=${LD} \
          --avg-rank-markers=all
done
