#!/bin/bash
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --mem=120gb
#SBATCH -J kinship_VanRaden
#SBATCH -o /home/hirschc1/burns756/empirical_sim/eo_files/%x_%j.out
#SBATCH -e /home/hirschc1/burns756/empirical_sim/eo_files/%x_%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=burns756@umn.edu
#SBATCH --no-requeue


module load R/4.0.4

# go to project folder
cd ~/empirical_sim/

# perform pca
Rscript scripts/kinship_VanRaden.R ${HMP} ${FOLDER}
