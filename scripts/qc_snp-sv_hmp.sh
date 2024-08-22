#!/bin/bash
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
#SBATCH --mem=180gb
#SBATCH -J qc_snp-sv_hmp
#SBATCH -o /home/hirschc1/burns756/empirical_sim/eo_files/%x_%j.out
#SBATCH -e /home/hirschc1/burns756/empirical_sim/eo_files/%x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=burns756@umn.edu
#SBATCH --no-requeue


module load R/3.6.0

# go to project folder
cd ~/emprirical_sim

# create folder to save results
mkdir -p ${FOLDER}
# qc hapmap file
Rscript scripts/qc_snp-sv_hmp.R ${HMP} ${SVS} ${FOLDER}
