#!/bin/bash
#SBATCH --time=5:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=120gb
#SBATCH -J get_markers_high-ld_sig-hit
#SBATCH -o /home/hirschc1/della028/projects/genomic_prediction/hybrids/analysis/MSI_dump/%x_%j.out
#SBATCH -e /home/hirschc1/della028/projects/genomic_prediction/hybrids/analysis/MSI_dump/%x_%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=della028@umn.edu
#SBATCH --no-requeue

module load R/3.6.0

# go to project folder
cd ~/projects/genomic_prediction/hybrids

# plot ld
Rscript scripts/get_markers_high-ld_sig-hit.R ${LDFILE} ${SIGHIT} ${FOLDER}
