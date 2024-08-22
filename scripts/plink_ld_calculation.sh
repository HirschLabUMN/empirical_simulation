#!/bin/bash
#SBATCH -J plink_ld_calculation
#SBATCH -o /home/hirschc1/della028/projects/genomic_prediction/hybrids/analysis/MSI_dump/%x_%j.out
#SBATCH -e /home/hirschc1/della028/projects/genomic_prediction/hybrids/analysis/MSI_dump/%x_%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=della028@umn.edu
#SBATCH --no-requeue

# transform % into whitespace of additional options
OPTS=$(echo $OPTS | tr "%" " ")

# go to project folder
cd ${TMPFOLDER}

# calculate LD
plink --file ${IN} --make-founders --r2 gz dprime with-freqs --ld-window-r2 0 --ld-window ${WINDOW}000 --ld-window-kb ${WINDOW} --geno ${FILTER} --out ${OUT} ${OPTS}
