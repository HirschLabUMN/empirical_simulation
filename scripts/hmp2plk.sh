#!/bin/bash
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=50gb
#SBATCH -J hmp2plk
#SBATCH -o /home/hirschc1/della028/projects/genomic_prediction/hybrids/analysis/MSI_dump/%x_%j.out
#SBATCH -e /home/hirschc1/della028/projects/genomic_prediction/hybrids/analysis/MSI_dump/%x_%j.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=della028@umn.edu
#SBATCH --no-requeue

# go to project folder
cd ~/projects/genomic_prediction/hybrids

# transform hmp into plink format
run_pipeline.pl -Xmx110g -importGuess ${IN} -export ${OUT} -exportType Plink
