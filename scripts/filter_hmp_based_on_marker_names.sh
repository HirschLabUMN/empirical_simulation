#!/bin/bash
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=150gb
#SBATCH -J filter_hmp_based_on_marker_names
#SBATCH -o /home/hirschc1/burns756/empirical_sim/eo_files/%x_%j.out
#SBATCH -e /home/hirschc1/burns756/empirical_sim/eo_files/%x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=burns756@umn.edu
#SBATCH --no-requeue

module load R/3.6.0

# go to project folder
cd ~/empirical_sim

# filter hapmap
run_pipeline.pl -Xmx150g -importGuess ${IN} \
                -includeSiteNamesInFile ${LIST} \
                -export ${OUT} \
                -exportType HapmapDiploid
