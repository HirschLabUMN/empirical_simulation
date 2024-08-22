#!/bin/bash
#SBATCH --time=8:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --mem=20gb
#SBATCH -J anova_sim_traits_hybrids_cor-curve
#SBATCH -o /home/hirschc1/burns756/empirical_sim/eo_files/%x_%j.out
#SBATCH -e /home/hirschc1/burns756/empirical_sim/eo_files/%x_%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=burns756@umn.edu
#SBATCH --no-requeue

module load R/4.2.2-openblas 

# go to project folder
cd ~/empirical_sim

for nmarkers in 1 3 5 10 20 30 40 50 60 70 80 90 100 150 200 250 300 350; do

  echo "${METHOD} - ${nmarkers} - ${EFFECT}"

  # define results folder for specific scenario
  RESULTS=${FOLDER}/traits/${METHOD}/n_markers_${nmarkers}/effects_${EFFECT}

  # check if sim file exists first
  FILE=$(echo ${RESULTS}/Simulated_Data*)

  if [[ -e ${FILE} ]]; then
    # run anova
    Rscript scripts/anova_sim_traits_hybrids.R ${RESULTS}
  else
    echo "  skip"
  fi

done
