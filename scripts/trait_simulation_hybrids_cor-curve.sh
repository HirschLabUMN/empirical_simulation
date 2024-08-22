#!/bin/bash
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=10gb
#SBATCH -J trait_simulation_hybrids_cor-curve
#SBATCH -o /home/hirschc1/burns756/empirical_sim/eo_files/%x_%j.out
#SBATCH -e /home/hirschc1/burns756/empirical_sim/eo_files/%x_%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=burns756@umn.edu
#SBATCH --no-requeue

module load R/4.2.2-openblas

# check for optional variables -- if not requested, set it blank
[[ -z ${QTNVAR} ]] && QTNVAR="" || QTNVAR="--${QTNVAR}"
[[ -z ${MINORQTLS} ]] && MINORQTLS="" || MINORQTLS="--${MINORQTLS}"
[[ -z ${EFFECTREDALL} ]] && EFFECTREDALL="" || EFFECTREDALL="--${EFFECTREDALL}"

# State the trait
echo ${TRAIT}
echo ${EFFECTRED}

# go to project folder
cd ~/empirical_sim

AVGRANK=all
for nmarkers in 1 3 5 10 20 30 40 50 60 70 80 90 100 150 200 250 300 350; do

  echo "  ${nmarkers}"
  echo ${SEED}
  # get gwas file
  GWASFILE=${FOLDER}/gwas_markers/gwas_top-${nmarkers}_markers.txt

  # create output folder
  OUTFOLDER=${FOLDER}/traits/michael_method/n_markers_${nmarkers}
  mkdir -p ${OUTFOLDER}

  # simulate traits
  Rscript scripts/trait_simulation_hybrids.R ${TRAIT} \
                                             ${GWASFILE} \
                                             ${HMP} ${SVS} ${RESCOR} \
                                             ${H2FILE} ${MEANSFILE} \
                                             ${OUTFOLDER}/effects_${EFFECTRED} \
                                             ${PHENOFILE} \
                                             --effect-size-reduction=${EFFECTRED} \
                                             --seed=${SEED} \
                                             --reps=${REPS} \
                                             --impute-type=${IMPUTETYPE} \
                                             --architecture=pleiotropic \
                                             ${QTNVAR} ${MINORQTLS} ${EFFECTREDALL}
  # change seed
  SEED=$((${SEED}+5))

done
