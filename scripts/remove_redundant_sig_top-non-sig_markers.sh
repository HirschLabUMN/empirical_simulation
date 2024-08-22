#!/bin/bash
#SBATCH --time=5:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=50gb
#SBATCH -J remove_redundant_sig_top-non-sig_markers
#SBATCH -o /home/hirschc1/burns756/empirical_sim/eo_files/%x_%j.out
#SBATCH -e /home/hirschc1/burns756/empirical_sim/eo_files/%x_%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=burns756@umn.edu
#SBATCH --no-requeue


module load R/3.6.0

# go to project folder
cd ~/empirical_sim

# create an array of trait names
traits=()
for trait in $(ls -d analysis/gwas/*/); do
  #echo $trait
  traits+=$(basename $trait)
  traits+=" "
done

for trait in ${traits[@]}; do
  for model in additive dominant; do
    echo "${trait} - ${model}"
    # get folder name
    FOLDER=analysis/gwas/${trait}/${model}_model/permutation
    # get input files
    GWAS=${FOLDER}/GAPIT.MLM.BLUE.GWAS.Results.csv
    LDFILE=${FOLDER}/ld_top${NONSIG}_non-sig-markers.ld.gz
    # get pvalue threshold used in gwas
    if [[ -e ${FOLDER}/MLM.BLUE.GWAS.Results.sig-hits.csv ]]; then
      PVAL=$(cut -f 7 -d , ${FOLDER}/MLM.BLUE.GWAS.Results.sig-hits.csv | sed 1d | uniq)
    else
      # if file doesn't exist, just set pval to a very low number
      PVAL=1e-100
    fi
    # remove redundant markers
    Rscript scripts/remove_redundant_sig_top-non-sig_markers.R \
            ${GWAS} ${LDFILE} ${FOLDER} \
            --n-top-non-sig=${NONSIG} \
            --pval-threshold=${PVAL} \
            --ld-threshold=${LD}
  done
done
