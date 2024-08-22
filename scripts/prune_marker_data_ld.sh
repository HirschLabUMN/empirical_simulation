#!/bin/bash
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=40gb
#SBATCH -J prune_marker_data_ld
#SBATCH -o /home/hirschc1/burns756/empirical_sim/eo_files/%x_%j.out
#SBATCH -e /home/hirschc1/burns756/empirical_sim/eo_files/%x_%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=burns756@umn.edu
#SBATCH --no-requeue

# go to project folder
cd ~/empirical_sim/

# set intermediate filenames
PLK=$(basename ${HMP} .hmp.txt)
PRUNNED=$(echo ${PLK}.marker-ids)

# transform hmp to plk
run_pipeline.pl -Xmx40g -importGuess ${HMP} \
                -export ${FOLDER}/${PLK} \
                -exportType Plink

# Load plink
module load plink/1.90b6.10
# prune plink file by ld
plink --file ${FOLDER}/${PLK}.plk \
      --indep-pairwise ${WINSIZE} 'kb' ${VARCOUNT} ${R2} \
      --geno ${GENOMISS} \
      --out ${FOLDER}/${PRUNNED} \
      --allow-extra-chr \
      --make-founders

# # filter file to have only pruned variants
# plink --file NAM_rils_projected-reseq-SNPs-only.all-RILs.chr-${CHR}.v8.plk --extract ${FOLDER}/${PRUNNED}.prune.in --make-bed --freq --out NAM_rils_projected-reseq-SNPs-only.all-RILs.chr-${CHR}.v8.prunned.plk --allow-extra-chr --make-founders

# filter hapmap
run_pipeline.pl -Xmx40g -importGuess ${HMP} \
                -includeSiteNamesInFile ${FOLDER}/${PRUNNED}.prune.in \
                -export ${OUT} \
                -exportType HapmapDiploid
