#!/bin/bash
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
#SBATCH --mem=50gb
#SBATCH -J gcta_pve_all-markers
#SBATCH -o /home/hirschc1/burns756/empirical_sim/eo_files/%x_%j.out
#SBATCH -e /home/hirschc1/burns756/empirical_sim/eo_files/%x_%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=burns756@umn.edu
#SBATCH --no-requeue

# go to project folder
cd ~/empirical_sim/

echo "Creating binary ped (.bed) file"

# set plink file name
PLK=$(basename ${HMP} .hmp.txt)

# transform hmp to plk
run_pipeline.pl -Xmx50g -importGuess ${HMP} \
                -export ${FOLDER}/${PLK} \
                -exportType Plink

# Load plink
module load plink/1.90b6.10

# create bed file (binary ped)
plink --file ${FOLDER}/${PLK}.plk \
      --make-bed \
      --out ${FOLDER}/${PLK}

echo "Calculating genetic relationship matrices"

# calculate additive genetic relationship matrix
gcta --bfile ${FOLDER}/${PLK} \
     --make-grm \
     --thread-num 10 \
     --out ${FOLDER}/${PLK}

# calculate dominance genetic relationship matrix
gcta --bfile ${FOLDER}/${PLK} \
     --make-grm-d \
     --thread-num 10 \
     --out ${FOLDER}/${PLK}

# get an array with trait names -- each subfolder of ${FOLDER} should be a trait
traits=($(ls -d ${FOLDER}/*/))
for (( i = 0; i < ${#traits[@]}; i++ )); do
  traits[${i}]=$(basename ${traits[${i}]})
done

echo "Estimating variance explained by markers for different traits"

for trait in ${traits[@]}; do
  echo "  ${trait}"
  # set variables
  PHENO=${FOLDER}/${trait}/${trait}.phen
  OUT=${FOLDER}/${trait}/pve_all-markers.${trait}
  # estimate variance explained by additive markers
  gcta --grm ${FOLDER}/${PLK} \
       --pheno ${PHENO} \
       --reml \
       --out ${OUT}.add \
       --thread-num 10
  # estimate variance explained by dominance markers
  gcta --grm ${FOLDER}/${PLK}.d \
       --pheno ${PHENO} \
       --reml \
       --out ${OUT}.dom \
       --thread-num 10
done
