#!/bin/bash
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
#SBATCH --mem=50gb
#SBATCH -J gcta_pve_sig-plus-non-sig-markers
#SBATCH -o /home/hirschc1/burns756/empirical_sim/eo_files/%x_%j.out
#SBATCH -e /home/hirschc1/burns756/empirical_sim/eo_files/%x_%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=burns756@umn.edu
#SBATCH --no-requeue

# go to project folder
cd ~/empirical_sim

# Load plink
module load plink/1.90b6.10

# get an array with trait names -- each subfolder of ${GCTAFOLDER} should be a trait
traits=($(ls -d ${GCTAFOLDER}/*/))
for (( i = 0; i < ${#traits[@]}; i++ )); do
  traits[${i}]=$(basename ${traits[${i}]})
done

# filter hmp and calculate pve per trait and genetic model
for trait in ${traits[@]}; do
  for model in additive dominant; do

    # get filename of list with significant markers
    LIST=${GWASFOLDER}/${trait}/${model}_model/permutation/list_top_non-redudant_gwas_markers.txt

    echo "${trait}: ${model} model"
    echo "  creating binary ped (.bed) file"

    # set plink file name
    PLK=$(basename ${HMP} .hmp.txt).top_non-redudant_markers.${model}-gwas

    # transform hmp to plk
    run_pipeline.pl -Xmx50g -importGuess ${HMP} \
                    -includeSiteNamesInFile ${LIST} \
                    -export ${GCTAFOLDER}/${trait}/${PLK} \
                    -exportType Plink 1>&2

    # create bed file (binary ped)
    plink --file ${GCTAFOLDER}/${trait}/${PLK}.plk \
          --make-bed \
          --out ${GCTAFOLDER}/${trait}/${PLK} 1>&2

    echo "  calculating genetic relationship matrices"

    if [[ ${model} == "additive" ]]; then
      # calculate additive genetic relationship matrix
      gcta --bfile ${GCTAFOLDER}/${trait}/${PLK} \
           --make-grm \
           --thread-num 10 \
           --out ${GCTAFOLDER}/${trait}/${PLK} 1>&2
    fi

    if [[ ${model} == "dominant" ]]; then
      # calculate dominance genetic relationship matrix
      gcta --bfile ${GCTAFOLDER}/${trait}/${PLK} \
           --make-grm-d \
           --thread-num 10 \
           --out ${GCTAFOLDER}/${trait}/${PLK} 1>&2
    fi

    echo "  estimating variance explained by markers for different traits"

    # define variables
    PHENO=${GCTAFOLDER}/${trait}/${trait}.phen
    OUT=${GCTAFOLDER}/${trait}/pve_top_non-redudant_markers.${trait}

    if [[ ${model} == "additive" ]]; then
      # estimate variance explained by additive markers
      gcta --grm ${GCTAFOLDER}/${trait}/${PLK} \
           --pheno ${PHENO} \
           --reml \
           --out ${OUT}.add \
           --thread-num 10 1>&2
    fi

    if [[ ${model} == "dominant" ]]; then
      # estimate variance explained by dominance markers
      gcta --grm ${GCTAFOLDER}/${trait}/${PLK}.d \
           --pheno ${PHENO} \
           --reml \
           --out ${OUT}.dom \
           --thread-num 10 1>&2
    fi

  done
done
