#!/bin/bash

function usage {
  echo
  echo "Usage: $(basename $0) [-d INT] [-e INT] [-f INT] [-h] [-i INT] [-p INT] [-P INT] [-s INT] [-t] [-w]" 2>&1
  echo
  echo '  -d    choose the starting iteration of marker dataset for prediction'
  echo '  -e    choose how to model add/dom effects'
  echo '  -f    number of folds for cross-validation'
  echo '  -g    choose how to model variance from random terms'
  echo '  -G    choose how to model variance from random terms for initial values'
  echo '  -h    print help'
  echo '  -i    number of cross-validation iterations'
  echo '  -o    name of output with commands to be run in parallel'
  echo '  -r    choose how to model variance from residual terms'
  echo '  -R    choose how to model variance from residual terms for initial values'
  echo '  -s    seed number'
  echo '  -t    choose between "sim" or "real" datatype'
  echo '  -w    add this flag to use weighted BLUEs from 1st stage'
  echo
  exit 1
}

function warning {
  echo
  echo "Argument for -$1 should be an integer" 2>&1
  echo
  exit 1
}

# set default options
NDATASET=1
MODELEFF=A
GSTR=idv
GSTRINIT=idv
NFOLDS=5
NITER=3
OUTFILE=commands.txt
RSTR=idv
RSTRINIT=idv
SEED=2021
DATATYPE=sim
WEIGHT=

# get command-line optional arguments
optstring=":d:e:f:g:G:hi:o:r:R:s:t:w"

while getopts ${optstring} arg; do
  case ${arg} in
    d) [[ $OPTARG =~ ^[0-9]+$ ]] && NDATASET=${OPTARG} || warning ${arg} ;;
    e) MODELEFF=${OPTARG} ;;
    f) [[ $OPTARG =~ ^[0-9]+$ ]] && NFOLDS=${OPTARG} || warning ${arg} ;;
    g) GSTR=${OPTARG} ;;
    G) GSTRINIT=${OPTARG} ;;
    h) usage ;;
    i) [[ $OPTARG =~ ^[0-9]+$ ]] && NITER=${OPTARG} || warning ${arg} ;;
    o) OUTFILE=${OPTARG} ;;
    r) RSTR=${OPTARG} ;;
    R) RSTRINIT=${OPTARG} ;;
    s) [[ $OPTARG =~ ^[0-9]+$ ]] && SEED=${OPTARG} || warning ${arg} ;;
    t) DATATYPE=${OPTARG} ;;
    w) WEIGHT=--envs-weight ;;

    :) echo "Error: -${OPTARG} requires an argument"; usage ;;
    *) echo "Unknown option: -${OPTARG}"; usage ;;

  esac
done

# shift argument index so positional arguments can start at $1
shift $(($OPTIND - 1))

# assert positional arguments
if [ $# -ge 1 ]; then
    echo -e "Ignoring positional arguments (none needed)\n"
fi

# # debug
# echo -e "positional args:\n  ${@}"
# echo "optional args:"
# for opt in NDATASET MODELEFF GSTR GSTRINIT NFOLDS NITER OUTFILE RSTR RSTRINIT SEED DATATYPE WEIGHT; do
#   echo "  ${opt}=${!opt}"
# done


# run predictions (stage 2)

# go to project folder
#cd ~/empirical_sim

# get genotypic data
GDATA=analysis/sim_traits/predictors/iter${NDATASET}/hybrids_geno.iter${NDATASET}.hmp.txt

for TRAIT in YLD EHT PHT Moisture; do

  # if sim data
  if [[ ${DATATYPE} == "sim" ]]; then

    for NMARKERS in 1 3 5 10 20 30 40 50 60 70 80 90 100 150 200 250 300 350; do
      for EFFECT in 1 0.1; do
        # get folder name
        FOLDER=analysis/sim_traits/${TRAIT}/traits/avg_rank_all/n_markers_${NMARKERS}/effects_${EFFECT}

        # check if blues file exists first
        INFILE=$(echo ${FOLDER}/blues_1st_stage.txt)
        if [[ -e ${INFILE} ]]; then

          # get pheno data and set output folder
          if [[ -z ${WEIGHT} ]]; then
            # if using blues without weights
            PDATA=${FOLDER}/blues_1st_stage.txt
            OUT=${FOLDER}/prediction_results/iter${NDATASET}/model-${MODELEFF}/Gstr-${GSTR}_Rstr-${RSTR}
            LOG=${OUT}/genomic_prediction_stage2_hybrids.log
            mkdir -p ${OUT}
          else
            # if using blues with weights
            PDATA=${FOLDER}/blues_1st_stage_weighted.txt
            OUT=${FOLDER}/prediction_results_weighted/iter${NDATASET}/model-${MODELEFF}/Gstr-${GSTR}_Rstr-${RSTR}
            LOG=${OUT}/genomic_prediction_stage2_hybrids_weighted.log
            mkdir -p ${OUT}
          fi

          # echo "    ${FOLDER}"
          # run prediction with cv1
          echo "Rscript scripts/genomic_prediction_stage2_hybrids.R ${GDATA} ${PDATA} ${OUT} --cv-type=CV1 --n-folds=${NFOLDS} --cv-iter=${NITER} --impute-effect=Add --impute-type=Middle --seed=${SEED} --model-effect=${MODELEFF} --model-Gstr=${GSTR} --model-Rstr=${RSTR} --model-Gstr-init=${GSTRINIT} --model-Rstr-init=${RSTRINIT} ${WEIGHT} > ${LOG%.log}.CV1.log" >> ${OUTFILE}
          # run prediction with cv2
          echo "Rscript scripts/genomic_prediction_stage2_hybrids.R ${GDATA} ${PDATA} ${OUT} --cv-type=CV2 --n-folds=${NFOLDS} --cv-iter=${NITER} --impute-effect=Add --impute-type=Middle --seed=${SEED} --model-effect=${MODELEFF} --model-Gstr=${GSTR} --model-Rstr=${RSTR} --model-Gstr-init=${GSTRINIT} --model-Rstr-init=${RSTRINIT} ${WEIGHT} > ${LOG%.log}.CV2.log" >> ${OUTFILE}


        fi

      done
    done

  fi

  # if real data
  if [[ ${DATATYPE} == "real" ]]; then
    # get folder name
    FOLDER=analysis/empirical_traits/${TRAIT}
    # get pheno data and set output folder
    if [[ -z ${WEIGHT} ]]; then
      # if using blues without weights
      PDATA=${FOLDER}/blues_1st_stage.txt
      OUT=${FOLDER}/prediction_results/iter${NDATASET}/model-${MODELEFF}/Gstr-${GSTR}_Rstr-${RSTR}
      LOG=${OUT}/genomic_prediction_stage2_hybrids.log
      mkdir -p ${OUT}
    else
      # if using blues with weights
      PDATA=${FOLDER}/blues_1st_stage_weighted.txt
      OUT=${FOLDER}/prediction_results_weighted/iter${NDATASET}/model-${MODELEFF}/Gstr-${GSTR}_Rstr-${RSTR}
      LOG=${OUT}/genomic_prediction_stage2_hybrids_weighted.log
      mkdir -p ${OUT}
    fi

    # echo "    ${FOLDER}"
    # run prediction with cv1
    echo "Rscript scripts/genomic_prediction_stage2_hybrids.R ${GDATA} ${PDATA} ${OUT} --cv-type=CV1 --n-folds=${NFOLDS} --cv-iter=${NITER} --impute-effect=Add --impute-type=Middle --seed=${SEED} --model-effect=${MODELEFF} --model-Gstr=${GSTR} --model-Rstr=${RSTR} --model-Gstr-init=${GSTRINIT} --model-Rstr-init=${RSTRINIT} ${WEIGHT} > ${LOG%.log}.CV1.log" >> ${OUTFILE}
    # run prediction with cv2
    echo "Rscript scripts/genomic_prediction_stage2_hybrids.R ${GDATA} ${PDATA} ${OUT} --cv-type=CV2 --n-folds=${NFOLDS} --cv-iter=${NITER} --impute-effect=Add --impute-type=Middle --seed=${SEED} --model-effect=${MODELEFF} --model-Gstr=${GSTR} --model-Rstr=${RSTR} --model-Gstr-init=${GSTRINIT} --model-Rstr-init=${RSTRINIT} ${WEIGHT} > ${LOG%.log}.CV2.log" >> ${OUTFILE}

  fi

done
