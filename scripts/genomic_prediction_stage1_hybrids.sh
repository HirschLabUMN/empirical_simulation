#!/bin/bash

#!/bin/bash

function usage {
  echo
  echo "Usage: $(basename $0) [-t] [-h]" 2>&1
  echo
  echo '  -t    choose between "sim" or "real" datatype'
  echo '  -h    print help'
  echo
  exit 1
}


# set default options
DATATYPE=sim

# get command-line optional arguments
optstring=":t:h"

while getopts ${optstring} arg; do
  case ${arg} in
    t) DATATYPE=${OPTARG} ;;
    h) usage ;;

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
# for opt in DATATYPE; do
#   echo "  ${opt}=${!opt}"
# done


# go to project folder
cd /home/candy/rafa/genomic_prediction/hybrids

echo "job started @ $(date) ---"

for TRAIT in YLD EHT PHT Moisture; do

  # if sim data
  if [[ ${DATATYPE} == "sim" ]]; then

    for NMARKERS in 1 3 5 10 20 30 40 50 60 70 80 90 100 150 200 250 300 350; do
      for EFFECT in 1 0.1; do
        # get folder name
        FOLDER=analysis/sim_traits/${TRAIT}/traits/avg_rank_all/n_markers_${NMARKERS}/effects_${EFFECT}/
        # check if sim file exists first
        INFILE=$(echo ${FOLDER}/Simulated_Data*)
        if [[ -e ${INFILE} ]]; then
          echo "  ${FOLDER}"
          # get blues -- without adding environmental weights
          Rscript scripts/genomic_prediction_stage1_hybrids.R ${INFILE} ${FOLDER} > ${FOLDER}/genomic_prediction_stage1_hybrids.log 2>&1 &
          # get blues -- adding environmental weights
          Rscript scripts/genomic_prediction_stage1_hybrids.R ${INFILE} ${FOLDER} --envs-weight > ${FOLDER}/genomic_prediction_stage1_hybrids_weighted.log 2>&1 &
          wait
        fi
      done
    done

  fi

  # if real data
  if [[ ${DATATYPE} == "real" ]]; then
    FOLDER=analysis/empirical_traits/${TRAIT}
    INFILE=data/NIFA_CompleteDataset.${TRAIT}.txt
    echo "  ${FOLDER}"
    # get blues -- without adding environmental weights
    Rscript scripts/genomic_prediction_stage1_hybrids.R ${INFILE} ${FOLDER} > ${FOLDER}/genomic_prediction_stage1_hybrids.log 2>&1 &
    # get blues -- adding environmental weights
    Rscript scripts/genomic_prediction_stage1_hybrids.R ${INFILE} ${FOLDER} --envs-weight > ${FOLDER}/genomic_prediction_stage1_hybrids_weighted.log 2>&1 &
    wait
  fi

done

echo "job finished @ $(date)"
