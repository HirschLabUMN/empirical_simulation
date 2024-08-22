#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
#SBATCH --mem=400gb
#SBATCH -J all_non-sig_ld_calc
#SBATCH -o /home/hirschc1/burns756/empirical_sim/eo_files/%x_%j.out
#SBATCH -e /home/hirschc1/burns756/empirical_sim/eo_files/%x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=burns756@umn.edu
#SBATCH --no-requeue

# Load plink
module load plink/1.90b6.10

# Set working directory
cd /home/hirschc1/burns756/empirical_sim

HMP=data/usda_hybrids_projected-SVs-SNPs.poly.low-missing.pruned-10.geno-miss-0.25.hmp.txt
# get folder name
FOLDER=analysis/gwas/
# get list of markers to calcukate LD
NONSIGMARKERS=${FOLDER}/list_sig_all_non-sig_markers.txt
# define plk filename
PLK=${FOLDER}/usda_hybrids_projected-SVs-SNPs.poly.low-missing.pruned-10.geno-miss-0.25.topall-non-sig-hits
# define ld output name
OUT=${FOLDER}/ld_topall_sig_non-sig-markers
# hmp2plk
run_pipeline.pl -Xmx40g -importGuess ${HMP} \
                        -includeSiteNamesInFile ${NONSIGMARKERS} \
                        -export ${PLK} \
                        -exportType Plink > /dev/null
# calculate LD within same chromosome
plink --file ${PLK}.plk \
        --make-founders \
        --r2 gz dprime with-freqs \
        --ld-window-r2 ${LD} \
        --ld-window 400000000 \
        --ld-window-kb 400000 \
        --out ${OUT}

# Print number of lines in ld file
echo $(wc -l ${OUT}.ld.gz)

# Reduce the file to the needed information - note that I kept the awk command for filtering for R2 > ${LD} just to be safe.
less ${FOLDER}/ld_topall_sig_non-sig-markers.ld.gz | awk -F ' +' '$10>${LD}' | awk -F ' +' '{print $4"\t"$8"\t"$10}' > ${FOLDER}/ld_topall_${LD}_cut_sig_non-sig-markers.ld

# Print the number of lines in the reduced file
echo $(wc -l ${FOLDER}/ld_topall_${LD}_cut_sig_non-sig-markers.ld)
