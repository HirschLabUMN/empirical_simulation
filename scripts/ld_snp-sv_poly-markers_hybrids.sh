#!/bin/bash
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=50gb
#SBATCH -J ld_snp-sv_poly-markers_${CHR}
#SBATCH -o /home/hirschc1/della028/projects/genomic_prediction/hybrids/analysis/MSI_dump/%x_%j.out
#SBATCH -e /home/hirschc1/della028/projects/genomic_prediction/hybrids/analysis/MSI_dump/%x_%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=della028@umn.edu
#SBATCH --no-requeue

# go to project folder
cd ~/projects/genomic_prediction/hybrids

# transform hmp into plink format
run_pipeline.pl -Xmx100g -importGuess data/usda_hybrids_projected-SVs-SNPs.poly.low-missing.chr${CHR}.hmp.txt -export ${SCRATCH}/usda_hybrids_projected-SVs-SNPs.poly.low-missing.chr${CHR} -exportType Plink

# go to scratch folder
cd ${SCRATCH}

# calculate LD
plink --file usda_hybrids_projected-SVs-SNPs.poly.low-missing.chr${CHR}.plk --make-founders --r2 gz dprime with-freqs --ld-window-r2 0 --ld-window 1000 --ld-window-kb 1 --geno 0.25 --out usda_hybrids_projected-SVs-SNPs.poly.low-missing.chr${CHR}

# keep only snp and sv r2 (excluding translocations)
zcat usda_hybrids_projected-SVs-SNPs.poly.low-missing.chr${CHR}.ld.gz | head -n 1 > ~/projects/genomic_prediction/hybrids/analysis/ld/ld_usda_hybrids_poly-snp-sv_only.chr${CHR}.ld
zcat usda_hybrids_projected-SVs-SNPs.poly.low-missing.chr${CHR}.ld.gz | awk '$3 ~ /^del|^dup|^ins|^inv|^tra/ && $7 !~ /^del|^dup|^ins|^inv|^tra/ || $3 !~ /^del|^dup|^ins|^inv|^tra/ && $7 ~ /^del|^dup|^ins|^inv|^tra/' - >> ~/projects/genomic_prediction/hybrids/analysis/ld/ld_usda_hybrids_poly-snp-sv_only.chr${CHR}.ld
