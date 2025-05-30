#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=10gb
#SBATCH -J counting_missing_data_per_site
#SBATCH -o /home/hirschc1/burns756/empirical_sim/eo_files/%x_%j.out
#SBATCH -e /home/hirschc1/burns756/empirical_sim/eo_files/%x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=burns756@umn.edu
#SBATCH --no-requeue

hapmap_file="data/usda_hybrids_projected-SVs-SNPs.poly.low-missing.hmp.txt"

# Initialize variables
line_count=0

# Create blank file
echo "" > 'data/missing_data_per_marker.txt'

# Process file line by line
while IFS= read -r line; do
  # Count occurrences of 'NN' in the line
  echo "$line" | grep -o 'NN' | wc -l >> 'data/missing_data_per_marker.txt'

  ((line_count++))
done < "$hapmap_file"
