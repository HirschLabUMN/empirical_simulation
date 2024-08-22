#!/bin/bash
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=20gb
#SBATCH -J select_top_1000_non-redundant_markers
#SBATCH -o /home/hirschc1/burns756/empirical_sim/eo_files/%x_%j.out
#SBATCH -e /home/hirschc1/burns756/empirical_sim/eo_files/%x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=burns756@umn.edu
#SBATCH --no-requeue

echo $trait
# Create a temporary file that is the first two columns of the LD file - this will be the list of potential markers and the markers they are in LD with
cut -f1,2 analysis/gwas/ld_topall_0.9_cut_sig_non-sig-markers.ld > analysis/gwas/${trait}marker_ld_temp.txt
# Create a temporary file for storing the potential top markers in
cut -f1 analysis/gwas/${trait}_markers_ranked.txt >  analysis/gwas/${trait}potential_top_markers_temp.txt
# Create an empty file to store top markers in
echo -n > analysis/gwas/${trait}_top_markers.txt
# While loop to run until the file has 1000 lines or until the file of potential markers in empty
while [ $(wc -l analysis/gwas/${trait}_top_markers.txt | cut -f1 -d" ") -lt 1000 ] && [ $(wc -l analysis/gwas/${trait}marker_ld_temp.txt | cut -f1 -d" ") -gt 0 ];do
    echo Currently on line: $(wc -l analysis/gwas/${trait}_top_markers.txt | cut -f1 -d" ")
    # Get the top marker from the file
    top_marker=$(head -n2 analysis/gwas/${trait}potential_top_markers_temp.txt | tail -n1 | cut -f1)
    # Save this marker to the file of top markers
    echo $top_marker >> analysis/gwas/${trait}_top_markers.txt
    # Remove this top marker from the list of potential markers
    grep -v $top_marker analysis/gwas/${trait}potential_top_markers_temp.txt > analysis/gwas/${trait}greptmp.txt
    mv analysis/gwas/${trait}greptmp.txt analysis/gwas/${trait}potential_top_markers_temp.txt
    # Grep this top marker from the file of ld markers (both columns)
    grep $top_marker analysis/gwas/${trait}marker_ld_temp.txt | cut -f1 > analysis/gwas/${trait}markers_2_remove.txt
    grep $top_marker analysis/gwas/${trait}marker_ld_temp.txt | cut -f2 >> analysis/gwas/${trait}markers_2_remove.txt
    # Inverse grep (-v) the list of markers in LD to the list of potential markers to remove anything in LD
    grep -v -f analysis/gwas/${trait}markers_2_remove.txt analysis/gwas/${trait}potential_top_markers_temp.txt > analysis/gwas/${trait}greptmp.txt
    mv analysis/gwas/${trait}greptmp.txt analysis/gwas/${trait}potential_top_markers_temp.txt
    # Delete the markers to remove file since it will be different each time
    rm analysis/gwas/${trait}markers_2_remove.txt
    # End while loop
done

# Remove temp files
rm analysis/gwas/${trait}marker_ld_temp.txt
rm analysis/gwas/${trait}markers_2_remove.txt
rm analysis/gwas/${trait}potential_top_markers_temp.txt
