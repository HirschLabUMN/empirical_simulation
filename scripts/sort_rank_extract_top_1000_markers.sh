# Create a unix for loop that reads in a file with the naming structure of gwas/*/*_model/permutation/GAPIT.MLM.BLUE.GWAS.Results.csv, sorts it numerically ascending by column 4, and then extracts the first column and adds the row number to each line. This should then be saved as a new file.
for trait in EHT; do
    for model in additive dominant; do
        for env in COR19 MIN19 MIN20 URB19; do
            echo $trait $env $model
            echo ${trait} ${env}_${model} > ${trait}_${env}_${model}_sorted.txt
            tail -n+2 analysis/gwas/${trait}_${env}/${model}_model/permutation/GAPIT.MLM.BLUE.GWAS.Results.csv | sort -n -k4 -t',' | cut -f1 -d',' | awk '{print $1, NR}' >> ${trait}_${env}_${model}_sorted.txt
        done
    done
done

sort -k1 -t' ' EHT_COR19_additive_sorted.txt

# Sort files based on the first column and then save over the file
for trait in EHT; do
    for model in additive dominant; do
        for env in COR19 MIN19 MIN20 URB19; do
            echo $trait $env $model
            sort -k1 -t' ' ${trait}_${env}_${model}_sorted.txt > ${trait}_${env}_${model}_resorted.txt
        done
    done
    # Merge the files together
    paste ${trait}_*_*_resorted.txt > ${trait}_merged.txt
done