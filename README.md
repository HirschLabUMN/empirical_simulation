# Optimizing trait simulations with GWAS

by Michael Burns, Rafael Della Coletta, Samuel Fernandes, Martin Bohn, Alex Lipka, and Candice Hirsch

> This project aims to determine the genetic architecture of traits by performing GWAS on empirical data obtained from maize hybrids, simulating traits mirroring these architectures, and running prediction models on both simulated and empirical data. For more details about this analysis, please read the manuscript: **ADD DOI HERE**


<!-- TOC START min:1 max:3 link:true asterisk:false update:true -->
- [Optimizing trait simulations with GWAS](#optimizing-trait-simulations-with-gwas)
  - [Important notes!](#important-notes)
  - [Genotypic data](#genotypic-data)
    - [QC](#qc)
  - [GWAS](#gwas)
    - [LD decay](#ld-decay)
    - [Create genotypic datasets](#create-genotypic-datasets)
    - [Calculate BLUEs](#calculate-blues)
    - [Run association analysis](#run-association-analysis)
  - [Simulations](#simulations)
    - [Correlation among environments](#correlation-among-environments)
    - [Phenotypic simulations](#phenotypic-simulations)
    - [QC simulations](#qc-simulations)
    - [Comparison to real data](#comparison-to-real-data)
    - [Correlation curve](#correlation-curve)
  - [Genomic prediction](#genomic-prediction)
    - [Create datasets](#create-datasets)
    - [Predict simulated traits](#predict-simulated-traits)
    - [Predict real traits](#predict-real-traits)
<!-- TOC END -->



## Important notes!

All the data necessary to replicate this analysis are available either at Data Repository for the University of Minnesota (DRUM) (https://hdl.handle.net/11299/252793) or in the supplemental materials in the publication. When submitting the files for publication, I changed their names for easier reference within the manuscript. However, please rename them to their original filenames before running any of the scripts according to this table:

| Source                                      | File | Name                                      |
| ------------------------------------------- | ---- | ----------------------------------------- |
| [DRUM](https://hdl.handle.net/11299/252793) | S2   | usda_rils_projected-SVs-SNPs.poly.hmp.txt |
| Manuscript supplemental material            | S3   | NIFA_CompleteDataset.xlsx                 |

All scripts necessary for data analysis are available at the `scripts` folder in this GitHub repository, and I recommend that you create a folder called `data` to keep the files above and another called `analysis` to keep any other files generated throughout this analysis.

All the analysis, unless indicated otherwise, were run using the Minnesota Supercomputing Institute (MSI) clusters, which uses the SLURM scheduler. You may need to adapt some scripts if you use a different scheduler.

Finally, here is a list of software necessary to run the analysis and their respective versions:

| Software | Version |
| -------- | ------- |
| R        | 3.6.0   |
| Python   | 3.6.6   |
| GNU bash | 4.2.46  |
| Tassel   | 5.2.56  |
| Plink    | 1.9     |




## Genotypic data

The genotypic data of the RILs `data/usda_rils_projected-SVs-SNPs.poly.hmp.txt` were obtained from the previous project for simulating inbred traits. I will then use information available on `data/NIFA_CompleteDataset.xlsx` to get pedigree of hybrids planted in 2020 USDA trials, and then create hybrid genotypes based on genotypes of parental (RIL) data.

```bash
# create file with svs ids
cut -f 1 data/usda_rils_projected-SVs-SNPs.poly.hmp.txt | grep -P "^del|^dup|^inv|^ins|^tra" > data/SVs_IDs_poly.txt

# split genotypic dataset into chr to speed things up
for chr in {1..10}; do
  head -n 1 data/usda_rils_projected-SVs-SNPs.poly.hmp.txt > data/usda_rils_projected-SVs-SNPs.chr${chr}.poly.hmp.txt
  awk -v chr="$chr" '$3 == chr' data/usda_rils_projected-SVs-SNPs.poly.hmp.txt >> data/usda_rils_projected-SVs-SNPs.chr${chr}.poly.hmp.txt
done

# create hybrid genotypes chr by chr
cd /home/hirschc1/della028/projects/genomic_prediction/hybrids/
for chr in {1..10}; do
  # define arguments for shell script
  HYBINFO=data/NIFA_CompleteDataset.xlsx
  RILGENO=data/usda_rils_projected-SVs-SNPs.chr${chr}.poly.hmp.txt
  OUT=data/usda_hybrids_projected-SVs-SNPs.chr${chr}.poly.hmp.txt
  HET=missing
  # create hybrids
  sbatch --export=HYBINFO=${HYBINFO},RILGENO=${RILGENO},OUT=${OUT},HET=${HET} scripts/create_hybrid_genotypes.sh
done

# remove monomorphic markers (maf < 0.05)
for chr in {1..10}; do
  run_pipeline.pl -Xmx40g \
                  -importGuess data/usda_hybrids_projected-SVs-SNPs.chr${chr}.poly.hmp.txt \
                  -FilterSiteBuilderPlugin \
                  -siteMinAlleleFreq 0.05 \
                  -endPlugin \
                  -export data/usda_hybrids_projected-SVs-SNPs.chr${chr}.poly.maf-filter,data/usda_hybrids_projected-SVs-SNPs.chr${chr}.poly.maf-filter2 -exportType HapmapDiploid
  # remove unnecessary file
  rm data/usda_hybrids_projected-SVs-SNPs.chr${chr}.poly.maf-filter2.json.gz
done

# merge chromosomes
cat data/usda_hybrids_projected-SVs-SNPs.chr1.poly.maf-filter.hmp.txt > data/usda_hybrids_projected-SVs-SNPs.poly.maf-filter.hmp.txt
for chr in {2..10}; do
  sed 1d data/usda_hybrids_projected-SVs-SNPs.chr${chr}.poly.maf-filter.hmp.txt >> data/usda_hybrids_projected-SVs-SNPs.poly.maf-filter.hmp.txt
done

# compress intermediate input files to save space
gzip data/usda_rils_projected-SVs-SNPs.*
gzip data/usda_hybrids_projected-SVs-SNPs.chr*
```



### QC

Prior running GWAS or genomic prediction, any missing data will be imputed by GAPIT. Thus, I just wanted to have an idea about how much having missing data the genotypic dataset currently has, and also plot marker distribution along the chromosomes.

```bash
# qc
sbatch --export=HMP=data/usda_hybrids_projected-SVs-SNPs.poly.maf-filter.hmp.txt,SVS=data/SVs_IDs_poly.txt,FOLDER=analysis/qc/snp-sv_hmp scripts/qc_snp-sv_hmp.sh
```

| marker type | missing filter | marker #  |
| ----------- | -------------- | --------- |
| SVs         | no filter      | 8,416     |
|             | 0.2            | 3,953     |
|             | 0.3            | 6,159     |
|             | 0.4            | 7,217     |
|             | 0.5            | 7,880     |
| SNPs        | no filter      | 8,510,450 |
|             | 0.2            | 3,939,659 |
|             | 0.3            | 4,381,741 |
|             | 0.4            | 4,508,709 |
|             | 0.5            | 4,818,301 |

Based on the information above and the distribution of markers along chromosomes before and after filtering (see plots in `analysis/qc/snp-sv_hmp`), we chose to keep markers with less than 50% missing data.

```bash
# create file with markers to keep
cat analysis/qc/snp-sv_hmp/SVs_to_keep.missing-threshold-0.5.txt analysis/qc/snp-sv_hmp/SNPs_to_keep.missing-threshold-0.5.txt > analysis/qc/snp-sv_hmp/markers_to_keep.missing-threshold-0.5.txt
# filter genotypic dataset
sbatch --export=IN=data/usda_hybrids_projected-SVs-SNPs.poly.maf-filter.hmp.txt,LIST=analysis/qc/snp-sv_hmp/markers_to_keep.missing-threshold-0.5.txt,OUT=data/usda_hybrids_projected-SVs-SNPs.poly.low-missing.hmp.txt scripts/filter_hmp_based_on_marker_names.sh
# compress previous hmp
gzip data/usda_hybrids_projected-SVs-SNPs.poly.maf-filter.hmp.txt

for chr in {1..10}; do
  head -n 1 data/usda_hybrids_projected-SVs-SNPs.poly.low-missing.hmp.txt > data/usda_hybrids_projected-SVs-SNPs.poly.low-missing.chr${chr}.hmp.txt
  awk -v chr="$chr" '$3 == chr' data/usda_hybrids_projected-SVs-SNPs.poly.low-missing.hmp.txt >> data/usda_hybrids_projected-SVs-SNPs.poly.low-missing.chr${chr}.hmp.txt
done
```




## GWAS

### LD decay

Visualizing LD decay is be important to define which LD window should be used in prediction, i.e. to find SNPs that are close enough to SVs to avoid LD due to population structure but also far enough so there's enough recombination happening between SNPs and SVs. From my simulation study with the inbred parents of hybrids, I found that LD decays very rapidly (~1kb) so I just wanted to check if the pattern holds for hybrids as well.

```bash
# create directory to store results
mkdir -p analysis/ld

# create directory to store tmp files
SCRATCH=/scratch.global/della028/hirsch_lab/genomic_prediction/hybrids/ld
mkdir -p ${SCRATCH}

# convert hapmap to plink format chr by chr using only markers from SNP chip
for chr in {1..10}; do
  run_pipeline.pl -Xmx10g -importGuess data/usda_hybrids_projected-SVs-SNPs.poly.low-missing.chr${chr}.hmp.txt \
                  -includeSiteNamesInFile <(sed 1d data/usda_22kSNPs_rils.hmp.txt | cut -f 1) \
                  -export ${SCRATCH}/usda_hybrids_SNP-chip-markers.chr${chr} \
                  -exportType Plink
done

# calculate LD
module load plink/1.90b6.10
WINDOW=1000
FILTER=0.25
for chr in {1..10}; do
  plink --file ${SCRATCH}/usda_hybrids_SNP-chip-markers.chr${chr}.plk \
        --make-founders --r2 gz dprime with-freqs --ld-window-r2 0 \
        --ld-window ${WINDOW}000 --ld-window-kb ${WINDOW} --geno ${FILTER} \
        --out ${SCRATCH}/usda_hybrids_SNP-chip-markers.chr${chr}.window-${WINDOW}kb.filter-${FILTER}
done

# merge ld files from each chromosome
zcat ${SCRATCH}/usda_hybrids_SNP-chip-markers.chr1.window-${WINDOW}kb.filter-${FILTER}.ld.gz > ${SCRATCH}/usda_hybrids_SNP-chip-markers.window-${WINDOW}kb.filter-${FILTER}.ld
for chr in {2..10}; do
  zcat ${SCRATCH}/usda_hybrids_SNP-chip-markers.chr${chr}.window-${WINDOW}kb.filter-${FILTER}.ld.gz | sed 1d >> ${SCRATCH}/usda_hybrids_SNP-chip-markers.window-${WINDOW}kb.filter-${FILTER}.ld
done

# plot decay
module load R/3.6.0
Rscript scripts/plot_ld_decay.R ${SCRATCH}/usda_hybrids_SNP-chip-markers.window-${WINDOW}kb.filter-${FILTER}.ld \
                                analysis/ld/ld_decay_hybrids_snps-chip-only_1000kb \
                                --unequal-windows
```

> LD decays very rapidly in hybrids too (as expected) and has very similar pattern to the LD decay of their inbred parents



### Create genotypic datasets

Prior GWAS, we need to prune the full SNP and SV dataset (`data/usda_hybrids_projected-SVs-SNPs.poly.low-missing.hmp.txt`) based on LD to decrease the number of perfectly correlated markers. I tested three different window sizes and three different missing data thresholds  to calculate LD to see how many markers would be removed.

```bash
# split by chromosome first to parallelize
for chr in {1..10}; do
  echo ${chr}
  head -n 1 data/usda_hybrids_projected-SVs-SNPs.poly.low-missing.hmp.txt > data/usda_hybrids_projected-SVs-SNPs.poly.low-missing.chr${chr}.hmp.txt
  awk -v chr="${chr}" '$3 == chr' data/usda_hybrids_projected-SVs-SNPs.poly.low-missing.hmp.txt >> data/usda_hybrids_projected-SVs-SNPs.poly.low-missing.chr${chr}.hmp.txt
done

VARCOUNT=1
R2=0.9
#for genomiss in 0.1 0.25 0.5; do
  for winsize in 10 100 1000; do
    genomiss=0.5 # There is a parallelization issue where you need to run each genomiss separately
    # folder to save intermediate files
    FOLDER=analysis/prune_markers_ld/winsize_${winsize}
    mkdir -p ${FOLDER}
    for chr in {1..10}; do
      # hapmap file to filter
      HMP=data/usda_hybrids_projected-SVs-SNPs.poly.low-missing.chr${chr}.hmp.txt
      # output name
      OUT=${FOLDER}/usda_hybrids_projected-SVs-SNPs.poly.low-missing.pruned-${winsize}.geno-miss-${genomiss}.chr${chr}.hmp.txt
      # run script
      sbatch --export=HMP=${HMP},FOLDER=${FOLDER},OUT=${OUT},WINSIZE=${winsize},VARCOUNT=${VARCOUNT},R2=${R2},GENOMISS=${genomiss} scripts/prune_marker_data_ld.sh
    done
  done
#done

# count how many markers were left after prunning
for genomiss in 0.1 0.25 0.5; do
  for winsize in 10 100 1000; do
    # svs
    echo "miss: ${genomiss} - win: ${winsize} - SVs"
    grep -P "^del|^dup|^ins|^inv|^tra" analysis/prune_markers_ld/winsize_${winsize}/usda_hybrids_projected-SVs-SNPs.poly.low-missing.pruned-${winsize}.geno-miss-${genomiss}.chr*.hmp.txt | wc -l
    # snps
    echo "miss: ${genomiss} - win: ${winsize} - SNPs"
    grep -v -P "^del|^dup|^ins|^inv|^tra" analysis/prune_markers_ld/winsize_${winsize}/usda_hybrids_projected-SVs-SNPs.poly.low-missing.pruned-${winsize}.geno-miss-${genomiss}.chr*.hmp.txt | wc -l
    echo ""
  done
done

# create new hapmap file with pruned markers
for genomiss in 0.1 0.25 0.5; do
  for winsize in 10 100 1000; do
    echo "${winsize} - ${genomiss}"
    cp analysis/prune_markers_ld/winsize_${winsize}/usda_hybrids_projected-SVs-SNPs.poly.low-missing.pruned-${winsize}.geno-miss-${genomiss}.chr1.hmp.txt data/usda_hybrids_projected-SVs-SNPs.poly.low-missing.pruned-${winsize}.geno-miss-${genomiss}.hmp.txt
    for chr in {2..10}; do
      sed 1d analysis/prune_markers_ld/winsize_${winsize}/usda_hybrids_projected-SVs-SNPs.poly.low-missing.pruned-${winsize}.geno-miss-${genomiss}.chr${chr}.hmp.txt >> data/usda_hybrids_projected-SVs-SNPs.poly.low-missing.pruned-${winsize}.geno-miss-${genomiss}.hmp.txt
    done
  done
done

# compress previous datasets
gzip data/usda_hybrids_projected-SVs-SNPs.poly.low-missing.chr*.hmp.txt
```

Number of markers remaining after pruning the genotypic dataset with different window sizes:

| missing data filter | window size | SVs remaining | SNPs remaining |
| ------------------- | ----------- | ------------- | -------------- |
| 0.1                 | 10 kb       | 81            | 208,916        |
| 0.1                 | 100 kb      | 61            | 57,527         |
| 0.1                 | 1000 kb     | 34            | 10,518         |
| 0.25                | 10 kb       | 3,110         | 359,990        |
| 0.25                | 100 kb      | 2,300         | 103,295        |
| 0.25                | 1000 kb     | 1,067         | 18,606         |
| 0.5                 | 10 kb       | 4,641         | 395,749        |
| 0.5                 | 100 kb      | 3,084         | 108,990        |
| 0.5                 | 1000 kb     | 1,267         | 17,920         |



### Calculate BLUEs

The next step is to prepare the phenotypic data. I wrote `scripts/get_BLUEs_from_empirical-data.R` to calculate the BLUEs for each environment and for each trait, as well their respective heritability, from the file `data/NIFA_CompleteDataset.xlsx`.

```bash
# Calculate BLUEs
# module load R/3.6.0 # If run on an HPC
# This will require asremlr
# The following R scripts were run locally
# EHT
Rscript scripts/get_BLUEs_from_empirical-data.R \
        data/NIFA_CompleteDataset.xlsx \
        EHT \
        data \
        analysis/pheno_BLUES_qc \
        --envs=COR19,MIN19,MIN20,URB19 \
        --checks=UIUC-CHECK1,UIUC-CHECK2,UIUC-CHECK3,UIUC-CHECK4,6049V2P

# PHT
Rscript scripts/get_BLUEs_from_empirical-data.R \
        data/NIFA_CompleteDataset.xlsx \
        PHT \
        data \
        analysis/pheno_BLUES_qc \
        --envs=COR19,COR20,MIN19,MIN20,URB19 \
        --checks=UIUC-CHECK1,UIUC-CHECK2,UIUC-CHECK3,UIUC-CHECK4,6049V2P

# Moisture
Rscript scripts/get_BLUEs_from_empirical-data.R \
        data/NIFA_CompleteDataset.xlsx \
        Moisture \
        data \
        analysis/pheno_BLUES_qc \
        --envs=BAY19,BEC-BL19,BEC-BL20,BEC-EP20,COR19,COR20,MIN19,MIN20,SYN19,SYN20,URB19 \
        --checks=UIUC-CHECK1,UIUC-CHECK2,UIUC-CHECK3,UIUC-CHECK4,6049V2P

# YLD
Rscript scripts/get_BLUEs_from_empirical-data.R \
        data/NIFA_CompleteDataset.xlsx \
        YLD \
        data \
        analysis/pheno_BLUES_qc \
        --envs=BEC-BL19,BEC-BL20,BEC-EP20,COR19,COR20,MIN19,MIN20,SYN19,SYN20,URB19 \
        --checks=UIUC-CHECK1,UIUC-CHECK2,UIUC-CHECK3,UIUC-CHECK4,6049V2P

# This part was run on MSI
# merge heritability files
Rscript scripts/merge_h2_files.R data 1stStage_heritabilities
```



### Run association analysis

We are going to use genotypic dataset markers pruned with window size of 10kb due to very fast LD decay of inbred parents and to maximize the number of both SNPs and SVs. We also chose the dataset filtered to have at most 25% of missing data to have more accurate LD calculations.

To run GWAS, I'm going to follow this three steps:

1. Numericalize genotypic data (use 0-1-2 for additive, or 0-1 for dominance)
2. Run PCA to determine how many PCs to use in GWAS
3. Create kinship matrix (VanRaden method)
3. Run GWAS with GAPIT (Q+K model)

```bash
# Make sure the following are installed
#r <- getOption("repos");
#r["CRAN"] <- "http://cran.rstudio.com/";
#options(repos=r);
#install.packages(c('BiocManager', "gplots", "ape", "EMMREML", "scatterplot3d", "bigmemory", "biganalytics"));


#if (!require("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install(c('multtest', 'genetics', 'snpStats', 'compiler'))

module load R/3.6.0

# step 1 - I had to download the gapit source code (found in `supp_materials` directory above) to my MSI directory to get this to work
for model in additive dominant; do
  # get input file
  HMP=data/usda_hybrids_projected-SVs-SNPs.poly.low-missing.pruned-10.geno-miss-0.25.hmp.txt
  # submit job
  sbatch --export=HMP=${HMP},MODEL=${model} scripts/hmp2numeric.sh # Uses R/4.0.4 to comply with GAPIT requirements
done

# step 2
for model in additive dominant; do
  # get input file
  HMP=data/usda_hybrids_projected-SVs-SNPs.poly.low-missing.pruned-10.geno-miss-0.25.${model}.num.txt
  # create folder
  FOLDER=analysis/pca_gwas/pruning_10/${model}_model
  mkdir -p ${FOLDER}
  # submit job
  sbatch --export=HMP=${HMP},FOLDER=${FOLDER} scripts/pca_prior_gwas.sh # Again, used R/4.0.4 to comply with GAPIT requirements
done

# step 3
# create kinship matrix
for model in additive dominant; do
  # get input file
  HMP=data/usda_hybrids_projected-SVs-SNPs.poly.low-missing.pruned-10.geno-miss-0.25.${model}.num.txt
  # create folder
  FOLDER=analysis/kinship/pruning_10/${model}_model
  mkdir -p ${FOLDER}
  # submit job
  sbatch --export=HMP=${HMP},FOLDER=${FOLDER} scripts/kinship_VanRaden.sh # Again, used R/4.0.4 to comply with GAPIT requirements
done

# step 4
# Be sure to install the following packages:
# r <- getOption("repos");
# r["CRAN"] <- "http://cran.rstudio.com/";
# options(repos=r);
# install.packages(c('doRNG', 'ggrepel'));

# based on PCAs above, I think 7 PCs is a good number to use in GWAS
PCS=7
GWASmodel=MLM
# run gwas (Q+K model) for each trait and environment
for pheno in data/1stStage_BLUEs.*_*.txt; do
  for model in additive dominant; do
    for pval_method in permutation; do
      # get input files
      HMP=data/usda_hybrids_projected-SVs-SNPs.poly.low-missing.pruned-10.geno-miss-0.25.${model}.num.txt
      MAP=data/usda_hybrids_projected-SVs-SNPs.poly.low-missing.pruned-10.geno-miss-0.25.${model}.map.txt
      KIN=analysis/kinship/pruning_10/${model}_model/GAPIT.Kin.VanRaden.csv
      PCA=analysis/pca_gwas/pruning_10/${model}_model/GAPIT.PCA.csv
      # create folder
      FOLDER=analysis/gwas/$(basename -a ${pheno} | cut -f 2 -d ".")/${model}_model/${pval_method}
      mkdir -p ${FOLDER}
      # submit job
      sbatch --export=PHENO=${pheno},HMP=${HMP},MAP=${MAP},KIN=${KIN},PCA=${PCA},FOLDER=${FOLDER},PCS=${PCS},GWASmodel=${GWASmodel},PVAL=${pval_method} scripts/gwas.sh # Again, used R/4.0.4 to comply with GAPIT requirements
    done
  done
done

# summarize GWAS results
# B73_RefGen_V4_chrm_info.txt is from https://ftp.ncbi.nlm.nih.gov/genomes/genbank/plant/Zea_mays/latest_assembly_versions/GCA_000005005.6_B73_RefGen_v4/GCA_000005005.6_B73_RefGen_v4_assembly_stats.txt
# centromeres_Schneider-2016-pnas_v4.bed is from https://www.biorxiv.org/content/10.1101/2023.02.28.530521v1
module load R/3.6.0
Rscript scripts/summarize_gwas_results.R \
        analysis/gwas \
        data/B73_RefGen_V4_chrm_info.txt \ 
        data/centromeres_Schneider-2016-pnas_v4.bed
```

> Note: In previous tests, when running FarmCPU model, we found that the marker `snp.1.306948638` was a significant hit in five GWAS, all of them associated with yield (`YLD_BEC-BL19`, `YLD_COR20`, `YLD_SYN19`, `YLD_SYN20`, and `YLD_URB19`) and had very high effect. This marker is located at the very end of chromosome 1 at position 306,948,638 on B73v4 reference (308,359,862 on B73v5). This marker didn't show up again in the Q+K model, and I won't have time to investigate this further.



#### Remove redundant significant hits

The Q+K model will generate peaks where more than one marker may be in high LD with the causative variant and with other markers. Thus, I need to calculate LD among significant markers (even across chromosomes) to see how much of the signal from the Q+K model is highly correlated due to LD.

```bash
# create an array with trait names
traits=()
for trait in $(ls -d analysis/gwas/*/); do
  #echo $trait
  traits+=$(basename $trait)
  traits+=" "
done

module load plink/1.90b6.10

# Create a list of significant GWAS results
for trait in ${traits[@]}; do
  echo "--- ${trait} ---"
  for model in additive dominant; do
    FOLDER=analysis/gwas/${trait}/${model}_model/permutation
    if [[ -e ${FOLDER}/MLM.BLUE.GWAS.Results.sig-hits.csv ]]; then
        # get folder name
        cut -f1 -d',' ${FOLDER}/MLM.BLUE.GWAS.Results.sig-hits.csv > ${FOLDER}/list_sig_markers.txt
        # get list of markers to calculate LD
        #SIGMARKERS=${FOLDER}/list_sig_markers.txt
    fi
  done
  echo ""
done

# calculate LD among sig markers
HMP=data/usda_hybrids_projected-SVs-SNPs.poly.low-missing.pruned-10.geno-miss-0.25.hmp.txt
for trait in ${traits[@]}; do
  echo "--- ${trait} ---"
  for model in additive dominant; do
    # get folder name
    FOLDER=analysis/gwas/${trait}/${model}_model/permutation
    # get list of markers to calculate LD
    SIGMARKERS=${FOLDER}/list_sig_markers.txt
    if [[ -e ${SIGMARKERS} ]]; then
      # define plk filename
      PLK=${FOLDER}/usda_hybrids_projected-SVs-SNPs.poly.low-missing.pruned-10.geno-miss-0.25.sig-hits
      # define ld output name
      OUT=${FOLDER}/ld_sig-markers
      # hmp2plk
      echo "Converting HMP to PLK"
      run_pipeline.pl -Xmx40g -importGuess ${HMP} \
                              -includeSiteNamesInFile ${SIGMARKERS} \
                              -export ${PLK} \
                              -exportType Plink > /dev/null
      # calculate LD across chromosomes
      echo "Calculating LD"
      plink --file ${PLK}.plk \
            --make-founders \
            --r2 dprime with-freqs inter-chr \
            --ld-window-r2 0 \
            --out ${OUT}
    fi
  done
  echo ""
done

# plot distribution of ld of sig markers
module load R/3.6.0
Rscript scripts/plot_ld_sig-markers.R analysis/gwas

# Create a list of the all significant and non-significant GWAS results 
cut -f1 data/usda_hybrids_projected-SVs-SNPs.poly.low-missing.pruned-10.geno-miss-0.25.hmp.txt | tail -n+2 > analysis/gwas/list_sig_all_non-sig_markers.txt
```
> Most of significant markers are in very high LD to each other.

After calculating LD, I need to (1) remove redundant significant markers (keep the one with highest effect), (2) remove non-significant markers in high LD to significant marker, and (3) remove non-redundant non-significant markers (keep the one with highest effect). The non-redundant set of markers will be used to simulate traits.

```bash
# Calculate LD among all markers and save only the SNP_A, SNP_B, and R2 columns, and only save the lines where R2 > 0.9
# State the hapmap file name
# determine LD threshold
LD=0.9

sbatch --export=${LD} scripts/all_non-sig_ld_calc.sh

# Determine the top markers to use from the GWAS
# Rank the markers across environments and models
module load R/3.6.0
for trait in EHT PHT Moisture YLD; do
  echo "---$trait---"
  Rscript scripts/rank_markers_across_env_and_models.R \
          analysis/gwas \
          ${trait} \
          effect \
          max
done

# Select the top 1000 non-redundant markers for each trait based on rank across environment and LD with other markers
for trait in EHT PHT Moisture YLD;do
  sbatch --export=trait=${trait} scripts/select_top_1000_non-redundant_markers.sh
done

# Select the marker results for the top 1000 non-redundant markers for each trait
# recreate an array with trait names
traits=()
for trait in $(ls -d analysis/gwas/*/); do
  #echo $trait
  traits+=$(basename $trait)
  traits+=" "
done

module load R/3.6.0
for trait in ${traits[@]}; do
  for model in additive dominant; do
    echo "${trait} - ${model}"
    # get folder name
    FOLDER=analysis/gwas/${trait}/${model}_model/permutation
    # get input files
    GWAS=${FOLDER}/GAPIT.MLM.BLUE.GWAS.Results.csv
    # get pvalue threshold used in gwas
    if [[ -e ${FOLDER}/MLM.BLUE.GWAS.Results.sig-hits.csv ]]; then
      PVAL=$(cut -f 7 -d , ${FOLDER}/MLM.BLUE.GWAS.Results.sig-hits.csv | sed 1d | uniq)
    else
      # if file doesn't exist, just set pval to a very low number
      PVAL=1e-100
    fi
    # Rscript to select the top markers and specify significance
    Rscript scripts/select_top_markers_in_each_env.R \
            ${GWAS} \
            ${PVAL} \
            analysis/gwas \
            ${trait} \
            ${model} \
            ${FOLDER}            
  done
done

# Combine all of the top markers gwas data into one file for various numbers of markers
# This will take a while to run, especially as the number of markers increases.
for nmarker in 1 3 5 10 20 30 40 50 60 70 80 90 100 150 200 250 300 350; do
  for trait in EHT PHT Moisture YLD; do
    # get output folder
    OUTFOLDER=analysis/sim_traits/${trait}/gwas_markers/
    # Create a tab separated header for the output file containing env, marker, chr, pos, maf, effect, type, pval, and signif
    echo -e "env\tmarker\tchr\tpos\tmaf\teffect\ttype\tpval\tsignif" > ${OUTFOLDER}/gwas_top-${nmarker}_markers.txt
    # State the iteration
    echo "${trait} - ${nmarker}"
    # Get list of markers
    head -n ${nmarker} analysis/gwas/${trait}_top_markers.txt > temp_markers.txt
    for marker in $(cat temp_markers.txt); do
      #echo "${trait} - ${nmarker} - ${marker}" # Debugging
      # get folder names for the trait (environments) to loop through
      for env in $(ls -d analysis/gwas/${trait}*/); do
        # Loop through each model type
        for model in additive dominant; do
          # get folder name with data
          FOLDER=${env}${model}_model/permutation
          # get input files
          GWASFILE=${FOLDER}/top_markers_gwas_data.txt
          # search for the listed markers in the relevant gwas results files
          grep -w ${marker} ${GWASFILE} >> ${OUTFOLDER}/gwas_top-${nmarker}_markers.txt
        done
      done
    done
  done
done
rm temp_markers.txt
```




#### Percent variance explained

To have a better idea how well the GWAS hits capture the phenotypic variation for each trait, I ran [GCTA](https://yanglab.westlake.edu.cn/software/gcta/) to calculate percent variance explained (PVE) for all markers used in GWAS and for all significant GWAS hits.

```bash
# cd ~/software
# # GCTA installation
# wget https://yanglab.westlake.edu.cn/software/gcta/bin/gcta_v1.94.0Beta_linux_kernel_2_x86_64.zip
# unzip gcta_v1.94.0Beta_linux_kernel_2_x86_64.zip
# # create a symbolic link to make it easier to write code
# cd ~/software/gcta_v1.94.0Beta_linux_kernel_2_x86_64
# ln -s gcta_v1.94.0Beta_linux_kernel_2_x86_64_static gcta
# # added gcta to PATH in ~/.bashrc to launch it anywhere in the command line

cd ~/empirical_sim

# create folder to keep files from GCTA analysis
mkdir -p analysis/gcta

# adjust phenotype files to gcta format
for pheno in data/1stStage_BLUEs.*_*.txt; do
  trait=$(basename ${pheno} .txt | cut -d '.' -f 2)
  mkdir -p analysis/gcta/${trait}
  sed 1d ${pheno} | sed -e 's/^/-9\t/' > analysis/gcta/${trait}/${trait}.phen
done

# calculate pve explained by all markers
HMP=data/usda_hybrids_projected-SVs-SNPs.poly.low-missing.pruned-10.geno-miss-0.25.hmp.txt
FOLDER=analysis/gcta
sbatch --export=HMP=${HMP},FOLDER=${FOLDER} scripts/gcta_pve_all-markers.sh
```

Overall, percent variance explained using all markers from GWAS showed that less than half of the phenotypic variation can be explained by these markers, suggesting that we are missing a lot of markers not in LD with causative variants. In addition, most of GWAS hits explain very little of phenotypic variation, which was expected: there were only a few GWAS hits that were significant and these traits are highly quantitative. Because of that, we will add the top 100, non-redudant non-significant markers to the significant markers and do the GCTA analysis again. The goal here is to see if adding these markers increase the percent variance explained.

```bash
# create an array with trait names
traits=()
for trait in $(ls -d analysis/gwas/*/); do
  #echo $trait
  traits+=$(basename $trait)
  traits+=" "
done

# get list of significant markers with top non-signficant markers
NMARKERS=100
for trait in ${traits[@]}; do
  echo ${trait}
  for model in additive dominant; do
    # get folder path and gwas filename
    FOLDER=analysis/gwas/${trait}/${model}_model/permutation
    GWASFILE=${FOLDER}/top-gwas-peaks_non-redundant_markers.txt
    OUT=${FOLDER}/list_top_non-redudant_gwas_markers.txt
    # create marker list
    sed 1d ${GWASFILE} | sort -g -k 6,6 | cut -f 1 | head -n ${NMARKERS} > ${OUT}
  done
done

# calculate pve explained by significant and top non-significant markers only
HMP=data/usda_hybrids_projected-SVs-SNPs.poly.low-missing.pruned-10.geno-miss-0.25.hmp.txt
GWASFOLDER=analysis/gwas
GCTAFOLDER=analysis/gcta
sbatch --export=HMP=${HMP},GWASFOLDER=${GWASFOLDER},GCTAFOLDER=${GCTAFOLDER} scripts/gcta_pve_sig-plus-non-sig-markers.sh

# plot PVE all markers vs significant markers
module load R/3.6.0
Rscript scripts/plot_pve_gcta.R analysis/gcta
```

Overall, the PVE increased the more markers we added, but in a few cases it actually decreased or increased in a non-additive way.



## Simulations

The main idea of this project is to simulate traits based on the genetic architecture of traits phenotyped in the field and then run genomic prediction models on them. Sharon will then run the same prediction models on the real data, and we will compare the results with the ones from simulated data. We hope that the results of simulated data will help us better interpret the results we see from the real data.

To simulate traits based on genetic architecture I need to:
- Get all top non-redudant GWAS hits (based on their average rank across environments) for each trait for all environments and effect types (i.e. addititive or dominant).
- Get correlation among envirnoments with `envRtype`.
- Get the mean values of BLUEs for each environment.
- Get the heritabilities for each environment.
- Give exact same or reduced effect size (from GWAS) and type (additive or dominant) to each causative variant at each env.
- Once traits are simulated, compare the correlation (Spearman because I'm more interested in the rank changes, not much the magnitude) and the distribution of simulated values with the real data.


### Correlation among environments

To simulate phenotypic data across multiple environments in simplePHENOTYPES, I need to use a correlation matrix among environments in the `cor_res` option of `create_phenotypes()` function.

Using the R package [EnvRtype](https://doi.org/10.1093/g3journal/jkab040), I gathered weather data from 11 environments (location-year) used in the 2019-2020 USDA hybrid trials as detailed in the document `data/NIFA_PhenotypicDatasetInfo.docx` and obtained a correlation matrix across these environments. The input file was manually generated on excel (table containing 6 columns: Location, Symbol, Latitude, Longitude, Start_date, End_date) and saved as `data/usda_sites_coord_for_simulation.csv`. The matrix was generated by `scripts/get_usda_env-types.R` and the output matrix `usda_envs_cor_matrix.txt` and correlation heatmap `usda_envs_cor_matrix.pdf` were saved in the `data` folder of this simulation project.

```bash
# I had to run this locally as compatibility issues arose with EnvRtype on MSI.
# get similarity across environments
Rscript scripts/get_usda_env-types.R data/usda_sites_coord_for_simulation.csv \
                                     data/usda_envs_cor_matrix.txt \
                                     --country=USA --heatmap
# remove temp files
rm USA*_msk_alt.*
```


### Phenotypic simulations

Get means of BLUEs for all traits:

```bash
echo -e "trait\tenv\tmean" > data/1st_stage_means.txt
for file in data/1stStage_BLUEs.*_*.txt; do
  trait=$(basename -a ${file} | cut -d "." -f 2 | cut -d "_" -f 1)
  env=$(basename -a ${file} | cut -d "." -f 2 | cut -d "_" -f 2)
  mean=$(sed 1d ${file} | awk '{ total += $2 } END { print total/NR }')
  echo -e "${trait}\t${env}\t${mean}" >> data/1st_stage_means.txt
done
```

Merge phenotypic data per environment and top GWAS results (non-redudant markers) per environment:

```bash
echo -e "trait\tenv\tmean" > data/1st_stage_means.txt
for file in data/1stStage_BLUEs.*_*.txt; do
  trait=$(basename -a ${file} | cut -d "." -f 2 | cut -d "_" -f 1)
  env=$(basename -a ${file} | cut -d "." -f 2 | cut -d "_" -f 2)
  mean=$(sed 1d ${file} | awk '{ total += $2 } END { print total/NR }')
  echo -e "${trait}\t${env}\t${mean}" >> data/1st_stage_means.txt
done

# get trait names (without env)
traits=($(ls -d analysis/gwas/*/))
for (( i = 0; i < ${#traits[@]}; i++ )); do
  traits[${i}]=$(basename ${traits[${i}]} | cut -d "." -f 2 | cut -d "_" -f 1)
done
# save current internal field separator
old_IFS=${IFS}
# keep only unique elements of array
IFS=" " read -r -a traits <<< "$(tr ' ' '\n' <<< "${traits[@]}" | sort -u | tr '\n' ' ')"
# restore previous IFS
IFS=${old_IFS}

module load R/3.6.0
# merge pheno per env
for trait in ${traits[@]}; do
  echo ${trait}
  Rscript scripts/merge_pheno_per_trait.R \
          ${trait} \
          BAY19,BEC-BL19,BEC-BL20,BEC-EP20,COR19,COR20,MIN19,MIN20,SYN19,SYN20,URB19 \
          data \
          1stStage_BLUEs
done

# merge gwas results per env
for trait in ${traits[@]}; do
  echo ${trait}
  Rscript scripts/merge_gwas_per_trait.R \
          ${trait} \
          BAY19,BEC-BL19,BEC-BL20,BEC-EP20,COR19,COR20,MIN19,MIN20,SYN19,SYN20,URB19 \
          analysis/gwas \
          top-gwas-peaks_non-redundant_markers.txt
done
```

Select the top markers based on highest average rank across envs:

```bash
# set variables
HMP=data/usda_hybrids_projected-SVs-SNPs.poly.low-missing.pruned-10.geno-miss-0.25.hmp.txt
SVS=data/SVs_IDs_poly.txt
RESCOR=data/usda_envs_cor_matrix.txt
H2FILE=data/1stStage_heritabilities.csv
MEANSFILE=data/1st_stage_means.txt
# optional variables
REPS=3
IMPUTETYPE=Middle
QTNVAR=QTN-variance
# allow minor qtls
MINORQTLS=minor-qtls

# define function to sample different seeds in the loop
get_seeded_random()
{
  seed="$1"
  openssl enc -aes-256-ctr -pass pass:"$seed" -nosalt < /dev/zero 2> /dev/null
}

##### keep GWAS effects #####

# starting seed
SEED=823162
# keep minor qtl effects with original gwas effects
EFFECTRED=1

# Make sure simplePHENOTYPES is installed
### 78193754
for trait in EHT PHT Moisture YLD; do

  # define folder to save results
  FOLDER=analysis/sim_traits/${trait}
  # set remaining variables
  PHENOFILE=data/1stStage_BLUEs.${trait}-per-env.txt
  # simulate traits
  sbatch --export=TRAIT=${trait},PHENOFILE=${PHENOFILE},HMP=${HMP},SVS=${SVS},RESCOR=${RESCOR},H2FILE=${H2FILE},MEANSFILE=${MEANSFILE},FOLDER=${FOLDER},EFFECTRED=${EFFECTRED},SEED=${SEED},REPS=${REPS},IMPUTETYPE=${IMPUTETYPE},QTNVAR=${QTNVAR},MINORQTLS=${MINORQTLS} scripts/trait_simulation_hybrids_cor-curve.sh # The r script for this will be in R/4.0.4 to comply with GAPIT and simplePHENOTYPES requirements
  # wait 3 seconds to avoid accessing the same file at the same time
  sleep 3
  # get a new seed
  SEED=($(shuf -i 1-100000 -n 1 --random-source=<(get_seeded_random ${SEED})))
done

##### reduce effects of all qtls #####

# starting seed
SEED=22144
# reduce effects of all qtls
EFFECTRED=0.1
EFFECTREDALL=reduce-effect-all

for trait in EHT PHT Moisture YLD; do
  # define folder to save results
  FOLDER=analysis/sim_traits/${trait}
  # set remaining variables
  PHENOFILE=data/1stStage_BLUEs.${trait}-per-env.txt
  # simulate traits
  sbatch --export=TRAIT=${trait},PHENOFILE=${PHENOFILE},HMP=${HMP},SVS=${SVS},RESCOR=${RESCOR},H2FILE=${H2FILE},MEANSFILE=${MEANSFILE},FOLDER=${FOLDER},EFFECTRED=${EFFECTRED},SEED=${SEED},REPS=${REPS},IMPUTETYPE=${IMPUTETYPE},QTNVAR=${QTNVAR},MINORQTLS=${MINORQTLS},EFFECTREDALL=${EFFECTREDALL} scripts/trait_simulation_hybrids_cor-curve.sh # The r script for this will be in R/4.0.4 to comply with GAPIT and simplePHENOTYPES requirements
  # wait 3 seconds to avoid accessing the same file at the same time
  sleep 3
  # get a new seed
  SEED=($(shuf -i 1-100000 -n 1 --random-source=<(get_seeded_random ${SEED})))
done
```

> After simulating traits, `scripts/trait_simulation_hybrids.R` change the names of the environments to `env1`, `env2`, etc. Here's a table with the conversion for future reference:
>
> | trait    | env number | env name |
> | -------- | ---------- | -------- |
> | EHT      | env1       | COR19    |
> |          | env2       | MIN19    |
> |          | env3       | MIN20    |
> |          | env4       | URB19    |
> | PHT      | env1       | COR19    |
> |          | env2       | COR20    |
> |          | env3       | MIN19    |
> |          | env4       | MIN20    |
> |          | env5       | URB19    |
> | YLD      | env1       | BEC-BL19 |
> |          | env2       | BEC-BL20 |
> |          | env3       | BEC-EP20 |
> |          | env4       | COR19    |
> |          | env5       | COR20    |
> |          | env6       | MIN19    |
> |          | env7       | MIN20    |
> |          | env8       | SYN19    |
> |          | env9       | SYN20    |
> |          | env10      | URB19    |
> | Moisture | env1       | BAY19    |
> |          | env2       | BEC-BL19 |
> |          | env3       | BEC-BL20 |
> |          | env4       | BEC-EP20 |
> |          | env5       | COR19    |
> |          | env6       | COR20    |
> |          | env7       | MIN19    |
> |          | env8       | MIN20    |
> |          | env9       | SYN19    |
> |          | env10      | SYN20    |
> |          | env11      | URB19    |




### QC simulations

Run ANOVA to confirm that simulated traits have GxE and also plot the percent variance explained by the QTLs controlling simulated trait variation.

```bash
METHOD=michael_method
for trait in EHT PHT Moisture YLD; do
  for effect in 1 0.1; do
    # define folder with simulated data
    FOLDER=analysis/sim_traits/${trait}
    # run anova and plot pve
    sbatch --export=FOLDER=${FOLDER},METHOD=${METHOD},EFFECT=${effect} scripts/anova_sim_traits_hybrids_cor-curve.sh # had to install MuMIn 1.40.4 through CRAN archive to stick with R version 3.6.0
  done
done
```

Display the significant markers used to simulate each trait

```bash
for trait in EHT PHT Moisture YLD; do
  echo "---${trait} Sig Hits---"
  awk '$9 == "TRUE"' analysis/sim_traits/${trait}/gwas_markers/gwas_top-350_markers.txt | cut -f 1-2
done
```


### Comparison to real data

Sharon shared all the empirical data with with me and can be found in `data/NIFA_CompleteDataset.xlsx`. This file was then transferred to my MSI account and then re-formated to a tab-delimited file for compatibility with other scripts using `scripts/format_empirical-data.R`. This script creates one file per trait with a particular set of environments: in this case, the environments match the order used when simulating data.

```bash
module load R/3.6.0

# EHT
Rscript scripts/format_empirical-data.R \
        data/NIFA_CompleteDataset.xlsx \
        EHT \
        --envs=COR19,MIN19,MIN20,URB19 \
        --checks=UIUC-CHECK1,UIUC-CHECK2,UIUC-CHECK3,UIUC-CHECK4,6049V2P

# PHT
Rscript scripts/format_empirical-data.R \
        data/NIFA_CompleteDataset.xlsx \
        PHT \
        --envs=COR19,COR20,MIN19,MIN20,URB19 \
        --checks=UIUC-CHECK1,UIUC-CHECK2,UIUC-CHECK3,UIUC-CHECK4,6049V2P

# YLD
Rscript scripts/format_empirical-data.R \
        data/NIFA_CompleteDataset.xlsx \
        YLD \
        --envs=BEC-BL19,BEC-BL20,BEC-EP20,COR19,COR20,MIN19,MIN20,SYN19,SYN20,URB19 \
        --checks=UIUC-CHECK1,UIUC-CHECK2,UIUC-CHECK3,UIUC-CHECK4,6049V2P

# Moisture
Rscript scripts/format_empirical-data.R \
        data/NIFA_CompleteDataset.xlsx \
        Moisture \
        --envs=BAY19,BEC-BL19,BEC-BL20,BEC-EP20,COR19,COR20,MIN19,MIN20,SYN19,SYN20,URB19 \
        --checks=UIUC-CHECK1,UIUC-CHECK2,UIUC-CHECK3,UIUC-CHECK4,6049V2P
```

Run ANOVA on real traits and calculate percent variance explained:

```bash
module load R/3.6.0

# run anova on real traits
for trait in EHT PHT Moisture YLD; do
  PHENO=data/NIFA_CompleteDataset.${trait}.txt
  sbatch --export=PHENO=${PHENO} scripts/anova_real_traits_hybrids.sh
done
```




### Correlation curve

Plot correlation between simulated and empirical traits and percent variance explained from ANOVA for each simulated trait:

```bash
module load R/4.2.2-openblas
# first summarize results
Rscript scripts/qc_sim_correlation_curve.R \
        analysis/sim_traits \
        EHT,PHT,Moisture,YLD \
        3,5,10,20,30,40,50,60,70,80,90,100,150,200,250,300,350 \
        0.1,1 \
        --avg-rank-markers=all

# plot results
Rscript scripts/plot_sim_correlation_curve.R \
        analysis/sim_traits \
        analysis/sim_traits/summary_cor_across-within_envs.txt \
        analysis/sim_traits/summary_stats_across_envs.txt \
        analysis/sim_traits/summary_stats_within_envs.txt \
        analysis/sim_traits/summary_pve.txt
```




## Genomic prediction

We will run genomic prediction with 500 random markers (sampled 3 times) to predict simulated and empirical traits.


### Create datasets

```bash
# function to get seed number
get_seeded_random()
{
  seed="$1"
  openssl enc -aes-256-ctr -pass pass:"$seed" -nosalt < /dev/zero 2> /dev/null
}

# select markers from pruned hybrids genotypic data
# bootstrap 10 times
HMP=data/usda_hybrids_projected-SVs-SNPs.poly.low-missing.pruned-10.geno-miss-0.25.hmp.txt
NMARKERS=500
SEED=72162
for iter in {1..10}; do
  # create folder to save iterations
  FOLDER=analysis/sim_traits/predictors/iter${iter}
  mkdir -p ${FOLDER}
  # simply randomly sample markers
  shuf -n ${NMARKERS} --random-source=<(get_seeded_random ${SEED}) <(cut -f 1 ${HMP} | sed 1d) > ${FOLDER}/marker_list.iter${iter}.txt
  # filter hmp to generate predictor dataset
  run_pipeline.pl -Xmx10g -importGuess ${HMP} \
                  -includeSiteNamesInFile ${FOLDER}/marker_list.iter${iter}.txt \
                  -export ${FOLDER}/hybrids_geno.iter${iter}.hmp.txt \
                  -exportType HapmapDiploid > /dev/null
  # change seed
  SEED=$((${SEED}+5))
done
```


### Predict simulated traits

I will be using ASREML-R to run predictions, but the license I have doesn't allow me to run things in the cloud (or MSI). So I will copy all the files from the simulations to our lab Linux server.

> Marker datasets and scripts were already transfered.

```bash
# Transfer data to local computer for asreml use
rsync -r burns756@agate.msi.umn.edu:/home/hirschc1/burns756/empirical_sim/analysis/sim_traits/ analysis/sim_traits
# simulated traits from GWAS and marker datasets for prediction
# used rsync -r to transfer data in analysis/sim_traits to local computer (~715 Mb)
# Don't forget to transfer the scripts!
```

I'm doing a two-stage genomic prediction analysis, where in the first stage I extract BLUEs from each environment separately:

```bash
### RUN LOCALLY ###
# get blues for each environment
# This will take a few minutes to complete
nohup bash scripts/genomic_prediction_stage1_hybrids.sh -t sim > analysis/sim_traits/genomic_prediction_stage1.log 2>&1 &
echo $! > analysis/sim_traits/nohup_pid_stage1.txt
# # if need to kill process, run:
# kill -9 `analysis/sim_traits/nohup_pid_stage1.txt`
# rm analysis/sim_traits/nohup_pid_stage1.txt
```

Then I use these BLUEs in a second-stage as the phenotypes for my genomic predictions across multiple environments with different variance structures. Based on prelimary analysis, simulated traits with many environments like yield and moisture took a long time to run this second stage. Thus, we decided to first remove environments that were very similar to each other based on their correlation matrix (`data/usda_envs_cor_matrix.pdf` - COR19, COR20, MIN19, MIN20, SYN19, SYN20) and then run the second stage for those traits.

> Environments kept for YLD were env4-9 and for Moisture were env5-10

```bash
# This runs very fast (< 10 seconds)
for trait in YLD Moisture; do
  for nmarkers in 5 10 20 30 40 50 60 70 80 90 100 150 200 250 300 350; do
    for effect in 1 0.1; do

      # get folder name
      FOLDER=analysis/sim_traits/${trait}/traits/michael_method/n_markers_${nmarkers}/effects_${effect}/
      # check if sim file exists first
      FILE=$(echo ${FOLDER}/blues_1st_stage.txt)
      if [[ -e ${FILE} ]]; then
        # select envs to keep
        [[ ${trait} == "YLD" ]] && ENVS=env4,env5,env6,env7,env8,env9
        [[ ${trait} == "Moisture" ]] && ENVS=env5,env6,env7,env8,env9,env10
        # filter yield and moisture simulated data to have only six most different environments
        echo ${FOLDER}
        Rscript scripts/filter_envs_1st_stage.R \
                ${FOLDER}/blues_1st_stage.txt \
                ${ENVS} &
        Rscript scripts/filter_envs_1st_stage.R \
                ${FOLDER}/blues_1st_stage_weighted.txt \
                ${ENVS} &
        wait
      fi

    done
  done
done
```

#### Additive or dominant models

Generate the commands that will be used to run the second stage of the genomic prediction.

```bash
# This takes a few seconds
# generate gnu parallel commands
for iter in {1..3}; do
  for effect in A D; do

    # BLUEs without env weights
    # G structure - diag / R structure - us
    bash scripts/genomic_prediction_stage2_hybrids.sh \
    -d ${iter} -e ${effect} -g fa -G diag -i 3 -r diag -R diag \
    -o scripts/commands_stage2-prediction.cor-curve.sim-traits.iter${iter}.no-AD.txt -t sim

  done
done
```

Run models with 3 different iterations of genotypic data:

```bash
# Be sure to have GNU parallel installed on your local computer:
# brew install parallel

# Transfer files to other computers for faster processing
# Get IP Address of computers to transfer files to (ifconfig | grep 'inet ')
# rsync -r analysis/sim_traits/ rafael@10.129.33.47:michael/empirical_sim/analysis/sim_traits
# rsync -r analysis/sim_traits/ hirschlab@10.129.183.73:michael/empirical_sim/analysis/sim_traits
# rsync -r analysis/sim_traits/ hirschlab@10.129.103.27:michael/empirical_sim/analysis/sim_traits
# Transer the scripts as well
# rsync -r scripts/ rafael@10.129.33.47:michael/empirical_sim/scripts
# rsync -r scripts/ hirschlab@10.129.183.73:michael/empirical_sim/scripts
# rsync -r scripts/ hirschlab@10.129.103.27:michael/empirical_sim/scripts

# run predictions across all environments using blups from previous stage
# Set iter value
# This will take a considerable amount of time (~2 days if running across 3 computers - one for each iteration)
iter=1
nohup parallel --jobs 2 --joblog analysis/sim_traits/genomic_prediction_stage2.cor-curve.sim-traits.iter${iter}.no-AD.parallel.log < scripts/commands_stage2-prediction.cor-curve.sim-traits.iter${iter}.no-AD.txt 2> analysis/sim_traits/genomic_prediction_stage2.cor-curve.sim-traits.iter${iter}.no-AD.log &

# to kill job -- get PID then send TERM signal twice
# use negative PID number to kill all processes
# ps xw | grep parallel
# kill -TERM -36087
# kill -TERM -36087
```

> To check progress of how many predictions were completed:

  ```bash
  # To check progress, use the following line of code (should contain one line for each GNU command script line (+1 for the header) when complete)
less analysis/sim_traits/genomic_prediction_stage2.cor-curve.sim-traits.iter${iter}.no-AD.parallel.log | wc -l
  ```


Finally, transfer the files back to MSI for analysis

```bash
# Transfer files back to MSI
# scp -r hirschlab@10.129.107.244:/Users/hirschlab/michael/empirical_sim/analysis/sim_traits/* analysis/sim_traits/
# scp -r hirschlab@10.129.4.197:/Users/hirschlab/michael/empirical_sim/analysis/sim_traits/* analysis/sim_traits/

folder_base=analysis/sim_traits
traits=EHT,PHT,Moisture,YLD
sim_n_markers=5,10,20,30,40,50,60,70,80,90,100,150,200,250,300,350
sim_effects=1,0.1
outfile_name=prediction_results.sim.iter1-3.no-AD.txt

# Summarize the prediction accuracies with a correlation curve
Rscript scripts/summarize_prediction_accuracies_cor-curve.R \
        ${folder_base} \
        ${traits} \
        ${sim_n_markers} \
        ${sim_effects} \
        ${outfile_name} \
        --pred-iters=3 \
        --models=A,D \
        --var-structures=Gstr-diag_Rstr-us

# plot results
for trait in YLD Moisture EHT PHT; do
  Rscript scripts/plot_prediction_accuracy_cor-curve.R \
          analysis/sim_traits/prediction_results.sim.iter1-3.no-AD.${trait}.summary.txt \
          analysis/sim_traits/prediction_plots \
          --error-bars=CI_accuracy
done
```



### Predict real traits

After running different prediction models for simulated data, I'll now run the same models on empirical data collected for the same environments.

I will be using ASREML-R to run predictions, but the license I have doesn't allow me to run things in the cloud (or MSI). So I will transfer the files above to our lab Linux server.

> Marker datasets were already transfered when analyzing simulated data.

```bash
# Do similar predictions for the real trait data
scp -r burns756@agate.msi.umn.edu:empirical_sim/data/NIFA_CompleteDataset.*.txt data/
```

I'm doing a two-stage genomic prediction analysis, where in the first stage I extract BLUEs from each environment separately:

```bash
# create folder to save results
mkdir -p analysis/empirical_traits/{EHT,PHT,YLD,Moisture}

# get blues for each environment
nohup bash scripts/genomic_prediction_stage1_hybrids.sh -t real > analysis/empirical_traits/genomic_prediction_stage1.log 2>&1 &
echo $! > analysis/empirical_traits/nohup_pid_stage1.txt
# # if need to kill process, run:
# kill -9 `analysis/empirical_traits/nohup_pid_stage1.txt`
# rm analysis/empirical_traits/nohup_pid_stage1.txt
```

Then I use these BLUEs in a second-stage as the phenotypes for my genomic predictions across multiple environments with different variance structures. Based on prelimary analysis, simulated traits with many environments like yield and moisture took a long time to run this second stage. Thus, we decided to first remove environments that were very similar to each other based on their correlation matrix (`data/usda_envs_cor_matrix.pdf` - COR19, COR20, MIN19, MIN20, SYN19, SYN20) and then run the second stage for those traits.

> Environments kept for YLD were env4-9 and for Moisture were env5-10

```bash
# filter yield and moisture simulated data to have only six most different environments
Rscript scripts/filter_envs_1st_stage.R analysis/empirical_traits/YLD/blues_1st_stage.txt \
                                        env4,env5,env6,env7,env8,env9
Rscript scripts/filter_envs_1st_stage.R analysis/empirical_traits/YLD/blues_1st_stage_weighted.txt \
                                        env4,env5,env6,env7,env8,env9
Rscript scripts/filter_envs_1st_stage.R analysis/empirical_traits/Moisture/blues_1st_stage.txt \
                                        env5,env6,env7,env8,env9,env10
Rscript scripts/filter_envs_1st_stage.R analysis/empirical_traits/Moisture/blues_1st_stage_weighted.txt \
                                        env5,env6,env7,env8,env9,env10
```

#### Additive or dominant models

Run models with 3 different iterations of genotypic data:

```bash
# generate gnu parallel commands
for iter in {1..3}; do
  for effect in A D; do

    # BLUEs without env weights
    # G structure - diag / R structure - us
    bash scripts/genomic_prediction_stage2_hybrids.sh \
    -d ${iter} -e ${effect} -g fa -G diag -i 3 -r diag -R diag \
    -o scripts/commands_stage2-prediction.cor-curve.real.iter${iter}.no-AD.txt -t real

  done
done

# Transfer files to macbooks
# rsync -r analysis/empirical_traits/ hirschlab@10.129.183.73:michael/empirical_sim/analysis/empirical_traits
# rsync -r analysis/empirical_traits/ hirschlab@10.129.103.27:michael/empirical_sim/analysis/empirical_traits

# Run the following commands on a local computer. Change iter value on each computer to run iter 1, 2, and 3.
# This will only run 17 jobs per iteration.
iter=1
nohup parallel --jobs 2 --joblog analysis/empirical_traits/genomic_prediction_stage2.cor-curve.real.iter${iter}.no-AD.parallel.log < scripts/commands_stage2-prediction.cor-curve.real.iter${iter}.no-AD.txt 2> analysis/empirical_traits/genomic_prediction_stage2.cor-curve.real.iter${iter}.no-AD.log &
```

Finally, summarize and plot prediction results:

```bash
# Transfer files back to personal macbook
# rsync -r hirschlab@10.129.103.27:michael/empirical_sim/analysis/ analysis 
# rsync -r hirschlab@10.129.183.73:michael/empirical_sim/analysis/ analysis 

# Transfer files back to MSI
# mkdir -p analysis/empirical_traits
# rsync -r analysis/ burns756@agate.msi.umn.edu:empirical_sim/analysis
module load R/4.3.0-openblas-rocky8
folder_base=analysis
traits=EHT,PHT,Moisture,YLD
sim_n_markers=5,10,20,30,40,50,60,70,80,90,100,150,200,250,300,350
sim_effects=1,0.1
outfile_name=prediction_results
Rscript scripts/summarize_prediction_accuracies_cor-curve.R \
        ${folder_base} \
        ${traits} \
        ${sim_n_markers} \
        ${sim_effects} \
        ${outfile_name} \
        --pred-iters=3 \
        --models=A,D \
        --var-structures=Gstr-fa_Rstr-diag

# plot results
for trait in YLD Moisture EHT PHT; do
  Rscript scripts/plot_prediction_accuracy_cor-curve.R \
          analysis/empirical_traits/prediction_results.real.iter1-3.no-AD.${trait}.summary.txt \
          analysis/empirical_traits/prediction_plots \
          --error-bars=CI_accuracy
done
```


## Scripts used for manuscript curation:
```bash
# To generate the tables (note: tables may not be named the exact same as their manuscript table ID if they change places):
Rscript scripts/manuscript_tables.R

# To generate the figures (note: figures may not be named the exact same as their manuscript figure ID if they change places):
Rscript scripts/manuscript_figures.R
```

