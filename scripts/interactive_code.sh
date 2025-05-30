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

# qc
sbatch --export=HMP=data/usda_hybrids_projected-SVs-SNPs.poly.maf-filter.hmp.txt,SVS=data/SVs_IDs_poly.txt,FOLDER=analysis/qc/snp-sv_hmp scripts/qc_snp-sv_hmp.sh

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

# Again run on MSI
# merge heritability files
Rscript scripts/merge_h2_files.R data 1stStage_heritabilities

# Make sure the following are installed
#r <- getOption("repos");
#r["CRAN"] <- "http://cran.rstudio.com/";
#options(repos=r);
#install.packages(c('BiocManager', "gplots", "ape", "EMMREML", "scatterplot3d", "bigmemory", "biganalytics"));


#if (!require("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install(c('multtest', 'genetics', 'snpStats', 'compiler'))

module load R/3.6.0

# step 1 - I had to download the gapit source code to my MSI directory to get this to work
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

# Remove redundant SNPs from GWAS results
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

#for trait in ${traits[@]}; do
#  echo "--- ${trait} ---"
#  for model in additive dominant; do
#    FOLDER=analysis/gwas/${trait}/${model}_model/permutation
#    if [[ -e ${FOLDER}/MLM.BLUE.GWAS.Results.sig-hits.csv ]]; then
#        # Add the header to the output file
#        head -n1 ${FOLDER}/GAPIT.MLM.BLUE.GWAS.Results.csv > ${FOLDER}/list_sig_all_non-sig_markers.txt
#        # Get a list of significant markers
#        cut -f1-3 -d',' ${FOLDER}/MLM.BLUE.GWAS.Results.sig-hits.csv > ${FOLDER}/tmp_sig_markers_list.csv
#        # Add significant markers to the output file
#        tail -n+2 ${FOLDER}/GAPIT.MLM.BLUE.GWAS.Results.csv | grep -f ${FOLDER}/tmp_sig_markers_list.csv | sort -g -t',' -k4 | cut -f1 -d',' >> ${FOLDER}/list_sig_all_non-sig_markers.txt
#        # Add non-significant markers to the output file
#        tail -n+2 ${FOLDER}/GAPIT.MLM.BLUE.GWAS.Results.csv | grep -v -f ${FOLDER}/tmp_sig_markers_list.csv | sort -g -t',' -k4 | cut -f1 -d',' >> ${FOLDER}/list_sig_all_non-sig_markers.txt
#        # Remove temporary file
#        rm ${FOLDER}/tmp_sig_markers_list.csv; else
#        # Add the header to the output file
#        head -n1 ${FOLDER}/GAPIT.MLM.BLUE.GWAS.Results.csv > ${FOLDER}/list_sig_all_non-sig_markers.txt
#        # Add non-significant markers to the output file
#        tail -n+2 ${FOLDER}/GAPIT.MLM.BLUE.GWAS.Results.csv | sort -g -t',' -k4 | cut -f1 -d',' >> ${FOLDER}/list_sig_all_non-sig_markers.txt
#    fi
#  done
#  echo ""
#done

# Check file sizes:
for trait in ${traits[@]}; do
  echo "--- ${trait} ---"
  for model in additive dominant; do
    #echo $(wc -l analysis/gwas/${trait}/${model}_model/permutation/list_top_5000_non-sig_markers.txt)
    #echo $(wc -l analysis/gwas/${trait}/${model}_model/permutation/list_top_10000_non-sig_markers.txt)
    echo $(wc -l analysis/gwas/${trait}/${model}_model/permutation/list_sig_all_non-sig_markers.txt)
  done
done

# You may need to re-run the traits list generation above depending on when you run this in relation to the last segment.
# calculate LD among sig and all non-sig markers
# THIS MAY BE REPLACED WITH A LOOPED SHELL SCRIPT IN THE FUTURE
#module load plink/1.90b6.10
#HMP=data/usda_hybrids_projected-SVs-SNPs.poly.low-missing.pruned-10.geno-miss-0.25.hmp.txt
#for trait in ${traits[@]}; do
#  echo "--- ${trait} ---"
#  for model in additive dominant; do
#    # get folder name
#    FOLDER=analysis/gwas/${trait}/${model}_model/permutation
#    # get list of markers to calcukate LD
#    NONSIGMARKERS=${FOLDER}/list_sig_all_non-sig_markers.txt
#    # define plk filename
#    PLK=${FOLDER}/usda_hybrids_projected-SVs-SNPs.poly.low-missing.pruned-10.geno-miss-0.25.topall-non-sig-hits
#    # define ld output name
#    OUT=${FOLDER}/ld_topall_non-sig-markers
#    # hmp2plk
#    run_pipeline.pl -Xmx40g -importGuess ${HMP} \
#                            -includeSiteNamesInFile ${NONSIGMARKERS} \
#                            -export ${PLK} \
#                            -exportType Plink > /dev/null
#    # calculate LD within same chromosome
#    plink --file ${PLK}.plk \
#          --make-founders \
#          --r2 gz dprime with-freqs \
#          --ld-window-r2 0 \
#          --ld-window 400000000 \
#          --ld-window-kb 400000 \
#          --out ${OUT}
#  done
#  echo ""
#  echo $(wc -l ${OUT}.ld.gz)
#done

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


# determine how many top total markers will be used
#TOTALMARK=1000
# determine LD threshold
#LD=0.9
# remove redundant markers
#sbatch --export=NONSIG=${TOTALMARK},LD=${LD} scripts/remove_redundant_sig_top-non-sig_markers.sh

# Calculate percent variance explained
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

# Run EnvRtype locally...

# Get means of blues for all traits
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

# create new folder for simulations
mkdir -p analysis/sim_traits

# first need to calculate LD between all markers
#HMP=data/usda_hybrids_projected-SVs-SNPs.poly.low-missing.pruned-10.geno-miss-0.25.hmp.txt
#module load plink/1.90b6.10
#for trait in EHT PHT Moisture YLD; do

#  # create new folder
#  FOLDER=analysis/sim_traits/${trait}/ld
#  mkdir -p ${FOLDER}
#  # get list of markers to calcukate LD
#  cut -f 2 analysis/gwas/gwas_top-peaks_${trait}-per-env.txt | sed 1d | sort -u > ${FOLDER}/list_top_markers.txt
#  # define plk filename
#  PLK=${FOLDER}/usda_hybrids_projected-SVs-SNPs.poly.low-missing.pruned-10.geno-miss-0.25.top-markers
#  # define ld output name
#  OUT=${FOLDER}/ld_top-markers
#  # hmp2plk
#  run_pipeline.pl -Xmx40g -importGuess ${HMP} \
#                          -includeSiteNamesInFile ${FOLDER}/list_top_markers.txt \
#                          -export ${PLK} \
#                          -exportType Plink > /dev/null
#  # calculate LD within same chromosome
#  plink --file ${PLK}.plk \
#        --make-founders \
#        --r2 gz dprime with-freqs \
#        --ld-window-r2 0 \
#        --ld-window 400000000 \
#        --ld-window-kb 400000 \
#        --out ${OUT}
#
#done

# select markers
#LD=0.9
#for trait in EHT PHT Moisture YLD; do
#  # define variables
#  GWASFILE=analysis/gwas/gwas_top-peaks_${trait}-per-env.txt
#  FOLDER=analysis/sim_traits/${trait}
#  LDFILE=${FOLDER}/ld/ld_top-markers.ld.gz
#  # submit job
#  sbatch --export=FOLDER=${FOLDER},GWASFILE=${GWASFILE},LDFILE=${LDFILE},LD=${LD} scripts/select_top_gwas_markers_across_envs.sh
#done

# Simulate phenotypes!
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

METHOD=michael_method
for trait in EHT PHT Moisture YLD; do
  for effect in 1 0.1; do
    # define folder with simulated data
    FOLDER=analysis/sim_traits/${trait}
    # run anova and plot pve
    sbatch --export=FOLDER=${FOLDER},METHOD=${METHOD},EFFECT=${effect} scripts/anova_sim_traits_hybrids_cor-curve.sh # had to install MuMIn 1.40.4 through CRAN archive to stick with R version 3.6.0
  done
done

# Display all of the significant markers for each trait
for trait in EHT PHT Moisture YLD; do
  echo "---${trait} Sig Hits---"
  awk '$9 == "TRUE"' analysis/sim_traits/${trait}/gwas_markers/gwas_top-350_markers.txt | cut -f 1-2
done

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

# run anova on real traits
for trait in EHT PHT Moisture YLD; do
  PHENO=data/NIFA_CompleteDataset.${trait}.txt
  sbatch --export=PHENO=${PHENO} scripts/anova_real_traits_hybrids.sh
done

module load R/4.3.0-openblas-rocky8
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

# Transfer data to local computer for asreml use
rsync -r burns756@agate.msi.umn.edu:/home/hirschc1/burns756/empirical_sim/analysis/sim_traits/ analysis/sim_traits
# simulated traits from GWAS and marker datasets for prediction
# used rsync -r to transfer data in analysis/sim_traits to local computer (~715 Mb)
# Don't forget to transfer the scripts!

### RUN LOCALLY ###
# get blues for each environment
# This will take a few minutes to complete
nohup bash scripts/genomic_prediction_stage1_hybrids.sh -t sim > analysis/sim_traits/genomic_prediction_stage1.log 2>&1 &
echo $! > analysis/sim_traits/nohup_pid_stage1.txt
# # if need to kill process, run:
# kill -9 `analysis/sim_traits/nohup_pid_stage1.txt`
# rm analysis/sim_traits/nohup_pid_stage1.txt

# This runs very fast (< 10 seconds)
#for trait in YLD Moisture; do
#  for nmarkers in 5 10 20 30 40 50 60 70 80 90 100 150 200 250 300 350; do
#    for effect in 1 0.1; do
#
#      # get folder name
#      FOLDER=analysis/sim_traits/${trait}/traits/michael_method/n_markers_${nmarkers}/effects_${effect}/
#      # check if sim file exists first
#      FILE=$(echo ${FOLDER}/blues_1st_stage.txt)
#      if [[ -e ${FILE} ]]; then
#        # select envs to keep
#        [[ ${trait} == "YLD" ]] && ENVS=env4,env5,env6,env7,env8,env9
#        [[ ${trait} == "Moisture" ]] && ENVS=env5,env6,env7,env8,env9,env10
#        # filter yield and moisture simulated data to have only six most different environments
#        echo ${FOLDER}
#        Rscript scripts/filter_envs_1st_stage.R \
#                ${FOLDER}/blues_1st_stage.txt \
#                ${ENVS} &
#        Rscript scripts/filter_envs_1st_stage.R \
#                ${FOLDER}/blues_1st_stage_weighted.txt \
#                ${ENVS} &
#        wait
#      fi
#
#    done
#  done
#done


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

# Install GNU parallel on local computer
# brew install parallel


# Transfer files to other computers for faster processing
# Get IP Address of computers to transfer files to (ifconfig | grep 'inet ')
# rsync -r analysis/sim_traits/ rafael@10.129.33.47:michael/empirical_sim/analysis/sim_traits
# rsync -r analysis/sim_traits/ hirschlab@10.129.2.191:michael/empirical_sim/analysis/sim_traits
# rsync -r analysis/sim_traits/ hirschlab@10.129.25.47:michael/empirical_sim/analysis/sim_traits
# Transer the scripts as well
# rsync -r scripts/ rafael@10.129.33.47:michael/empirical_sim/scripts
# rsync -r scripts/ hirschlab@10.129.2.191:michael/empirical_sim/scripts
# rsync -r scripts/ hirschlab@10.129.25.47:michael/empirical_sim/scripts

# run predictions across all environments using blups from previous stage
# Set iter value
# This will run 536 jobs per iteration.
iter=1
nohup parallel --jobs 2 --joblog analysis/sim_traits/genomic_prediction_stage2.cor-curve.sim-traits.iter${iter}.no-AD.parallel.log < scripts/commands_stage2-prediction.cor-curve.sim-traits.iter${iter}.no-AD.txt 2> analysis/sim_traits/genomic_prediction_stage2.cor-curve.sim-traits.iter${iter}.no-AD.log &
# [1] 36087

# # to kill job -- get PID then send TERM signal twice
# # use negative PID number to kill all processes
# ps xw | grep parallel
# kill -TERM -36087
# kill -TERM -36087

# To check progress, use the following line of code (should be 672 when complete)
less analysis/sim_traits/genomic_prediction_stage2.cor-curve.sim-traits.iter${iter}.no-AD.parallel.log | wc -l

# Transfer files back to MSI
# scp -r hirschlab@10.129.107.244:/Users/hirschlab/michael/empirical_sim/analysis/sim_traits/* analysis/sim_traits/
# scp -r hirschlab@10.129.4.197:/Users/hirschlab/michael/empirical_sim/analysis/sim_traits/* analysis/sim_traits/
folder_base=analysis/sim_traits
traits=EHT,PHT,Moisture,YLD
sim_n_markers=5,10,20,30,40,50,60,70,80,90,100,150,200,250,300,350
sim_effects=1,0.1
outfile_name=prediction_results.sim.iter1-3.no-AD.txt
Rscript scripts/summarize_prediction_accuracies_cor-curve.R \
        ${folder_base} \
        ${traits} \
        ${sim_n_markers} \
        ${sim_effects} \
        ${outfile_name} \
        --pred-iters=3 \
        --models=A,D \
        --var-structures=Gstr-diag_Rstr-us

mkdir -p analysis/sim_traits_gwas/cor_curve
Rscript scripts/qc_pred_correlation_curve.R ....

# plot results
for trait in YLD Moisture EHT PHT; do
  Rscript scripts/plot_prediction_accuracy_cor-curve.R \
          analysis/sim_traits/prediction_results.sim.iter1-3.no-AD.${trait}.summary.txt \
          analysis/sim_traits/prediction_plots \
          --error-bars=CI_accuracy
done

# Do similar predictions for the real trait data
scp -r burns756@agate.msi.umn.edu:empirical_sim/data/NIFA_CompleteDataset.*.txt data/

# create folder to save results
mkdir -p analysis/empirical_traits/{EHT,PHT,YLD,Moisture}

# get blues for each environment
nohup bash scripts/genomic_prediction_stage1_hybrids.sh -t real > analysis/empirical_traits/genomic_prediction_stage1.log 2>&1 &
echo $! > analysis/empirical_traits/nohup_pid_stage1.txt
# # if need to kill process, run:
# kill -9 `analysis/empirical_traits/nohup_pid_stage1.txt`
# rm analysis/empirical_traits/nohup_pid_stage1.txt

Rscript scripts/filter_envs_1st_stage.R analysis/empirical_traits/YLD/blues_1st_stage.txt \
                                        env4,env5,env6,env7,env8,env9
Rscript scripts/filter_envs_1st_stage.R analysis/empirical_traits/YLD/blues_1st_stage_weighted.txt \
                                        env4,env5,env6,env7,env8,env9
Rscript scripts/filter_envs_1st_stage.R analysis/empirical_traits/Moisture/blues_1st_stage.txt \
                                        env5,env6,env7,env8,env9,env10
Rscript scripts/filter_envs_1st_stage.R analysis/empirical_traits/Moisture/blues_1st_stage_weighted.txt \
                                        env5,env6,env7,env8,env9,env10

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

mkdir -p analysis/sim_traits_gwas/cor_curve
Rscript scripts/qc_pred_correlation_curve.R ....

# plot results
for trait in YLD Moisture EHT PHT; do
  Rscript scripts/plot_prediction_accuracy_cor-curve.R \
          analysis/empirical_traits/prediction_results.real.iter1-3.no-AD.${trait}.summary.txt \
          analysis/empirical_traits/prediction_plots \
          --error-bars=CI_accuracy
done










# Future Cleaning:
# rm -dr analysis/sim_traits/*/traits/avg_rank_all/n_markers_1
# rm -dr analysis/sim_traits/*/traits/avg_rank_all/n_markers_3