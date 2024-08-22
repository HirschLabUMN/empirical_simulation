# Set working directory
cd /home/hirschc1/burns756/empirical_sim

traits=()
for trait in $(ls -d analysis/gwas/*/); do
  #echo $trait
  traits+=$(basename $trait)
  traits+=" "
done

echo $traits

for trait in ${traits[@]}; do
    for model in additive dominant; do
        sbatch --export=trait=${trait},model=${model} scripts/all_non-sig_ld_calc.sh
    done
done