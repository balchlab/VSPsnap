jobname=`echo $1 | cut -f 1 -d '.'`

echo "#!/bin/csh"
echo "#SBATCH --job-name=${jobname}"
echo '#SBATCH --nodes=1'
echo '#SBATCH --ntasks=1'
#echo '#SBATCH --partition=highmem'
echo '#SBATCH --mem=150GB'
echo '#SBATCH --time=240:00:00'
#echo '#SBATCH --cpus-per-task=4'
echo "#SBATCH --output=${jobname}_err.log"
echo
echo 'cd $SLURM_SUBMIT_DIR'
echo 'module load python/3.8.3'
echo
echo "python $1"
echo
