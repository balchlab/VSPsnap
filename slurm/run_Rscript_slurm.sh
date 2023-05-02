jobname=`echo $1 | tail -c 14`

echo "#!/bin/csh"
echo "#SBATCH --job-name=${jobname}"
echo '#SBATCH --nodes=1'
echo '#SBATCH --ntasks=1'
#echo '#SBATCH --mem=16gb'
echo '#SBATCH --time=48:00:00'
#echo '#SBATCH --cpus-per-task=8'
echo "#SBATCH --output=${jobname}_err.log"
echo
echo 'cd $SLURM_SUBMIT_DIR'
echo 'module load R'
echo
echo "R CMD BATCH --no-restore --no-save $1 ${jobname}.Rout"
echo
