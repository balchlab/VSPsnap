echo "#!/bin/csh"
echo "#SBATCH --job-name=$1"
echo '#SBATCH --nodes=1'
echo '#SBATCH --ntasks=1'
#echo '#SBATCH --mem=16gb'
echo '#SBATCH --time=120:00:00'
#echo '#SBATCH --cpus-per-task=8'
echo "#SBATCH --output=$1_err.log"
echo
echo 'cd $SLURM_SUBMIT_DIR'

echo "zip -r $1.zip $1"