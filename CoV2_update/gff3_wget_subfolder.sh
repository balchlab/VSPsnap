echo "#!/bin/csh"
echo "#SBATCH --job-name=$1"
echo '#SBATCH --nodes=1'
echo '#SBATCH --ntasks=1'
#echo '#SBATCH --mem=16gb'
echo '#SBATCH --time=120:00:00'
#echo '#SBATCH --cpus-per-task=8'
echo "#SBATCH --output=$1_wget_err.log"
echo
echo 'cd $SLURM_SUBMIT_DIR'

echo 'mkdir $1'
echo 'cd $1'

echo "wget -r -l 1 -nH -nd ftp://download.big.ac.cn/GVM/Coronavirus/gff3/$1/"
echo "wget --no-remove-listing ftp://download.big.ac.cn/GVM/Coronavirus/gff3/$1/"