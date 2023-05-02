for sample in ${@:2}; do

./$1 ${sample} | sbatch

done

