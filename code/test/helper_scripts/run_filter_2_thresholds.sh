#!/bin/bash

#SBATCH --time=0-1
#SBATCH -n 1
#SBATCH -o "/tigress/tcomi/aclark4_temp/results/thresh_%A"

export PYTHONPATH=/home/tcomi/projects/aclark4_introgression/code/

module load anaconda3
conda activate introgression3

ARGS="_test .001 viterbi 10000 .025 10000 .025 10000 .025 10000 .025 unknown 1000 .01"

python ${PYTHONPATH}analyze/filter_2_thresholds_main.py $ARGS
