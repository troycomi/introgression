#!/bin/bash

#SBATCH --time=0-1
#SBATCH -n 1
#SBATCH -o "/tigress/tcomi/aclark4_temp/results/filter_%A"

export PYTHONPATH=/home/tcomi/projects/aclark4_introgression/code/

module load anaconda
conda activate introgression

ARGS="p4e2 .001 viterbi 10000 .025 10000 .025 10000 .025 10000 .025 unknown 1000 .01"

python ${PYTHONPATH}analyze/summarize_strain_states_main.py $ARGS
