#!/bin/bash

#SBATCH --time=0-4
#SBATCH --array=0,2-5
#SBATCH -n 1
#SBATCH -o "/tigress/tcomi/aclark4_temp/results/summarize_%A_%a"

module load anaconda3
conda activate introgression3

config=/home/tcomi/projects/aclark4_introgression/code/config.yaml

introgression \
    --config $config \
    --log-file test.log \
    -vvvv \
    summarize-regions \
        --state N_45 \
        --state CBS432
