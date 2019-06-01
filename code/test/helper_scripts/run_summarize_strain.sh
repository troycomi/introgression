#!/bin/bash

#SBATCH --time=0-1
#SBATCH -n 1
#SBATCH -o "/tigress/tcomi/aclark4_temp/results/filter_%A"

config=/home/tcomi/projects/aclark4_introgression/code/config.yaml

module load anaconda3
conda activate introgression3

introgression \
    --config $config \
    --log-file test.log \
    -vv \
    summarize-strains
