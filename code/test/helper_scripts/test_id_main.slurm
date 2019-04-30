#!/bin/bash

## SBATCH --array=1
#SBATCH --time=2-0

#SBATCH -n 1
#SBATCH -o "/tigress/tcomi/aclark4_temp/results/id_%A"

config=/home/tcomi/projects/aclark4_introgression/code/config.yaml

module load anaconda3
conda activate introgression3

introgression \
    --config $config \
    -vvvv \
    --log-file test.log \
    id-regions
