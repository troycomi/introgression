#!/bin/bash

#SBATCH --time=6-0
#SBATCH -n 1
#SBATCH -o "/tigress/tcomi/aclark4_temp/results/predict_%A"

config=/home/tcomi/projects/aclark4_introgression/code/config.yaml
#Make sure chrms is set to only I

module load anaconda3
conda activate introgression3

introgression \
    --config $config \
    --log-file test.log \
    -vv \
    predict
