#!/bin/bash
#SBATCH -c 1
#SBATCH --constraint=skylake

bash /project/spidercourse/Data/ecoli-analysis/run-variant-calling.sh 
