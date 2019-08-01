#!/bin/bash
#SBATCH -c 1
#SBATCH --constraint=skylake
#SBATCH --reservation=spidercourse_2

bash /project/spidercourse/Data/ecoli-analysis/run-variant-calling.sh 
