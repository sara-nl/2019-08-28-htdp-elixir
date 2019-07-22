#!/bin/bash
#SBATCH -c 1
#SBATCH --constraint=skylake

bash /project/spidercourse/Data/ecoli-analysis/run_variant_calling.sh 
