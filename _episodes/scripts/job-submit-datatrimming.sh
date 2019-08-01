#!/bin/bash
#SBATCH -c 1
#SBATCH --constraint=skylake
#SBATCH --reservation=spidercourse_2

#As data manager uncomment the following line
#bash /project/spidercourse/Data/ecoli-analysis/data_qc.sh 

#As a regular user uncomment the following line
#bash $HOME/ecoli-analysis/data_qc.sh 
