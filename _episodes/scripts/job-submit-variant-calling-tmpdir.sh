#!/bin/bash
#SBATCH -c 1
#SBATCH --constraint=skylake

mkdir "$TMPDIR"/var-calling
cd "$TMPDIR"/var-calling

cp $HOME/ecoli-analysis/run-variant-calling-tmpdir.sh .

time bash run-variant-calling-tmpdir.sh 
