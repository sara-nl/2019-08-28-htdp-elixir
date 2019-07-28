#!/bin/bash
#SBATCH -t 10:00
#SBATCH -c 1
#SBATCH --constraint=skylake

echo "================================================================="
echo "                        Welcome to Spider"
echo "================================================================="
echo ""
echo "The status of the worker nodes is as displayed below:"
sinfo
echo ""
echo "The current running jobs are listed below:"
squeue
echo ""
echo "You just ran your first job on" $HOSTNAME " with a job ID " $SLURM_JOBID
echo ""
sleep 15s
