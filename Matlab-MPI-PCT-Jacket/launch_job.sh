#!/bin/bash
#PBS -l nodes=63:ppn=12

# Join the STDOUT and STDERR streams into STDOUT
# PBS -j oe

# Write STDOUT to output.txt
#PBS -N matlab-gpu-pblm

export WORK_DIR=/home13/jburkhar/software_projects/plbm/Matlab-MPI-PCT-Jacket
cd $WORK_DIR

# Initialize and clean Environment Modules 
. /usr/local/packages/Modules/current/init/bash
module purge
module load matlab

date; matlab -display null -nosplash < launch_job.m ; date;
