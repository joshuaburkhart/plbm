#!/bin/bash
#PBS -l nodes=cn27:ppn=3+cn28:ppn=3+cn29:ppn=3+cn30:ppn=3+cn31:ppn=3+cn32:ppn=3+cn33:ppn=3+cn34:ppn=3+cn35:ppn=3+cn36:ppn=3+cn61:ppn=3+cn62:ppn=3+cn63:ppn=3+cn64:ppn=3+cn65:ppn=3+cn66:ppn=3+cn67:ppn=3+cn68:ppn=3+cn69:ppn=3+cn70:ppn=3+cn71:ppn=3+cn72:ppn=3+cn97:ppn=3+cn98:ppn=3+cn99:ppn=3+cn100:ppn=3+cn101:ppn=3+cn102:ppn=3+cn103:ppn=3+cn104:ppn=3+cn105:ppn=3+cn106:ppn=3+cn107:ppn=3


# Join the STDOUT and STDERR streams into STDOUT
# PBS -j oe

# Write STDOUT to output.txt
#PBS -N matlab-gpu-pblm
#PBS -l walltime=80:0:00

export WORK_DIR=/home11/ozog/School/CIS-555-Computational_Science/CIS555-Project/LargerDataSet-GPU-pfor-MPI
cd $WORK_DIR

# Initialize and clean Environment Modules 
. /usr/local/packages/Modules/current/init/bash
module purge
module load matlab

date; matlab -display null -nosplash < launch_job.m ; date;
