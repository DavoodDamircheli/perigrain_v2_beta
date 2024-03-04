#!/bin/bash
# No shell commands before PBS is set up.
#
# "workq" is the default job queue.
#PBS -q workq 
#
# Set the appropriate project allocation code
# PBS -A hpc_periwh01 
#PBS -A hpc_perigrain1
#
# Set number of nodes and number of processors on each node
# to be used. See cluster user guide for corresponding ppn number
#PBS -l nodes=4:ppn=20
# PBS -l nodes=1:ppn=20
# PBS -l nodes=4:ppn=20
# PBS -l nodes=5:ppn=20
#
# Set time job is allowed to run in hh:mm:ss
# PBS -l walltime=20:00:00 
#PBS -l walltime=21:00:00 
# PBS -l walltime=7:00:00 
# PBS -l walltime=0:20:00 
#
# Send stdout to a named file: move to external command line
#
# Merge stderr messages with stdout
#PBS -j oe 
#
# Give the job a name for easier tracking
#PBS -N 3dbulk-$meshsize_ratio-$contact_rad_factor
#
# Shell commands may begin here

#export PROJ_DIR=/work/$USER/peri-wheel
export PROJ_DIR=$HOME/peri-wheel
cd $PROJ_DIR

# -H for hpc=yes
examples/3d-bulktest/run.sh -H -h $meshsize_ratio -c $contact_rad_factor -s $shape
