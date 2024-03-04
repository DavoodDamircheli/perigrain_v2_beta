#!/bin/bash
# No shell commands before PBS is set up.
#
# "workq" is the default job queue.
#PBS -q workq 
#
# Set the appropriate project allocation code
#PBS -A hpc_periwh01 
#
# Set number of nodes and number of processors on each node
# to be used. See cluster user guide for corresponding ppn number
#PBS -l nodes=4:ppn=20
#
# Set time job is allowed to run in hh:mm:ss
#PBS -l walltime=00:30:00 
#
# Send stdout to a named file: move to external command line
#
# Merge stderr messages with stdout
#PBS -j oe 
#
# Give the job a name for easier tracking
#PBS -N wheeljob-$shape_$theta_$meshsize
#
# Shell commands may begin here

#export PROJ_DIR=/work/$USER/peri-wheel
export PROJ_DIR=$HOME/peri-wheel
cd $PROJ_DIR

examples/hpc_test/run.sh $shape $theta $meshsize $nodes
