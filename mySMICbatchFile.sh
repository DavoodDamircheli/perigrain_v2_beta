#!/bin/bash
#PBS -q workq
#PBS -A ?????
#PBS -l nodes=1:ppn4
#PBS -l walltime=00:05:00

cd $PBS_O_WORKDIR
mpirun ./examples/3d-coll-multi_plus/run.sh
