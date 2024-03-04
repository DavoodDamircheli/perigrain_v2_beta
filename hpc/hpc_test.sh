#!/bin/bash
## calling the pbs script with specified logfile, jobname, and variables
## shape, theta, meshsize

nodes=$1
shape=n4
meshsize=80e-3
for theta in 1
do
    qsub -o hpc_test_v${nodes}_${shape}_${theta}_${meshsize}.log -N ${shape}_${theta}_${meshsize} -v shape="$shape",theta="$theta",meshsize="$meshsize",nodes="$nodes" base_test_n$nodes.sh
done
