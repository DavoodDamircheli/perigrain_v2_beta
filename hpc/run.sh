#!/bin/bash
## calling the pbs script with specified logfile, jobname, and variables
## shape, theta, meshsize

#shape=ring
#meshsize=80e-3
#for theta in 0.4 0.6
#do
    #qsub -o ${shape}_${theta}_${meshsize}.log -N ${shape}_${theta}_${meshsize} -v shape="$shape",theta="$theta",meshsize="$meshsize" base_pbs.sh
#done

#shape=L
#meshsize=80e-3
#for theta in 0.4 0.6 0.8
#do
    #qsub -o ${shape}_${theta}_${meshsize}.log -N ${shape}_${theta}_${meshsize} -v shape="$shape",theta="$theta",meshsize="$meshsize" base_pbs.sh
#done

#shape="circ"
#meshsize=80e-3
#for theta in 1
#do
    #qsub -o ${shape}_${theta}_${meshsize}.log -N ${shape}_${theta}_${meshsize} -v shape="$shape",theta="$theta",meshsize="$meshsize" base_pbs.sh
#done

#shape="pertdisk_inc"
#meshsize=125e-3
#for theta in 0.4
#do
    #qsub -o ${shape}_${theta}_${meshsize}.log -N ${shape}_${theta}_${meshsize} -v shape="$shape",theta="$theta",meshsize="$meshsize" base_pbs.sh
#done


#shape="plus_abs"	# the parameter of plus is controlled substituted by theta
#meshsize=125e-3
##for theta in 0.4 0.6 0.8 1
#for theta in 0.4 0.6 0.8 1
#do
    #qsub -o ${shape}_${theta}_${meshsize}.log -N ${shape}_${theta}_${meshsize} -v shape="$shape",theta="$theta",meshsize="$meshsize" base_pbs.sh
#done

#######################################################################
# frac

#frac="frac"

#shape=ring
#meshsize=80e-3
#for theta in 0.8
#do
    #qsub -o ${shape}_${theta}_${meshsize}.log -N ${shape}_${theta}_${meshsize} -v shape="$shape",theta="$theta",meshsize="$meshsize",frac="$frac" base_pbs.sh
#done

#shape=plus
#meshsize=80e-3
#for theta in 0.4
#do
    #qsub -o ${shape}_${theta}_${meshsize}.log -N ${shape}_${theta}_${meshsize} -v shape="$shape",theta="$theta",meshsize="$meshsize",frac="$frac" base_pbs.sh
#done

#shape=n4
#meshsize=80e-3
#for theta in 0.4 1
#do
    #qsub -o ${frac}_${shape}_${theta}_${meshsize}.log -N ${frac}_${shape}_${theta}_${meshsize} -v shape="$shape",theta="$theta",meshsize="$meshsize",frac="$frac" base_pbs.sh
#done


#######################################################################

frac="frac"
#finermsz=0.3
meshsize=80e-3
# shape=plus
shape=ring
# shape=circ
# theta=0.4
theta=0.6
# theta=1
# for finermsz in 0.5 0.3 0.25
#for finermsz in 0.3
for finermsz in 1.2
do
    qsub -o ${frac}_${shape}_${theta}_${meshsize}_${finermsz}.log -N ${frac}_${shape}_${theta}_${meshsize}_${finermsz} -v shape="$shape",theta="$theta",meshsize="$meshsize",frac="$frac",finermsz="$finermsz" base_pbs.sh
done


