#!/bin/bash
path="`dirname \"$0\"`"
echo "path = $path"

# examples/3d-bulktest/run.sh -h $meshsize_ratio -c $contact_rad_factor

meshsize_ratio=3
contact_rad_factor=3
# shape="sphere"
shape='disk_w_hole_3d'

logpath="/ddnA/work/debdeep/peri-wheel-output/logfiles"
name="3dbulk-${meshsize_ratio}-${contact_rad_factor}"
dateval=$(date +%Y-%m-%d-%H-%M-%S)

mkdir $logpath
balance

qsub -o $logpath/$dateval-$name.log -N $name -v shape="$shape",meshsize_ratio="$meshsize_ratio",contact_rad_factor="$contact_rad_factor" $path/base_pbs.sh
