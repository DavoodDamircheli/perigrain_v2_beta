#!/bin/bash
path="`dirname \"$0\"`"
echo "path = $path"

logpath="/ddnA/work/debdeep/peri-wheel-output/logfiles"
name="chain-catch"
dateval=$(date +%Y-%m-%d-%H-%M-%S)

mkdir $logpath
balance

# qsub -o $logpath/$dateval-$name.log -N $name -v shape="$shape",meshsize_ratio="$meshsize_ratio",contact_rad_factor="$contact_rad_factor" $path/base_pbs.sh
qsub -o $logpath/$dateval-$name.log -N $name $path/base_pbs.sh
