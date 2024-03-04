#!/bin/bash
path="`dirname \"$0\"`"

#for seed in 1 2 3 4 5 6 7 8 
for seed in 1
do
    #$path/run.sh -L 100e-3 -l 50e-3 -P 50e-3
    $path/run.sh -P 170e-3
done
