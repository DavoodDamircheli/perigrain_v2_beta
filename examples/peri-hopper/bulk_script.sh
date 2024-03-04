#!/bin/bash
path="`dirname \"$0\"`"

# for seed in 1 2 3 4 5 6 7 8 
# do
#     $path/run.sh -r 6e-3 -s $seed &
# done

# with fracture: -f to enable fracture; -g to specify fracture toughness scaling
$path/run.sh -r 6e-3 -s 1 -f -g 0.001

