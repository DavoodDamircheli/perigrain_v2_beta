#!/bin/bash
path="`dirname \"$0\"`"

# parent_dir=/work/$USER/peri-wheel-output/wheel

## older fixed_angular data
parent_dir=/ddnA/work/debdeep/pier-wheel-output-old


### frac
### only shape left: circ (in progress on smic)
#meshsize=80e-3
#theta=1
#for subdir in frac_n4_1_80e-3 frac_plus_0.4_80e-3_0.5 frac_ring_0.6_80e-3_0.3
#do
#    dir="${parent_dir}/${subdir}"
#    python3 $HOME/peri-wheel/examples/wheel/extract.py --data_dir $dir
#done



## nofrac
#prefix='shapes'
meshsize=80e-3
# for subdir in nofrac_n4_0.4_$meshsize nofrac_n4_1_$meshsize nofrac_ring_0.8_$meshsize nofrac_plus_0.4_$meshsize nofrac_circ_1_$meshsize
#do
#    dir="${parent_dir}/${subdir}"
#    python3 $HOME/peri-wheel/examples/wheel/extract.py --data_dir $dir
#done

##prefix='shapes_125'
### We don't have rings, or circles
#for subdir in nofrac_n4_1_125e-3 nofrac_pertdisk_inc_1_125e-3 nofrac_pertdisk_inc_0.6_125e-3 nofrac_plus_abs_0.4_125e-3 
#do
#    dir="${parent_dir}/${subdir}"
#    python3 $HOME/peri-wheel/examples/wheel/extract.py --data_dir $dir
#done
