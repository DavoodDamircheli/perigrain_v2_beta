#!/bin/bash
path="`dirname \"$0\"`"

# peri-wheel location
cd $HOME/peri-wheel

## dir where all subdirs reside
parent_dir=/ddnA/work/debdeep/pier-wheel-output-old
## output plot location
plt_dir=$parent_dir/plots
mkdir -p $plt_dir

quantities=(
disp_x
#disp_y
vel_x
#vel_y
#vel_angular
#acc_x
#acc_y
#currpos_x
#currpos_y
slip
#applied_fd_x
#applied_fd_y
#applied_torque_density
#applied_f_x
#applied_f_y
#applied_torque
)

##for timeavg in "" "--timeavg" "--cumulative" #"--regression_only"
#for timeavg in ""
#do
#    for qindex in "${!quantities[@]}"
#    do
#	quantity="${quantities[$qindex]}"
#	fc=1
#	lc=550
#	modulo=5
#
#	prefix="fixed_angular_${fc}_${timeavg}"
#
#	meshsize=80e-3
#	subdirs="nofrac_circ_1_$meshsize nofrac_ring_0.8_$meshsize nofrac_plus_0.4_$meshsize nofrac_n4_1_$meshsize"
#
#	labels="Disk Ring Plus Square"
#	timeavg_xlabel="--timeavg_xlabel Slip"
#
#	python3 $HOME/peri-wheel/examples/wheel/multiplot.py -p $parent_dir --subdirs $subdirs --labels $labels --quantity $quantity --prefix $prefix --img_dir $plt_dir $timeavg --fc $fc --lc $lc --timeavg_drawline $timeavg_xlabel --xquantity time&
#    done
#done



### shapes, with fracture

for timeavg in ""
do
    for qindex in "${!quantities[@]}"
    do
	quantity="${quantities[$qindex]}"
	fc=1
	lc=550
	modulo=5

	prefix="frac_fixed_angular_${fc}_${timeavg}"

	meshsize=80e-3
	subdirs=" frac_ring_0.6_80e-3_0.3 frac_plus_0.4_80e-3_0.5 frac_n4_1_80e-3"
	labels="Disk Ring Plus Square"
	timeavg_xlabel="--timeavg_xlabel Slip"

	python3 $HOME/peri-wheel/examples/wheel/multiplot.py -p $parent_dir --subdirs $subdirs --labels $labels --quantity $quantity --prefix $prefix --img_dir $plt_dir $timeavg --fc $fc --lc $lc --timeavg_drawline $timeavg_xlabel --xquantity time&
    done
done
