# peri-wheel location
cd $HOME/peri-wheel

## output plot location
#plt_dir=/bigpart/plots
## local
plt_dir=$HOME/local_wheelplots

mkdir -p $plt_dir

fc=1
#lc=600
lc=460
modulo=5

#avg=""
avg="--avg"
avg_window=5

#meshsize=80e-3
meshsize=125e-3

quantities=(
disp_x
vel_x
vel_y
acc_x
currpos_x
slip
#cumulative_torque_density
)

titles=(
"Horizontal displacement"
"Horizontal velocity" 
"Vertical velocity" 
"Horizontal acceleration" 
"Horizontal position"
"Slip"
#"Total torque density applied on the wheel to maintain angular velocity"
)

for quantity in "${!quantities[@]}"
#for quantity in "slip"
do

    ####
    # wheel size to gravel size ratio
    parent_dir=examples_output/fixed_angular

    shape="n4"
    prefix=sizeratio_${shape}

    dirlist=""
    labellist=""
    for meshsize in 160e-3 175e-3 190e-3
    do
	dirname=nofrac_${shape}_1_${meshsize}_1
	dirlist="${dirlist} ${dirname}"
	labellist="${labellist} ${meshsize}"
    done
    #labellist="36mm 42mm 43mm"
    # python3 wheel_attrib_multi.py -p $parent_dir --list $dirlist --labels $labellist --img_dir $plt_dir --prefix $prefix --quantity "${quantities[$quantity]}" --title "${titles[$quantity]}" --fc $fc --lc $lc --modulo $modulo $avg --avg_window $avg_window
    python3 wheel_attrib_multi.py -p $parent_dir --list $dirlist --labels $labellist --img_dir $plt_dir --prefix $prefix --quantity "${quantities[$quantity]}" --title "${titles[$quantity]}" --fc $fc --lc $lc --modulo $modulo

    prefix="size_ratio_timeavg_${shape}"
    python3 wheel_attrib_multi.py -p $parent_dir --list $dirlist --labels $labellist --img_dir $plt_dir --prefix $prefix --quantity "${quantities[$quantity]}" --title "${titles[$quantity]}" --fc $fc --lc $lc --modulo $modulo --plottype timeavg

    #prefix="n4"
    #shape="plus"
    #meshsize=80e-3
    #prefix=${shape}_${meshsize}
    #python3 wheel_attrib_multi.py -p /bigpart/peri-wheel-output/fixed_angular --list nofrac_${shape}_0.4_${meshsize} nofrac_${shape}_0.6_${meshsize} nofrac_${shape}_0.8_${meshsize} nofrac_${shape}_1_${meshsize} --labels "0.4" "0.6" "0.8" "1" --img_dir $plt_dir --prefix $prefix --quantity "${quantities[$quantity]}" --title "${titles[$quantity]}" --fc $fc --lc $lc --modulo $modulo $avg --avg_window $avg_window

    #prefix='shapes'
    #meshsize=80e-3
    #python3 wheel_attrib_multi.py -p /bigpart/peri-wheel-output/fixed_angular --list nofrac_n4_0.4_$meshsize nofrac_n4_1_$meshsize nofrac_ring_0.8_$meshsize nofrac_L_0.6_$meshsize nofrac_plus_0.4_$meshsize --labels rectangle square ring L plus --img_dir $plt_dir --prefix $prefix --quantity "${quantities[$quantity]}" --title "${titles[$quantity]}" --fc $fc --lc $lc --modulo $modulo $avg --avg_window $avg_window

    #prefix='shapes_125'
    #python3 wheel_attrib_multi.py -p /bigpart/peri-wheel-output/fixed_angular --list nofrac_n4_1_125e-3 nofrac_pertdisk_inc_1_125e-3 nofrac_pertdisk_inc_0.6_125e-3 nofrac_plus_abs_0.4_125e-3 --labels square circle ring plus --img_dir $plt_dir --prefix $prefix --quantity "${quantities[$quantity]}" --title "${titles[$quantity]}" --fc $fc --lc $lc --modulo $modulo $avg --avg_window $avg_window


    #prefix='wheel_size'
    #python3 wheel_attrib_multi.py -p "examples_output" --list "gravel_5_angular/n4_nofrac_grainsize_100" "gravel_5_angular_wheel_2.5/n4_nofrac_grainsize_100" "gravel_5_angular_wheel_2/n4_nofrac_grainsize_100" --labels "h/3" "h/2.5" "h/2" --img_dir $plt_dir --prefix $prefix --quantity "${quantities[$quantity]}" --title "${titles[$quantity]}" --fc $fc --lc $lc --modulo $modulo


done
