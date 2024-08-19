#!/bin/bash
path="`dirname \"$0\"`"
hpc="no"
#hpc="yes"

non_hpc_cores=4

frac="nofrac"
resume="no"





echo "#####################################################"

# paths
#dir=/home/davood/Downloads/3d-grain/N1000
dir=/home/davood/Downloads/3d-grain/N40T64/N40T64_8

#dir=/home/davood/Downloads/3d-grain/N30C2F1T8/N2 
gen_plot(){

    echo 'gen_plot:Generating plots'
    dotsize=0.5
    bond_linewidth=0.1
    fcval=24

    # qq=force
    qq=damage
    cmap='viridis'
    # 
    #python3 plot_potato_2.py --all_dir $dir --setup_file $dir/setup.h5 --dotsize 30   
    #python3 plot3d_visual.py --all_dir $dir --setup_file $dir/setup.h5 --dotsize 30   
    python3  plot3d_boxhighlight.py --all_dir $dir --setup_file $dir/setup.h5 --dotsize 30
   }
echo "#####################################################"



gen_vid(){
    echo 'Generating video'
    ./gen_vid.sh $dir
}
gen_plot
