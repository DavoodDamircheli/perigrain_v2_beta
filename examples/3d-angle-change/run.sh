#!/bin/bash
path="`dirname \"$0\"`"
hpc="no"
#hpc="yes"

non_hpc_cores=4

frac="nofrac"
resume="no"





echo "#####################################################"

# paths
# dir=/home/davood/Downloads/kalthoff_different_angle/OUT_RESULT
#dir=/home/davood/Downloads/kalthoff_different_angle/one_image
#
dir=/home/davood/Downloads/plus3d_col
gen_plot(){

    echo 'gen_plot:Generating plots'
    dotsize=0.5
    bond_linewidth=0.1
    fcval=24

    # qq=force
    qq=damage
    # cmap='Greys_r'
    #cmap='Greys'
    cmap='viridis'
    #------------------kalthoff-------------------
    #python3 plot_kalthoff.py --all_dir $dir --setup_file $dir/setup.h5 --dotsize .5 --view_az 10 --view_angle 30  
    # 
    #python3 plot_kalthoff_filter.py --all_dir $dir --setup_file $dir/setup.h5 --dotsize 10.0 --view_az 20 --view_angle 210  
    #####python3 plot_kalthoff_filter.py --all_dir $dir --setup_file $dir/setup.h5 --dotsize 15.0 --view_az 11 --view_angle 145  
    #python3 plot_kalthoff_plotly.py --all_dir $dir --setup_file $dir/setup.h5 --dotsize 10.0 --view_az 5 --view_angle 30  
    #------------------plus -3d----------------------

    python3 plot_3d_paper_filter.py --all_dir $dir --setup_file $dir/setup.h5 --dotsize 10.0 --view_az 11 --view_angle 0  
    #sxiv $dir/*.png &
}
echo "#####################################################"



gen_vid(){
    echo 'Generating video'
    ./gen_vid.sh $dir
}
gen_plot
