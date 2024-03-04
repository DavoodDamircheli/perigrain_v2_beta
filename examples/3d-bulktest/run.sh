#!/bin/bash
path="`dirname \"$0\"`"
echo "path = $path"

hpc="no"
# hpc="yes"
non_hpc_cores=4

frac="nofrac"

# input from command line arg
arglist=""
# while getopts af:ph:r:L:x:h:y:t:R:G:K:n:v: flag
while getopts a:c:f:G:h:HK:L:pr:R:s:t:m:n:v::x:y: flag
do
    case "${flag}" in
        a) 
	    acc="--acc_val ${OPTARG}"
	    arglist="${arglist}_${flag}-${OPTARG}";;
        c) 
	    contact_rad_factor="--contact_rad_factor ${OPTARG}"
	    arglist="${arglist}_${flag}-${OPTARG}";;
	f)
	    frac="frac";;
        G) 
	    G_scale="--G_scale ${OPTARG}"
	    arglist="${arglist}_${flag}-${OPTARG}";;
        h) 
	    mshfac="--meshsize_factor ${OPTARG}"
	    arglist="${arglist}_${flag}-${OPTARG}";;
        H) 
	    hpc="yes"
	    arglist="${arglist}_${flag}-${OPTARG}";;
        K) 
	    K_scale="--K_scale ${OPTARG}"
	    arglist="${arglist}_${flag}-${OPTARG}";;
        L) 
	    L="--L ${OPTARG}"
	    arglist="${arglist}_${flag}-${OPTARG}";;
        m) 
	    delta_factor="--delta_factor ${OPTARG}"
	    arglist="${arglist}_${flag}-${OPTARG}";;
        n) 
	    non_hpc_cores=${OPTARG};;
        p) 
	    plot="--plot";;
        r) 
	    particle_rad="--particle_rad ${OPTARG}"
	    arglist="${arglist}_${flag}-${OPTARG}";;
        R) 
	    rho_scale="--rho_scale ${OPTARG}"
	    arglist="${arglist}_${flag}-${OPTARG}";;
        s) 
	    shape="--shape ${OPTARG}"
	    arglist="${arglist}_${flag}-${OPTARG}";;
        t) 
	    wallh_ratio="--wallh_ratio ${OPTARG}"
	    arglist="${arglist}_${flag}-${OPTARG}";;
        v) 
	    vel="--vel_val ${OPTARG}"
	    arglist="${arglist}_${flag}-${OPTARG}";;
        x) 
	    nx="--nx ${OPTARG}"
	    arglist="${arglist}_${flag}-${OPTARG}";;
        y) 
	    ny="--ny ${OPTARG}"
	    arglist="${arglist}_${flag}-${OPTARG}";;
    esac
done

echo "hpc = $hpc"

#resume="yes"
resume="no"

gen_setup="yes"
#gen_setup="no"


# dirname="3dtest"
pwd="$(cd `dirname $0` && pwd)"
dirname="$(basename $pwd)"
# echo "saving to dirname $dirname"
# exit

if [ "$hpc" = "yes" ]
then
    data_output_loc="/work/$USER/peri-wheel-output/$dirname"
else
    data_output_loc="examples_output/$dirname"
fi

# create subdirectory
dir=$data_output_loc/${frac}_${arglist}

config=$dir/main.conf
sfile=$dir/setup.h5
logfile=$dir/output.log
echo "logfile: $logfile"

create_env(){
    mkdir -p $dir
    #clear logfile
    echo '' > $logfile
    echo $(cd $path; git show -s) > $logfile
    echo "arglist: $arglist" >> $logfile

    cp $path/base.conf $config
    cp $path/setup.py  $dir/

    # enable fracture of not at a given timestep
    # if [ "$frac" = "frac" ]
    # then
	# echo "enable_fracture = 0" >> $config
	# echo "set_fracture_timestep = $wheel_pos_timestep" >> $config
    #
	# echo "self_contact = 0" >> $config
	# echo "set_self_contact_timestep = $wheel_pos_timestep" >> $config
    # else
	# echo "enable_fracture = 0" >> $config
	# echo "self_contact = 0" >> $config
    # fi
    #
    if [ "$resume" = "yes" ]
    then
	## get the last index
	last=$(ls $dir/tc_*.h5 | tail -1) # Get the largest indices
	last=${last##*/} # strip the path
	last="${last%.*}" # strip the extension
	last=$(echo $last | awk -F '_' '{ print $2 }')

	echo "do_resume = 1" >> $config
	echo "resume_ind = $last" >> $config
	echo "wall_resume = 1" >> $config
    else
	echo ''
    fi
}

gen_setup(){
    echo "Meshing boundary box interior"
    if [ "$hpc" = "yes" ]
    then
	gmshpath=$HOME/.local/bin/gmsh
    else
	gmshpath=gmsh
    fi

    $gmshpath $path/../../meshdata/3d/closepacked_cube.geo -3
    echo "Generating experiment setup"
    # echo "python3 $path/setup.py --setup_file $sfile $hopper_ratio $particle_rad $nx $ny $L $wallh_ratio $vel $acc $rho_scale $K_scale $G_scale $mshfac $plot >> $logfile"
    echo "python3 $path/setup.py --setup_file $sfile $hopper_ratio $particle_rad $nx $ny $L $wallh_ratio $vel $acc $rho_scale $K_scale $G_scale $mshfac $delta_factor $contact_rad_factor $plot $shape" >> $logfile
    python3 $path/setup.py --setup_file $sfile $hopper_ratio $particle_rad $nx $ny $L $wallh_ratio $vel $acc $rho_scale $K_scale $G_scale $mshfac $delta_factor $contact_rad_factor $plot $shape >> $logfile
    # cp setup.png $dir/
}


# run code
run()
{
    echo 'running'

    if [ "$hpc" = "yes" ]
    then
	echo 'On hpc'
	export NPROCS=`wc -l $PBS_NODEFILE |gawk '//{print $1}'`
	echo "NPROCS=" 
	echo $NPROCS
	mpirun -machinefile $PBS_NODEFILE -np $NPROCS bin/simulate3d -c $config -o $dir -i $sfile  >> $logfile
    else
	echo "On non-hpc with ${non_hpc_cores} cores"
	mpirun -n ${non_hpc_cores} bin/simulate3d -c $config -o $dir -i $sfile  >> $logfile
    fi
}

gen_plot(){
    echo 'Generating plots'

    dotsize=0.5
    alpha=0.4
    bond_linewidth=0.1
    fcval=24

    # qq=force
    qq=damage
    # cmap='Greys_r'
    cmap='Greys'
    # cmap='viridis'
    # python3 plot3d_timestep.py --data_dir $dir --img_dir $dir --setup_file $sfile --dotsize $dotsize --quantity $qq --nocolorbar --colormap $cmap #--plot_bonds --bond_linewidth ${bond_linewidth} --fc $fcval
    python3 plot3d_timestep.py --all_dir $dir --alpha $alpha
    echo "$dir"
    sxiv $dir/*.png &
}

extract(){
    echo 'extracting data'
    python3 $path/extract.py --data_dir $dir

    # two directories to keep quantities.h5 data for each particle
    newdir=$dir/P_A
    mkdir $newdir
    python3 $path/extract.py --data_dir $dir --wheel_ind 0 --output_file $newdir/quantities.h5

    newdir=$dir/P_B
    mkdir $newdir
    python3 $path/extract.py --data_dir $dir --wheel_ind 1 --output_file $newdir/quantities.h5
}

multiplot(){
    echo 'plotting data'
    # quantity=vel_y
    quantity=vel_z

    for quantity in vel_x vel_y vel_z
    do
	python3 $path/multiplot.py -p $dir --subdirs P_A P_B --quantity $quantity --img_dir $dir
    done
}

#gen_vid(){
    #echo 'Generating video'
    #./gen_vid.sh $dir
#}

create_env
gen_setup
run
gen_plot
# extract
# multiplot
#
# gen_vid
