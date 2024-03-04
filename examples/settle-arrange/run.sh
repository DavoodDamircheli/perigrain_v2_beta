#!/bin/bash
path="`dirname \"$0\"`"

hpc="no"
#hpc="yes"
non_hpc_cores=4


# input from command line arg
arglist=""
while getopts h:r:s:z:m:M:R:K:G:p:L:l:P: flag
do
    case "${flag}" in
        h) 
	    shape="--shape ${OPTARG}"
	    arglist="${arglist}_${flag}-${OPTARG}";;
        s) 
	    seed="--seed ${OPTARG}"
	    arglist="${arglist}_${flag}-${OPTARG}";;
        g) 
	    meshsize="--grain_meshsize ${OPTARG}"
	    arglist="${arglist}_${flag}-${OPTARG}";;
        P) 
	    P_meshsize="--P_meshsize ${OPTARG}"
	    arglist="${arglist}_${flag}-${OPTARG}";;
        m) 
	    minrad="--minrad ${OPTARG}"
	    arglist="${arglist}_${flag}-${OPTARG}";;
        L) 
	    LL="--L ${OPTARG}"
	    arglist="${arglist}_${flag}-${OPTARG}";;
        l) 
	    hh="--hh ${OPTARG}"
	    arglist="${arglist}_${flag}-${OPTARG}";;
        M) 
	    maxrad="--maxrad ${OPTARG}"
	    arglist="${arglist}_${flag}-${OPTARG}";;
        R) 
	    rho_scale="--rho_scale ${OPTARG}"
	    arglist="${arglist}_${flag}-${OPTARG}";;
        K) 
	    K_scale="--K_scale ${OPTARG}"
	    arglist="${arglist}_${flag}-${OPTARG}";;
        G) 
	    G_scale="--G_scale ${OPTARG}"
	    arglist="${arglist}_${flag}-${OPTARG}";;
        p) 
	    plot="--plot";;
        n) 
	    non_hpc_cores=${OPTARG};;
    esac
done


#resume="yes"
resume="no"

dotsize=1

## specify fracture via 4th argument
frac="nofrac"
#frac="frac"

dirname="settle_arrange"

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
setuppy=$dir/setup.py
setuppng=$dir/setup.png
echo "logfile: $logfile"

create_env(){
    mkdir -p $dir
    #clear logfile
    echo '' > $logfile
    echo "arglist: $arglist" >> $logfile

    cp $path/base.conf $config
    cp $path/setup.py  $setuppy

    # enable fracture of not at a given timestep
    if [ "$frac" = "frac" ]
    then
	echo "enable_fracture = 0" >> $config
	#echo "set_fracture_timestep = $wheel_pos_timestep" >> $config

	echo "self_contact = 0" >> $config
	#echo "set_self_contact_timestep = $wheel_pos_timestep" >> $config
    else
	echo "enable_fracture = 0" >> $config
	echo "self_contact = 0" >> $config
    fi

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
    echo "Generating experiment setup"

    ##setupcmd="python3 $path/setup.py --setup_file $sfile --pngfile $setuppng $shape $seed $grain_meshsize $minrad $maxrad $rho_scale $K_scale $G_scale $plot"

    setupcmd="python3 $path/setup.py --setup_file $sfile --pngfile $setuppng $shape $seed ${grain_meshsize} ${P_meshsize} $minrad $maxrad $rho_scale $K_scale $G_scale $LL $hh $plot"
    echo "$setupcmd"
    echo "$setupcmd" >> $logfile
    eval `$setupcmd >> $logfile`
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
	mpirun -machinefile $PBS_NODEFILE -np $NPROCS bin/simulate2d -c $config -o $dir -i $sfile  >> $logfile
    else
	echo "On non-hpc with ${non_hpc_cores} cores"
	mpirun -n ${non_hpc_cores} bin/simulate2d -c $config -o $dir -i $sfile  >> $logfile
    fi
}

gen_plot(){
    echo 'Generating plots'
    python3 plot_timestep.py --data_dir $dir --img_dir $dir --setup_file $sfile --dotsize $dotsize --quantity force --nocolorbar --colormap 'Greys_r'
    echo "$dir"
    sxiv $dir/*.png &
}

#gen_vid(){
    #echo 'Generating video'
    #./gen_vid.sh $dir
#}

create_env
gen_setup
run
gen_plot
#gen_vid
