#!/bin/bash
path="`dirname \"$0\"`"

hpc="no"
#hpc="yes"
non_hpc_cores=4

### default values of dt and timesteps
#dtval=1e-6
#tsval=400000

## with fracture: default values of dt and timesteps
dtval=5e-7
tsval=800000

## specify fracture via 4th argument
frac="nofrac"
#frac="frac"


# input from command line arg
arglist=""
while getopts fg:G:h:K:L:n:pr:R:s:t:T:w:x:y:z: flag
do
    case "${flag}" in
        f) 
	    frac="frac"
	    arglist="${arglist}_${flag}";;
        h) 
	    hopper_ratio="--hopper_ratio ${OPTARG}"
	    arglist="${arglist}_${flag}-${OPTARG}";;
        r) 
	    particle_rad="--particle_rad ${OPTARG}"
	    arglist="${arglist}_${flag}-${OPTARG}";;
        s) 
	    seed="--seed ${OPTARG}"
	    arglist="${arglist}_${flag}-${OPTARG}";;
        t) 
	    dtval="${OPTARG}"
	    arglist="${arglist}_${flag}-${OPTARG}";;
        T) 
	    tsval="${OPTARG}"
	    arglist="${arglist}_${flag}-${OPTARG}";;
        L) 
	    L="--L ${OPTARG}"
	    arglist="${arglist}_${flag}-${OPTARG}";;
        x) 
	    nx="--nx ${OPTARG}"
	    arglist="${arglist}_${flag}-${OPTARG}";;
        y) 
	    ny="--ny ${OPTARG}"
	    arglist="${arglist}_${flag}-${OPTARG}";;
        z) 
	    partmeshfac="--particle_meshsize_factor ${OPTARG}"
	    arglist="${arglist}_${flag}-${OPTARG}";;
        w) 
	    wallmeshfac="--wall_meshsize_factor ${OPTARG}"
	    arglist="${arglist}_${flag}-${OPTARG}";;
        R) 
	    rho_scale="--rho_scale ${OPTARG}"
	    arglist="${arglist}_${flag}-${OPTARG}";;
        K) 
	    K_scale="--K_scale ${OPTARG}"
	    arglist="${arglist}_${flag}-${OPTARG}";;
        g) 
	    Gnot_scale="--Gnot_scale ${OPTARG}"
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

gen_setup="yes"
#gen_setup="no"

dotsize=1
#dotsize=0.5

dirname="hopper"

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
setuppng=$dir/setup.png
echo "logfile: $logfile"

create_env(){
    mkdir -p $dir
    #clear logfile
    echo '' > $logfile
    echo "arglist: $arglist" >> $logfile

    cp $path/base.conf $config
    cp $path/setup.py  $dir/


    # dt and timestep specified here instead of the config file
    echo "timesteps = $tsval" >> $config
    echo "dt = $dtval " >> $config

    # enable fracture of not at a given timestep
    if [ "$frac" = "frac" ]
    then
	echo "enable_fracture = 1" >> $config
	echo "self_contact = 1" >> $config
	echo "nl_bdry_only = 0" >> $config
    else
	echo "enable_fracture = 0" >> $config
	echo "self_contact = 0" >> $config
	echo "nl_bdry_only = 1" >> $config
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
    echo "python3 $path/setup.py --setup_file $sfile $hopper_ratio $particle_rad $nx $ny $L $rho_scale $K_scale $G_scale $Gnot_scale $plot --img_filename $setuppng ${partmeshfac} $wallmeshfac >> $logfile"
    python3 $path/setup.py --setup_file $sfile $hopper_ratio $particle_rad $seed $nx $ny $L $rho_scale $K_scale $G_scale $Gnot_scale $plot --img_filename $setuppng ${partmeshfac} $wallmeshfac >> $logfile
    #cp setup.png $dir/
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
    # plotquantity='force'
    # colormap='Greys_r'

    plotquantity='damage'
    colormap='viridis'

    echo 'Generating plots'
    python3 plot_timestep.py --data_dir $dir --img_dir $dir --setup_file $sfile --dotsize $dotsize --quantity ${plotquantity} --nocolorbar --colormap $colormap
    #python3 plot_timestep.py --data_dir $dir --img_dir $dir --setup_file $sfile --dotsize $dotsize --quantity force --nocolorbar --colormap 'Greys_r' --nogrid --dpi 500
    echo "$dir"
    sxiv $dir/*.png &
}

compute_kinetic_energy(){
    echo 'Computing kinetic energy'
    python3 $path/timedata.py --directory $dir --quantity kinetic_energy --channel_y 0.0
}

#gen_vid(){
    #echo 'Generating video'
    #./gen_vid.sh $dir
#}
#
create_env
gen_setup
run
gen_plot
#gen_vid
compute_kinetic_energy
