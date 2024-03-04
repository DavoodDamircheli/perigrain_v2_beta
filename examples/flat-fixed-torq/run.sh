#!/bin/bash
path="`dirname \"$0\"`"
logfile=$path/output.log

resume="no"
#resume="yes"

cohesion="no"
#cohesion="yes"
cohesion_scaling="0.2"

shapes=(
#0
#0.1
#0.2
#0.3
0.4
#0.6
#ring
#plus
#ring0.2
#ring0.4
)

#clear logfile
echo '' > $logfile

function run {
    for i in "${!shapes[@]}"
    do 
	shape=${shapes[i]}
	echo "#####################################################"
	echo "Using shape: $shape"

	# create subdirectory
	dir=${path}/data/${shape}_${1}
	mkdir -p $dir

	base=$path/flat_base.conf

	config=$dir/main.conf
	sfile=$dir/setup.h5

	## generate experiment setup
	python3 $path/setup.py --shape $shape --setup_file $sfile >> $logfile
	#python3 $path/setup.py --shape $shape --setup_file $sfile  --plot >> $logfile

	# construct config file
	cp $base $config

	# enable fracture of not
	if [ "$1" = "frac" ]
	then
	    echo 'Enable fracture and self_contact.'

	    echo "enable_fracture = 1" >> $config
	    echo "self_contact = 1" >> $config
	else
	    echo 'Disable fracture and self_contact.'

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

	    #echo "# set a given paticle (index = 0) to movable (it was previously not movable)" >> $config
	    #echo "set_movable_index = 0" >> $config
	    #echo "set_movable_timestep = $last" >> $config
	    #echo "set_stoppable_index = 0" >> $config
	    #echo "set_stoppable_timestep = $last" >> $config
	else
	    echo ""
	    #echo "# set a given paticle (index = 0) to movable (it was previously not movable)" >> $config
	    #echo "set_movable_index = 0" >> $config
	    #echo "set_movable_timestep = 15000" >> $config
	    #echo "set_stoppable_index = 0" >> $config
	    #echo "set_stoppable_timestep = 15000" >> $config

	    #echo "# reset particle zero position to bulk height at a given timestep" >> $config
	    #echo "reset_partzero_y = 1" >> $config
	    #echo "reset_partzero_y_timestep = 15000" >> $config
	fi

	if [ "$cohesion" = "yes" ]
	then
	    echo "enable_cohesion = 1" >> $config
	    echo "cohesion_scaling = $cohesion_scaling" >> $config
	fi
	
	# run code
	echo 'running'
	bin/simulate2d -c $config -o $dir -i $sfile  >> $logfile

	echo 'generation images'
	python3 plot_timestep.py --data_dir $dir --img_dir $dir --setup_file $sfile

	# angular velocity
	python3 wheel_attrib.py --data_dir $dir --img_dir $dir --setup_file $sfile

	# generate video
	#./gen_vid.sh $dir
    done
}

# call function
run 'frac'
#run ''
