#!/bin/bash

shape=grains
path="`dirname \"$0\"`"
hpc="no"


if [ "$hpc" = "yes" ]
then
	output_path="/work/$USER/XXXX/$dirname"
else
	output_path="examples_output/compress/hollow_sphere"
fi

echo "dir is the output path for results"
dir="$output_path/N2"

echo "$dir"



#-----------------copyimg some config files in dir------
config=$dir/main.conf
sfile=$dir/setup.h5
logfile=$dir/output.log
cp $path/base.conf $config
cp $path/setup.py $dir/
mkdir -p  $dir



#------------------------Generate setup----------------------------
gen_setup(){
	echo "we are here"
	python3 -u $path/setup.py --shape=$shape --setup_file $dir/setup.h5  

}
#----------------------MPIRUN-----------------------------

run(){

    mpirun -n 3 bin/simulate3d -c $path/base.conf -i $dir/setup.h5 -o $dir
}




#---------------Genrate plots----------------
gen_plot(){

	python3 plot3d_timestep.py --all_dir $dir --dotsize 30 
}
gen_setup
run
gen_plot





