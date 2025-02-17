#ifndef READ_CONF
#define READ_CONF 

#include <iostream>
#include <fstream>
#include <algorithm>

//#include "particle/timeloop.h"
//#include "particle/contact.h"

//using namespace std;

//template <unsigned dim>
class ConfigVal
{
public:
    unsigned timesteps;
    unsigned modulo;
    double dt;
    bool is_parallel;

    bool do_resume, wall_resume;
    unsigned resume_ind;
    bool save_file;
    bool enable_fracture;
    bool nl_bdry_only;

    bool allow_damping, allow_friction;
    double  normal_stiffness; 
    double damping_ratio, friction_coefficient;
    double new_snot;

    bool self_contact;
    double self_contact_rad;

      bool enable_cohesion=0;
      double cohesion_scaling=1;

    // torque
    bool enable_torque;
    int set_torque_timestep=-1;
    int set_fracture_timestep=-1;
    int set_self_contact_timestep=-1;

    int set_zero_wall_speed_timestep=-1;

    int set_damping_timestep=-1;

    // turn on movable tag at a given timestep
    int set_movable_index = -1;
    int set_movable_timestep = -1;
    // turn on stoppable tag at a given timestep
    int set_stoppable_index = -1;
    int set_stoppable_timestep = -1;
    // reset wheel center to the top of the bulk
    bool reset_partzero_y = 0;
    unsigned reset_partzero_y_timestep = 0;
    double wheel_rad = 0;

    //fix velocity
    int particle_w_fixed_velocity = -1;
    int set_fix_velocity_timestep = -1;
    //bool fix_velocity_x = 0;
    //bool fix_velocity_y = 0;
    //bool fix_velocity_z = 0;
    //bool fix_velocity_all = 0;
    //bool fix_velocity_mean = 0;
    bool fix_velocity_angular = 0;
    double prescribed_angular_vel = -999;
    bool fix_velocity_mean_x = 0;
    double prescribed_velocity_mean_x = -999;

    // messing with the wall
    double wall_left, wall_right, wall_top, wall_bottom;
    double wall_x_min, wall_y_min, wall_z_min, wall_x_max, wall_y_max, wall_z_max;

    double speed_wall_left, speed_wall_right, speed_wall_top, speed_wall_bottom;
    double speed_wall_x_min, speed_wall_y_min, speed_wall_z_min, speed_wall_x_max, speed_wall_y_max, speed_wall_z_max;

    //brazil nut
    bool brazil = 0;
    int brazil_reset_wall_bottom_freq;
    int brazil_nonzero_wall_bottom_speed_freq;
    float brazil_reset_wall_bottom_y;


    // force field applied
    int	forcefield_type;
    double forcefield_scaling;

    bool gradient_extforce = 0;
    unsigned extforce_maxstep = 0;

    // wall force
    bool compute_wall_reaction = 0;

    // obtained via command line
    string setup_filename = "data/hdf5/all.h5";
    string output_dir = "output/hdf5/";
    string plotinfo_filename;

    ConfigVal (){
	//Timeloop TL(8000, 100);
	//TL.dt = 0.02/1e5;
	//
	timesteps = 8000;
	modulo = 100;
	dt = 0.02/1e5;
	is_parallel = 1;

	// TL
	do_resume = 0;
	wall_resume = 0;
	resume_ind = 50;
	save_file = 1;
	enable_fracture = 0;
	nl_bdry_only = 0;

	// CN
	allow_damping = 1;
	allow_friction = 1;

	// self_contact
	self_contact = 1;
	self_contact_rad = -1;

	enable_torque = 0;

	// -1 means not specified
	normal_stiffness = -1;
	damping_ratio = -1; 
	friction_coefficient = -1;
	new_snot = -1;
	
	// default value implies not loaded from the file
       wall_left   = -999;
       wall_right  = -999;
       wall_top    = -999;
       wall_bottom = -999;
       wall_x_min = -999;
       wall_y_min = -999;
       wall_z_min = -999;
       wall_x_max = -999;
       wall_y_max = -999;
       wall_z_max = -999;

       // default value is zero, which is correct, if not defined
       speed_wall_left   = 0;
       speed_wall_right  = 0;
       speed_wall_top    = 0;
       speed_wall_bottom = 0;

       brazil_reset_wall_bottom_freq = -1;
       brazil_nonzero_wall_bottom_speed_freq = -1;
       brazil_reset_wall_bottom_y = -999;

       speed_wall_x_min = 0;
       speed_wall_y_min = 0;
       speed_wall_z_min = 0;
       speed_wall_x_max = 0;
       speed_wall_y_max = 0;
       speed_wall_z_max = 0;

       // force field applied
       forcefield_type = 0;
       forcefield_scaling = 1e9;
    };

    // Read values from a file: name = value format, (whitespace is removed), # for comment
    void read_file(std::string filename){


	std::ifstream cFile (filename);
	if (cFile.is_open())
	{
	    std::cout << "Reading config file: " << filename << std::endl;
	    std::string line;
	    while(getline(cFile, line)){
		line.erase(std::remove_if(line.begin(), line.end(), ::isspace), line.end());
		if(line[0] == '#' || line.empty())
		    continue;
		auto delimiterPos = line.find("=");
		auto name = line.substr(0, delimiterPos);
		auto value = line.substr(delimiterPos + 1);

		if (value == "") {
		    std::cout << "[Error] Missing value for option: " << name << std::endl;
		}

		if (name == "timesteps") {
		    timesteps = std::stoi(value);
		}
		else if (name == "modulo") {
		    modulo = std::stoi(value);
		}
		else if (name == "dt") {
		    dt = std::stod(value);
		}
		else if (name == "is_parallel") {
		    is_parallel = std::stoi(value);
		}
		else if (name == "do_resume") {
		    do_resume = std::stoi(value);
		}
		else if (name == "wall_resume") {
		    wall_resume = std::stoi(value);
		}
		else if (name == "resume_ind") {
		    resume_ind = std::stoi(value);
		}
		else if (name == "save_file") {
		    save_file = std::stoi(value);
		}
		else if (name == "enable_fracture") {
		    enable_fracture = std::stoi(value);
		}
		else if (name == "nl_bdry_only") {
		    nl_bdry_only = std::stoi(value);
		}

		else if (name == "particle_w_fixed_velocity") {
		    particle_w_fixed_velocity = std::stoi(value);
		}
		else if (name == "set_fix_velocity_timestep") {
		    set_fix_velocity_timestep = std::stoi(value);
		}
		//else if (name == "fix_velocity_x") {
		    //fix_velocity_x = std::stoi(value);
		//}
		//else if (name == "fix_velocity_y") {
		    //fix_velocity_y = std::stoi(value);
		//}
		//else if (name == "fix_velocity_z") {
		    //fix_velocity_z = std::stoi(value);
		//}
		//else if (name == "fix_velocity_mean") {
		    //fix_velocity_mean = std::stoi(value);
		//}
		else if (name == "fix_velocity_angular") {
		    fix_velocity_angular = std::stoi(value);
		}
		else if (name == "prescribed_angular_vel") {
		    prescribed_angular_vel = std::stod(value);
		}
		else if (name == "fix_velocity_mean_x") {
		    fix_velocity_mean_x = std::stoi(value);
		}
		else if (name == "prescribed_velocity_mean_x") {
		    prescribed_velocity_mean_x = std::stod(value);
		}

		else if (name == "allow_damping") {
		    allow_damping = std::stoi(value);
		}
		else if (name == "allow_friction") {
		    allow_friction = std::stoi(value);
		}
		else if (name == "damping_ratio") {
		    damping_ratio = std::stod(value);
		}
		else if (name == "friction_coefficient") {
		    friction_coefficient = std::stod(value);
		}
		else if (name == "normal_stiffness") {
		    normal_stiffness = std::stod(value);
		}
		else if (name == "new_snot") {
		    new_snot = std::stod(value);
		}

		else if (name == "forcefield_type") {
		    forcefield_type = std::stoi(value);
		}
		else if (name == "forcefield_scaling") {
		    forcefield_scaling = std::stod(value);
		}

		else if (name == "self_contact") {
		    self_contact = std::stoi(value);
		}
		else if (name == "self_contact_rad") {
		    self_contact_rad = std::stod(value);
		}

		else if (name == "enable_torque") {
		    enable_torque = std::stoi(value);
		}
		else if (name == "set_torque_timestep") {
		    set_torque_timestep = std::stoi(value);
		}
		else if (name == "set_fracture_timestep") {
		    set_fracture_timestep = std::stoi(value);
		}
		else if (name == "set_self_contact_timestep") {
		    set_self_contact_timestep = std::stoi(value);
		}

		else if (name == "set_damping_timestep") {
		    set_damping_timestep = std::stoi(value);
		}

		else if (name == "set_zero_wall_speed_timestep") {
		    set_zero_wall_speed_timestep = std::stoi(value);
		}

		else if (name == "enable_cohesion") {
		    enable_cohesion = std::stoi(value);
		}
		else if (name == "cohesion_scaling") {
		    cohesion_scaling = std::stod(value);
		}


		else if (name == "set_movable_index") {
		    set_movable_index = std::stoi(value);
		}
		else if (name == "set_movable_timestep") {
		    set_movable_timestep = std::stoi(value);
		}
		else if (name == "set_stoppable_index") {
		    set_stoppable_index = std::stoi(value);
		}
		else if (name == "set_stoppable_timestep") {
		    set_stoppable_timestep = std::stoi(value);
		}

		else if (name == "reset_partzero_y") {
		     reset_partzero_y = std::stoi(value);
		}
		else if (name == "reset_partzero_y_timestep") {
		     reset_partzero_y_timestep = std::stoi(value);
		}
		else if (name == "wheel_rad") {
		     wheel_rad = std::stod(value);
		}

		else if (name == "wall_left") {
		    wall_left = std::stod(value);
		}
		else if (name == "wall_right") {
		    wall_right = std::stod(value);
		}
		else if (name == "wall_top") {
		    wall_top = std::stod(value);
		}
		else if (name == "wall_bottom") {
		    wall_bottom = std::stod(value);
		}
		else if (name == "speed_wall_left") {
		    speed_wall_left = std::stod(value);
		}
		else if (name == "speed_wall_right") {
		    speed_wall_right = std::stod(value);
		}
		else if (name == "speed_wall_top") {
		    speed_wall_top = std::stod(value);
		}
		else if (name == "speed_wall_bottom") {
		    speed_wall_bottom = std::stod(value);
		}

		// for brazil nut exp
		else if (name == "brazil") {
		    brazil = std::stoi(value);
		}
		else if (name == "brazil_reset_wall_bottom_freq") {
		    brazil_reset_wall_bottom_freq = std::stoi(value);
		}
		else if (name == "brazil_nonzero_wall_bottom_speed_freq") {
		    brazil_nonzero_wall_bottom_speed_freq = std::stoi(value);
		}
		else if (name == "brazil_reset_wall_bottom_y") {
		    brazil_reset_wall_bottom_y = std::stod(value);
		}

		else if (name == "wall_x_min") {
		    wall_x_min = std::stod(value);
		}
		else if (name == "wall_y_min") {
		    wall_y_min = std::stod(value);
		}
		else if (name == "wall_z_min") {
		    wall_z_min = std::stod(value);
		}
		else if (name == "wall_x_max") {
		    wall_x_max = std::stod(value);
		}
		else if (name == "wall_y_max") {
		    wall_y_max = std::stod(value);
		}
		else if (name == "wall_z_max") {
		    wall_z_max = std::stod(value);
		}
		else if (name == "speed_wall_x_min") {
		    speed_wall_x_min = std::stod(value);
		}
		else if (name == "speed_wall_y_min") {
		    speed_wall_y_min = std::stod(value);
		}
		else if (name == "speed_wall_z_min") {
		    speed_wall_z_min = std::stod(value);
		}
		else if (name == "speed_wall_x_max") {
		    speed_wall_x_max = std::stod(value);
		}
		else if (name == "speed_wall_y_max") {
		    speed_wall_y_max = std::stod(value);
		}
		else if (name == "speed_wall_z_max") {
		    speed_wall_z_max = std::stod(value);
		}

		else if (name == "gradient_extforce") {
		    gradient_extforce = std::stoi(value);
		}
		else if (name == "extforce_maxstep") {
		    extforce_maxstep = std::stoi(value);
		}

		else if (name == "compute_wall_reaction") {
		    compute_wall_reaction = std::stoi(value);
		}


		else{
		    std::cerr << "[Error]: Wrong config name: " << name << " !!\n";
		}
	    }
	    
	    std::cout << "Reading config file: Done" << std::endl;
	}
	else {
	    std::cerr << "Couldn't open config file for reading.\n";
	}
    };

    void print(){
	std::cout << "timesteps: " << timesteps << std::endl;
	std::cout << "modulo: "  << modulo << std::endl;
	std::cout << "dt: " << dt << std::endl;
	std::cout << "is_parallel: " << is_parallel << std::endl;
	std::cout << "do_resume: " << do_resume << std::endl;
	std::cout << "wall_resume: " << wall_resume << std::endl;
	std::cout << "resume_ind: " << resume_ind << std::endl;
	std::cout << "save_file: " << save_file << std::endl;
	std::cout << "enable_fracture: " << enable_fracture << std::endl;
	std::cout << "nl_bdry_only: " << nl_bdry_only << std::endl;
	std::cout << "allow_damping: " << allow_damping << std::endl;
	std::cout << "allow_friction: " << allow_friction << std::endl;
	std::cout << "damping_ratio: " << damping_ratio << std::endl;
	std::cout << "friction_coefficient: " << friction_coefficient << std::endl;
	std::cout << "normal_stiffness: " << normal_stiffness << std::endl;
	std::cout << "new_snot" << new_snot << std::endl;

	//std::cout << "speed_wall_top: " << speed_wall_top << std::endl;
	std::cout << "----------------------------------------------------------------------" << std::endl;
    };

private:
    /* data */
};


#endif /* ifndef READ_CONF */
