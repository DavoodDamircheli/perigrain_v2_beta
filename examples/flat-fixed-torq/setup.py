import numpy as np
# import time
from random import seed
from random import random

import sys, os
sys.path.append(os.getcwd())

import shape_dict, material_dict
# from genmesh import genmesh
from exp_dict import ShapeList, Wall, Contact, Experiment, plot_setup
from arrangements import get_incenter_mesh_loc

import argparse
# Instantiate the parser
parser = argparse.ArgumentParser(description='Optional app description')
# Optional argument
parser.add_argument('--shape', type=str, help='shape of grain', default='pertdisk')
parser.add_argument('--shape_param_1', type=float, help='parameter for grain geometry', default=0.4)
parser.add_argument('--setup_file', type=str, help='output setup directory', default='data/hdf5/all.h5')
parser.add_argument('--wheel_rad', type=str, help='output filename containing wheel radius', default='output/wheel_rad')
parser.add_argument('--plot', action='store_true', help='whether plot or not')
# finish parsing
args = parser.parse_args()

print('plot', args.plot)
print('saving experiment setup to', args.setup_file)

def wheel_on_inclined():
    """ Wheel on inclined plane
    """

    # whether or not to add particles
    wheel_rad = 3e-3
    inner_circle_ratio = 0.7

    torque_val = -17

    # delta = 1e-3
    # delta = 10e-3
    delta = 1e-3
    # meshsize = delta/2
    meshsize = delta/7
    contact_radius = delta/4

    # blade_l = 5e-3
    # blade_s = 0.5e-3
    blade_l = 12e-3
    blade_s = 1e-3

    # inclination 
    # blade_angle = -np.pi/10
    blade_angle = 0

    # wheel location
    l_offset = -2e-3
    s_offset = 0e-3
    # particle toughness

    # wheel properties
    wheel_rho_scale = 1
    # modify plank properties
    plank_rho_scale = 0.1
    plank_G_scale =0.1
    plank_K_scale = 0.1
    plank_Gnot_scale = 0.7

    # at least this much space to leave between particles. Will get halved while selecting max radius
    min_particle_spacing = contact_radius*1.001

    # wall info
    # cl = 6e-3
    # cl = 10e-3 # 10 cm
    cl = blade_l*1.5

    wall_left   = -cl
    wall_right  = cl
    wall_top    = cl
    wall_bottom = -cl
    # wall_bottom = -5e-3
    wall = Wall(1, wall_left, wall_right, wall_top, wall_bottom)

    # particle generation boundary
    c = min_particle_spacing/2
    P = np.array([
        [wall_left + c, wall_bottom + c ],
        [wall_right - c, wall_bottom + c],
        # [wall_right - c, wall_top - c],
        # [wall_left + c, wall_top - c]
        ## only bottom half
        [wall_right - c, 0 - blade_s - c],
        [wall_left + c, 0-blade_s - c]

        ])


    # Create a list of shapes
    SL = ShapeList()

    # wheel
    # SL.append(shape=shape_dict.small_disk(scaling=wheel_rad) , count=1, meshsize=meshsize, material=material_dict.peridem_deformable(delta))
    # SL.append(shape=shape_dict.small_disk(scaling=wheel_rad) , count=1, meshsize=meshsize, material=material_dict.peridem_deformable(delta))
    # SL.append(shape=shape_dict.pygmsh_geom_test(scaling=wheel_rad, meshsize=meshsize) , count=1, meshsize=meshsize, material=material_dict.peridem_deformable(delta))
    # SL.append(shape=shape_dict.gmsh_test(scaling=wheel_rad, meshsize=meshsize/2) , count=1, meshsize=meshsize, material=material_dict.peridem_deformable(delta))
    SL.append(shape=shape_dict.wheel_annulus(scaling=wheel_rad, inner_circle_ratio=inner_circle_ratio, meshsize=meshsize, filename_suffix='00') , count=1, meshsize=meshsize, material=material_dict.peridem_deformable(delta, rho_scale=wheel_rho_scale, K_scale=1, G_scale=1))
    # SL.append(shape=shape_dict.annulus() , count=1, meshsize=meshsize, material=material_dict.peridem_deformable(delta))
    # SL.append(shape=shape_dict.wheel_ring(scaling=2e-3) , count=1, meshsize=meshsize, material=material_dict.peridem_deformable(delta))
    SL.append(shape=shape_dict.plank(l=blade_l, s=blade_s) , count=1, meshsize=meshsize, material=material_dict.peridem_deformable(delta, rho_scale=plank_rho_scale, K_scale=plank_K_scale, G_scale=plank_G_scale))


    # generate the mesh for each shape
    particles = SL.generate_mesh(dimension = 2, contact_radius = contact_radius, plot_mesh=False, plot_shape=False, shapes_in_parallel=False)

    print('Here')

    # plank, rotation
    plank = particles[1][0]
    plank.rotate(blade_angle)

    ## want deformable
    particles[1][0].movable =0
    particles[1][0].stoppable =0
    ## clamped the bottom
    # for i in range(len(plank.pos)):
        # bottom_edge = -blade_s
        # if np.abs(plank.pos[i][1] - bottom_edge) < 1e-10 :
            # plank.clamped_nodes.append(i)


    radial_dir = np.array([-np.cos(-blade_angle) , -np.sin(-blade_angle)])
    tang_dir = np.array([-np.sin(-blade_angle) , np.cos(-blade_angle)])
    shift_vec = (blade_l + l_offset)  * radial_dir + (blade_s + wheel_rad + contact_radius+s_offset) * tang_dir
    # shift_vec = blade_l  * tang_dir
    print('shift_val', shift_vec)
    particles[0][0].shift( shift_vec)


    ## applying torque
    wheel = particles[0][0]
    wheel.breakable = 0
    wheel.stoppable = 1

    #torque density specified
    wheel.torque_val = torque_val * wheel.material.rho
    wheel.axis = 2 # z-axis

    #gravity
    g_val = -1e4
    wheel.extforce += [0, g_val * wheel.material.rho]

    # initial angular velcity
    # 1000 RPM (pretty high for burr coffee grinder) = 1000 * 2 * pi / 60 (rad/s) ~ 104 rad/s
    # w_val = -600
    w_val = 0
    perp = np.array([[0, -1], [1, 0]])
    centroid = np.mean(wheel.pos + wheel.disp, axis=0)
    print('centroid', centroid)
    for j in range(len(wheel.pos)):
        r_vec = (wheel.pos[j] - centroid)
        r = np.sqrt( np.sum(r_vec**2) )
        u_dir = r_vec / r
        t_dir = perp @ u_dir
        wheel.vel[j] += w_val * r * t_dir


    # contact properties
    normal_stiffness = 18 * material_dict.peridem(delta).bulk_modulus /( np.pi * np.power(delta,4));
    damping_ratio = 0.8
    friction_coefficient = 0.8
    contact  = Contact(contact_radius, normal_stiffness, damping_ratio, friction_coefficient)

    if args.plot:
        plot_setup(particles, dotsize=1, contact_radius=contact_radius)

    return Experiment(particles, wall, contact)

#######################################################################

exp = wheel_on_inclined()
# save the data
print('saving experiment setup to', args.setup_file)
exp.save(args.setup_file)
