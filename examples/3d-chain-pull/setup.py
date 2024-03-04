import numpy as np
# import time
from random import seed
from random import random

import sys, os
sys.path.append(os.getcwd())

import shape_dict, material_dict
# from genmesh import genmesh
from exp_dict import ShapeList, Wall3d, Contact, Experiment, plot3d_setup, GridDist
from arrangements import get_incenter_mesh_loc

from shape_params import Param

import argparse
# Instantiate the parser
parser = argparse.ArgumentParser(description='Optional app description')
# Optional argument
# parser.add_argument('--shape', type=str, help='shape of each particle', default='small_disk')
parser.add_argument('--particle_rad', type=float, help='radius of each particle', default=8e-3)
parser.add_argument('--L', type=float, help='half-length of container', default=100e-3)
parser.add_argument('--wallh_ratio', type=float, help='half-length of wall height vs wall width', default=2)
parser.add_argument('--nx', type=int, help='number of particles in x dir', default=10)
parser.add_argument('--ny', type=int, help='number of particles in y dir', default=10)
parser.add_argument('--vel_val', type=float, help='initial velocity', default=-10)
parser.add_argument('--acc_val', type=float, help='initial velocity', default=-10)

parser.add_argument('--G_scale', type=float, help='shear modulus scaling', default=0.5)
parser.add_argument('--Gnot_scale', type=float, help='shear modulus scaling', default=1e-4)
parser.add_argument('--K_scale', type=float, help='bulk modulus scaling', default=0.5)
parser.add_argument('--rho_scale', type=float, help='density scaling', default=1)

parser.add_argument('--meshsize_factor', type=float, help='meshsize factor compared to radius', default=4)
parser.add_argument('--delta_factor', type=float, help='delta factor compared to radius', default=3)
parser.add_argument('--contact_rad_factor', type=float, help='contact radius factor compared to radius', default=4)

parser.add_argument('--setup_file', type=str, help='output setup directory', default='data/hdf5/all.h5')
parser.add_argument('--plot', action='store_true', help='whether to show plot or not')
# finish parsing
args = parser.parse_args()

print('plot', args.plot)
print('saving experiment setup to', args.setup_file)

""" Two particles colliding in 3D
"""

rad = 1e-3

delta = rad/args.delta_factor
meshsize = rad/args.meshsize_factor
contact_radius = rad/args.contact_rad_factor    # conserves momentum better (than delta/3)

SL = ShapeList()

# shape=shape_dict.sphere_3d(rad=rad, meshsize=meshsize)
# shape = shape=shape_dict.plus_small_3d()
# shape = shape_dict.disk_w_hole_3d()
shape = shape_dict.ring_sharp_3d()
# shape = shape_dict.ring_sharp_finer_3d()

material = material_dict.peridem_3d(delta)

# material = material_dict.kalthoff3d(delta)

# count = 3
count = 4
SL.append(shape=shape, count=count, meshsize=meshsize, material=material)
# SL.append(shape=shape_dict.sphere_small_3d(), count=2, meshsize=meshsize, material=material_dict.peridem_3d(delta))
# SL.append(shape=shape_dict.sphere_small_3d(), count=2, meshsize=meshsize, material=material_dict.peridem_3d(delta))
# SL.append(shape=shape_dict.disk_w_hole_3d(), count=2, meshsize=meshsize, material=material_dict.peridem_3d(delta))
# SL.append(shape=shape_dict.plus_small_3d(), count=2, meshsize=meshsize, material=material_dict.peridem_3d(delta))

material.print()

# wall info
L = 8e-3
x_min = -L
y_min = -L
z_min = -L
x_max = L
y_max = L
z_max = L
wall = Wall3d(1, x_min, y_min, z_min, x_max, y_max, z_max)

# particles = SL.generate_mesh(dimension = 3, contact_radius=contact_radius, plot_node_text=False)
particles = SL.generate_mesh(dimension = 3, contact_radius=contact_radius, plot_node_text=False, plot_mesh=False)

# particles[0][0].movable = 0
# offset for the first particle
p0_y = -3e-3
# p0_z = 3e-3
# p0_y = 0e-3
p0_z = 0e-3
particles[0][0].shift([0, p0_y, p0_z])

dr = 0.8
for i in range(1, count):
    if i % 2 == 1:
        particles[0][i].rotate3d('y', -np.pi/2)
    particles[0][i].shift([0, p0_y + (2*i - i*dr)*rad, p0_z])

particles[0][0].vel += [0, args.vel_val, 0]
particles[0][count-1].vel += [0, -args.vel_val, 0]
# particles[0][i].acc += [0, 0, args.acc_val]


# contact properties
# normal_stiffness = 18 * material_dict.peridem(delta).bulk_modulus /( np.pi * np.power(delta,5));
# normal_stiffness = material.cnot / contact_radius
normal_stiffness = material_dict.peridem_3d(contact_radius).cnot / contact_radius
# normal_stiffness = 15 * mat.E /( np.pi * np.power(delta,5) * (1 - 2*mat.nu));

damping_ratio = 0.8
friction_coefficient = 0.8

contact  = Contact(contact_radius, normal_stiffness, damping_ratio, friction_coefficient)

plot3d_setup(particles, dotsize=15, wall=wall, show_particle_index=True, delta=delta, contact_radius=contact_radius)

exp = Experiment(particles, wall, contact)

#######################################################################

# save the data
print('saving experiment setup to', args.setup_file)
exp.save(args.setup_file)
