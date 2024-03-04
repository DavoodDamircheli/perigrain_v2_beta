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
# shape = shape_dict.ring_sharp_3d()
# shape = shape_dict.torus_3d()
torus_thickness_ratio = 0.2
r2 = rad * torus_thickness_ratio/2
shape = shape_dict.torus_inscribed_3d(xyz=[0,0,0], r1=rad, r2=r2, meshsize=0.1e-3, meshdata_dir='meshdata', filename_suffix='00')
# shape = shape_dict.ring_sharp_finer_3d()

material = material_dict.peridem_3d(delta)

# material = material_dict.kalthoff3d(delta)

material.print()

# wall info
gcount = 10
# gcount = 10

safe_dist = 4 * rad - 4 * torus_thickness_ratio * rad - 2 * contact_radius * 1.001

# L = 8e-3
# L = (2 * gcount - 1) * 2 * rad + 2 * contact_radius * 1.01
L = (gcount - 1) * safe_dist  + 2 * rad + 2 * contact_radius * 1.01

# p0_x = -L + 2 * rad
# p0_y = -L + 2 * rad
p0_x = -L/2 + rad + contact_radius * 1.01
p0_y = -L/2 + rad + contact_radius * 1.01
p0_z = 0

x_min = -L/2
y_min = -L/2
z_min = -L/2
x_max = L/2
y_max = L/2
z_max = L/2
wall = Wall3d(1, x_min, y_min, z_min, x_max, y_max, z_max)

# number in the x or y direction
# gcount = 5 
count = gcount*gcount + 2 * gcount * (gcount-1) 
SL.append(shape=shape, count=count, meshsize=meshsize, material=material)

# particles = SL.generate_mesh(dimension = 3, contact_radius=contact_radius, plot_node_text=False)
particles = SL.generate_mesh(dimension = 3, contact_radius=contact_radius, plot_node_text=False, plot_mesh=False)

# particles[0][0].movable = 0
# offset for the first particle
# p0_x = -3e-3
# p0_y = -3e-3
# p0_z = -3e-3
# particles[0][0].shift([p0_x, p0_y, p0_z])

## Smaller value = more spread out
# dr = 0.8
# dr = 0.6
# dr = 0

# safe_dist = 4 * rad - 4 * torus_thickness_ratio * rad - 2 * contact_radius * 1.001
# safe_dist = 4 * rad - 4 * torus_thickness_ratio * rad
# safe_dist = 2 * (2 - dr) * rad

unmovable_p_ind = []

p_ind = 0
# distribute on x-y plane
for i in range(gcount):
    x_shift = p0_x +  safe_dist * i
    for j in range(gcount):
        # first and last layer clamped
        # if ((i==0) or (i==gcount-1)) and (j==0 or j==(gcount-1)):
        if (i in [0, gcount-1]) or (j in [0, gcount-1]):
            unmovable_p_ind.append(p_ind)
        y_shift = p0_y + safe_dist * j
        # p_ind = i * gcount + j
        print('p_ind', p_ind)
        particles[0][p_ind].shift([x_shift, y_shift, p0_z])
        p_ind += 1



# connecting rings in y-dir
for i in range(gcount):
    x_shift = p0_x + safe_dist * i  
    for j in range(gcount-1):
        y_shift = p0_y + safe_dist * j + safe_dist/2
        print('p_ind', p_ind)
        particles[0][p_ind].rotate3d('y', -np.pi/2)
        particles[0][p_ind].shift([x_shift, y_shift, p0_z])
        p_ind += 1

# connecting rings in x-dir
for i in range(gcount-1):
    x_shift = p0_x + safe_dist * i  + safe_dist/2 
    for j in range(gcount):
        y_shift = p0_y + safe_dist * j 
        print('p_ind', p_ind)
        particles[0][p_ind].rotate3d('x', -np.pi/2)
        particles[0][p_ind].shift([x_shift, y_shift, p0_z])
        p_ind += 1


for i in range(count):
    if i in unmovable_p_ind:
        pass
    else:
        particles[0][i].vel += [0, 0, args.vel_val]
        particles[0][i].acc += [0, 0, args.acc_val]

print('unmovable_p_ind', unmovable_p_ind)

for i in unmovable_p_ind:
    particles[0][i].movable = 0
    # particles[0][i].vel += [0, 0, 0]
    # particles[0][i].acc += [0, 0, 0]
    # particles[0][i].stoppable = 0
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
