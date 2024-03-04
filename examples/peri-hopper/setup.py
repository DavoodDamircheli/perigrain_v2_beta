import numpy as np
# import time
from random import seed
import random

import sys, os
sys.path.append(os.getcwd())

import shape_dict, material_dict
# from genmesh import genmesh
from exp_dict import ShapeList, Wall, Contact, Experiment, plot_setup, GridDist
from arrangements import get_incenter_mesh_loc

import argparse
# Instantiate the parser
parser = argparse.ArgumentParser(description='Optional app description')
# Optional argument
# parser.add_argument('--shape', type=str, help='shape of each particle', default='small_disk')
parser.add_argument('--hopper_ratio', type=float, help='ratio of hopper gap to particle diameter', default=2.5)
parser.add_argument('--particle_rad', type=float, help='radius of each particle', default=8e-3)
parser.add_argument('--seed', type=int, help='seed for random number generator', default=1)
parser.add_argument('--L', type=float, help='half-length of container', default=120e-3)
parser.add_argument('--nx', type=int, help='number of particles in x dir', default=10)
parser.add_argument('--ny', type=int, help='number of particles in y dir', default=10)

parser.add_argument('--particle_meshsize_factor', type=float, help='particle meshsize factor wrt particle radius', default=3)
parser.add_argument('--wall_meshsize_factor', type=float, help='wall meshsize factor wrt particle radius', default=3)

parser.add_argument('--Gnot_scale', type=float, help='fracture toughness scale', default=0.1)
parser.add_argument('--G_scale', type=float, help='shear modulus scaling', default=0.5)
parser.add_argument('--K_scale', type=float, help='bulk modulus scaling', default=0.5)
parser.add_argument('--rho_scale', type=float, help='density scaling', default=1)

parser.add_argument('--setup_file', type=str, help='output setup directory', default='data/hdf5/all.h5')
parser.add_argument('--img_filename', type=str, help='name with path of the experiment setup image', default='setup.png')
parser.add_argument('--plot', action='store_true', help='whether to show plot or not')
# finish parsing
args = parser.parse_args()

print('plot', args.plot)
print('saving experiment setup to', args.setup_file)

particle_rad = args.particle_rad
# note: delta = (8/3)e-3 and dt = 1e-6 is stable
delta = particle_rad/2
particle_meshsize = particle_rad/args.particle_meshsize_factor
contact_radius = particle_rad/3

wall_meshsize = particle_rad/args.wall_meshsize_factor

# wall info
# L = 200e-3
L = args.L

wall_left   = -L
wall_right  = L
wall_top    = 2*L
wall_bottom = -2*L

hopper_ratio = args.hopper_ratio
hopper_gap = hopper_ratio * 2 * particle_rad 
l = (wall_right - wall_left - hopper_gap)/4
s = 1e-3

nx = args.nx
ny = args.ny


acc_val = -10
#material = material_dict.peridem_deformable(delta, G_scale=args.G_scale, K_scale=args.K_scale, Gnot_scale=0.1, rho_scale=args.rho_scale)
material = material_dict.peridem_deformable(delta, G_scale=args.G_scale, K_scale=args.K_scale, Gnot_scale=args.Gnot_scale, rho_scale=args.rho_scale)

# show_plot = False
show_plot = args.plot
dotsize = 1

# borders of particle bulk
region_left = wall_left+contact_radius/2
region_right = wall_right-contact_radius/2
region_top = wall_top-contact_radius/2
region_bottom = 0+contact_radius/2

# space around the particle
xspace = (region_right - region_left)/nx/2 - particle_rad - contact_radius/2
yspace = (region_top - region_bottom)/ny/2 - particle_rad - contact_radius/2

print('Maximum space around particle to allow perturbation without activating contact', xspace, yspace)

poslist = GridDist(region_left, region_right, region_top, region_bottom, nx, ny).gridloc()

SL = ShapeList()

# boundary shapes: 0
shape = shape_dict.plank(l=l, s=s)
SL.append(shape=shape, count=2, meshsize=wall_meshsize, material=material)

# particle shapes: 1
shape  = shape_dict.small_disk(scaling=particle_rad)
SL.append(shape=shape, count=len(poslist), meshsize=particle_meshsize, material=material)

# particles = SL.generate_mesh(dimension=2, contact_radius=contact_radius, plot_mesh=False, plot_shape=False, shapes_in_parallel=True)
particles = SL.generate_mesh(dimension=2, contact_radius=contact_radius, plot_mesh=False, plot_shape=False, shapes_in_parallel=False)

## Boundary shapes
particles[0][0].shift([-L+l,-s])
particles[0][1].shift([L-l,-s])

particles[0][0].movable = 0
particles[0][1].movable = 0
particles[0][0].stoppable = 0
particles[0][1].stoppable = 0

seed(args.seed)
print('seed', args.seed)
for i in range(len(poslist)):
    part = particles[1][i]
    part.shift(poslist[i])

    # random perturbation
    part.pos[:,0] += (-xspace) + 2 * xspace * random.uniform(0, 1)
    part.pos[:,1] += (-yspace) + 2 * yspace * random.uniform(0, 1)

    part.acc += [0, acc_val]
    part.extforce += [0, acc_val * part.material.rho]

wall = Wall(1, wall_left, wall_right, wall_top, wall_bottom)

# contact properties
normal_stiffness = 18 * material.bulk_modulus /( np.pi * np.power(delta,4));
damping_ratio = 0.9
friction_coefficient = 0.9
contact  = Contact(contact_radius, normal_stiffness, damping_ratio, friction_coefficient)

# plot
plot_setup(particles, wall=wall, contact_radius=contact_radius, delta=delta, show_plot=show_plot, dotsize=dotsize, save_filename=args.img_filename)

exp = Experiment(particles, wall, contact)

#######################################################################

# save the data
print('saving experiment setup to', args.setup_file)
exp.save(args.setup_file)
