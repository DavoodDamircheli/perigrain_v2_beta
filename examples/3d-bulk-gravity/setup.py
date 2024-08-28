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
# import pdb
# pdb.set_trace()
parser = argparse.ArgumentParser(description='Optional app description')

""" A bulk of particles, locations are generated from a mesh
Each particle mesh is generated using scaling factor, hence the neighborhood and boundary nodes are computed separately
"""
# Optional argument
# parser.add_argument('--shape', type=str, help='shape of each particle', default='small_disk')
parser.add_argument('--particle_rad', type=float, help='radius of each particle', default=8e-3)
parser.add_argument('--L', type=float, help='half-length of container', default=100e-3)
parser.add_argument('--wallh_ratio', type=float, help='half-length of wall height vs wall width', default=2)
parser.add_argument('--nx', type=int, help='number of particles in x dir', default=10)
parser.add_argument('--ny', type=int, help='number of particles in y dir', default=10)
parser.add_argument('--vel_val', type=float, help='initial velocity', default=-20)
parser.add_argument('--acc_val', type=float, help='initial velocity', default=-10)

# 'hollow_sphere', 'plus3d', 'kalthoff'
parser.add_argument('--shape', type=str, help='shape of grain', default='cube')
parser.add_argument('--msh_path', type=str, help='shape of grain', default='/home/davood/projects/beta_perigrain_v2/grain-data')
parser.add_argument('--G_scale', type=float, help='shear modulus scaling', default=0.5)
parser.add_argument('--Gnot_scale', type=float, help='shear modulus scaling', default=1e-4)
parser.add_argument('--K_scale', type=float, help='bulk modulus scaling', default=0.5)
parser.add_argument('--rho_scale', type=float, help='density scaling', default=1)
parser.add_argument('--meshsize_factor', type=float, help='meshsize factor compared to radius', default=50)
parser.add_argument('--delta_factor', type=float, help='delta factor compared to radius', default=3)
parser.add_argument('--contact_rad_factor', type=float, help='contact radius factor compared to radius', default=5)
parser.add_argument('--setup_file', type=str, help='output setup directory', default='data/hdf5/all.h5')
parser.add_argument('--plot', action='store_true', help='whether to show plot or not')
# finish parsing
args = parser.parse_args()
args.plot=1
print('plot', args.plot)
print('saving experiment setup to', args.setup_file)

""" Two particles colliding in 3D
"""

rad = 100e-3
lx = rad
#n_x, n_y, n_z = (4, 4, 4)
n_x, n_y, n_z = (1,1,2)
mf = args.meshsize_factor
mf = 4

count = 2
delta = rad/args.delta_factor
meshsize = rad/mf
contact_radius = rad/args.contact_rad_factor    # conserves momentum better (than delta/3)

L = n_x * (rad + contact_radius)
Ly = n_y * (rad + contact_radius)
Lz = n_z * (rad + contact_radius)

gap = rad + contact_radius
x = np.linspace(-L+gap, L-gap, n_x)
y = np.linspace(-Ly+gap, Ly-gap, n_y)
z = np.linspace(-Lz+gap, Lz-gap, n_z)
X, Y, Z = np.meshgrid(x, y, z)
xf = X.flatten()
yf = Y.flatten()
zf = Z.flatten()
shifts = np.c_[xf, yf, zf]

SL = ShapeList()

material = material_dict.peridem_3d(delta)

material.print()

x_min = -L
y_min = -Ly
z_min = -Lz
x_max = L
y_max = Ly
z_max = Lz
wall = Wall3d(1, x_min, y_min, z_min, x_max, y_max, z_max)
#------------------max count is 400--------------
if args.shape == 'grains':
    for i in range(1, count + 1):
        print(i)
        ind = int(i)
        msh_file = args.msh_path + '/mesh_' + str(ind) + '.msh'
        shape = shape_dict.Shape(P=None, nonconvex_interceptor=None, msh_file=msh_file, scale_mesh_to=rad, centroid_origin=True)
        SL.append(shape=shape, count=1, meshsize=meshsize, material=material, plot_shape=False, shape_type=args.shape)

elif args.shape == 'plus3d':
    shape = shape_dict.plus3d(l=lx, meshsize=meshsize)
    SL.append(shape=shape, count=count, meshsize=meshsize, material=material, plot_shape=False,shape_type=args.shape)

elif args.shape == 'cube':
    shape = shape_dict.cube_3d(lx=lx, ly=lx, lz=lx, meshsize=meshsize)
    SL.append(shape=shape, count=count, meshsize=meshsize, material=material, plot_shape=False,shape_type=args.shape)

elif args.shape == 'sphere':
    shape = shape_dict.sphere_3d(rad=lx, meshsize=meshsize)
    SL.append(shape=shape, count=count, meshsize=meshsize, material=material, plot_shape=False,shape_type=args.shape)

elif args.shape == 'hollow_sphere':
    shape = shape_dict.hollow_sphere(R=lx, r=lx/2, meshsize=meshsize)
    SL.append(shape=shape, count=count, meshsize=meshsize, material=material, plot_shape=False,shape_type=args.shape)


particles = SL.generate_mesh(dimension = 3, contact_radius=contact_radius, plot_node_text=False, plot_shape=False, plot_mesh=False)

## Apply transformation
seed(1)
g_val = -5e4
for i in range(count):
    #particles[0][i].scale(scaling_list[i])

    ## generate random rotation
    particles[0][i].rotate3d('x',0 + (random() * (2*np.pi - 0)) )
    particles[0][i].shift(shifts[i])
    particles[0][i].acc += [0, 0,g_val]
    particles[0][i].extforce += [0,0, g_val * particles[0][i].material.rho]

normal_stiffness = material_dict.peridem_3d(contact_radius).cnot / contact_radius



damping_ratio = 0.8
friction_coefficient = 0.8
contact  = Contact(contact_radius, normal_stiffness, damping_ratio, friction_coefficient)

if args.plot:
    plot3d_setup(particles, dotsize=15, wall=wall, show_particle_index=True, delta=delta, contact_radius=contact_radius, save_filename='setup.png')

exp = Experiment(particles, wall, contact)

#######################################################################

# save the data
print('saving experiment setup to', args.setup_file)
exp.save(args.setup_file)
