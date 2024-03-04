import numpy as np
import sys, os
sys.path.append(os.getcwd())

import shape_dict,material_dict
from exp_dict import ShapeList, Wall3d, Contact, Experiment, plot3d_setup,plot_node_bond_distribution, plot_bonds,plot_nodes, plot_bonds_nonConvex

import pdb
#pdb.set_trace()
#--------kalthoff plate---------
import argparse
parser = argparse.ArgumentParser()

parser.add_argument('--delta_frac',type=float, help='cancer',default=4)

parser.add_argument('--plot',action='store_true', help='cancer2')

args = parser.parse_args()
lx = 200e-3
ly = 100e-3
lz = 9e-3
nl = 50e-3                    # notch_lenght 
na = 15               #notch_angle 
radianA =  np.deg2rad(na)
ncor = 75e-3     #notch_dist_corner 

#----- we need some error for convex tol
rate_cnx = 2
l_trig = 2 * nl *np.tan(radianA/2) 
convex_tol = (lx-(2*(l_trig+ncor)))/rate_cnx

convex_tol = 3e-3 
meshSize  = lx/3
delta = meshSize*3

#delta      = meshSize*delta_frac
Gnot_scale =1
rho_scale  =1

sl = ShapeList()

#shape1 = shape_dict.kalthoff_3d_plate_beta0(lx,ly,lz, nl, na,ncor,convex_tol, meshSize) 

#def kalthoff_3d_plate_beta(ix=6e-3, iy=2e-3, iz=.2e-3, Ln=1e-3, Wn=0.5e-3,distance_from_edge=1e-3, meshsize=.3e-3, meshdata_dir='meshdata', filename_suffix='00'):
notch_lenght = 50e-3
notch_w = 15.5e-3
notch_dist_corner = 70e-3
ncor = notch_dist_corner


shape1 = shape_dict.kalthoff_3d_plate_beta(lx,ly,lz, notch_lenght,notch_w, notch_dist_corner,meshSize)


#shape1 = shape_dict.kalthoff_3d_plate()
material=material_dict.kalthoff(delta, Gnot_scale=Gnot_scale, rho_scale=rho_scale)
material.print()

sl.append(shape=shape1, count=1, meshsize=meshSize,material=material )
particles = sl.generate_mesh(dimension = 3, plot_mesh=False) 

#plot_bonds_nonConvex(particles)
plot_bonds(particles)
#plot_nodes(particles)   # this one does nt work for new appoach of removing 

#plot3d_setup(particles)
