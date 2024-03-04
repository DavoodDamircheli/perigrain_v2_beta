## Add the project root as a path
import sys, os
sys.path.append(os.getcwd())

import exp_dict

import numpy as np
import matplotlib.pyplot as plt

import shape_dict

# wall info
wall_left   = -10e-3
wall_right  = 10e-3
wall_top    = 10e-3
wall_bottom = 0e-3

# particle generation boundary
# P = np.array([
    # [wall_left, wall_bottom],
    # [wall_right, wall_bottom],
    # [wall_right, wall_top],
    # [wall_left, wall_top]
    # ])

# shape=shape_dict.small_disk(steps=30)
#shape=shape_dict.small_disk_fromfile()
# shape=shape_dict.test()
# shape=shape_dict.vase(l=10e-3, s=5e-3, steps=20, n_pi=3, phase_1=0, phase_2=np.pi/2)
#shape = shape_dict.star(steps = 20, scaling=1e-3, incirc_ratio=0.3)

# shape = shape_dict.wheel_annulus(scaling=10e-3, meshsize=1e-3, inner_circle_ratio=0.5, filename_suffix='00', nci_steps=2)
shape = shape_dict.wheel_annulus_local(scaling=10e-3, meshsize=1e-3, inner_circle_ratio=0.5, filename_suffix='00', nci_steps=2)
## doesn't work
# shape = shape_dict.gear_inscribed_par(scaling = 10e-3, teethn=5, inner_circle_ratio=0.4, filename_suffix='00', meshsize=1e-3, hollow=False, hollow_inner_rad=0.3)
## doesn't work
# shape = shape_dict.belt_reverse_gear(scaling = 10e-3, teethn=5, half_thickness_ratio=0.8, inner_circle_ratio=0.7, filename_suffix='00', meshsize=1e-3, tooth_gap_ratio=0.6)

msh_file = shape.msh_file
# particle generation spacing
P_meshsize = 3.5e-3

# print(shape)
# P = shape.P
# print(P)
# plt.plot(P[:,0], P[:,1])
# plt.show()

# msh = exp_dict.get_incenter_mesh_loc(P, P_meshsize, msh_file=None, modify_nodal_circles= False, gen_bdry = False )

# for shapes with no P and only .msh
msh = exp_dict.get_incenter_mesh_loc(None, P_meshsize, msh_file=msh_file, modify_nodal_circles= False, gen_bdry = False )

# shape = shape_dict.test()
# msh_file = shape.msh_file
# print(msh_file)
# msh = exp_dict.get_incenter_mesh_loc(P=[], meshsize=[], msh_file=msh_file, dimension=2, modify_nodal_circles= False, gen_bdry = False )

savefile='output/arr_no_modify.png'
# def plot(self, do_plot = True, plot_edge = True, plot_mesh = True, remove_axes = True, save_file = None):
msh.plot(save_file = savefile, plot_mesh=True, remove_axes=True, plot_circles=False)
msh.info()
# plt.show()
plt.close()


# msh = exp_dict.get_incenter_mesh_loc(P, P_meshsize, modify_nodal_circles= True, gen_bdry = False )
# msh.plot(save_file = '/home/debdeep/gdrive/work/peridynamics/granular/data/incircles_nodal.png' )
# plt.close()

# msh = exp_dict.get_incenter_mesh_loc(P, P_meshsize, modify_nodal_circles= True, gen_bdry = True )
msh = exp_dict.get_incenter_mesh_loc(None, P_meshsize, msh_file=msh_file, modify_nodal_circles=True, gen_bdry = True )
savefile='output/arr_modified.png'
msh.plot(save_file = savefile, plot_mesh=False, remove_axes=True, plot_circles=True)
msh.info()
# plt.show()
plt.close()
