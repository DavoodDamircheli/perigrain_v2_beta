import sys, os
sys.path.append(os.getcwd())

import shape_dict, material_dict
from exp_dict import ShapeList, Wall3d, Contact, Experiment, plot3d_setup, GridDist
from arrangements import get_incenter_mesh_loc
from shape_params import Param
from genmesh import genmesh

L = 20e-3
meshsize = 1e-3
meshsize_sc = 10
# meshsize_sc = 5
R = 1e-3

# shape = shape_dict.cube_3d(L=L, meshsize=meshsize)
# shape = shape_dict.cube_extrude_3d(L=L, meshsize=meshsize)
# shape = shape_dict.sheet_rect_w_gap_3d(L=L, gap_l_ratio=0.5, meshsize=meshsize)
# shape = shape_dict.disk_w_hole_prog_3d(R=3e-3, height_sc=0.2, gamma=0.2, meshsize=R/meshsize_sc)
# shape = shape_dict.rect_w_cut_3d(R=10e-3, h_sc=0.05, meshsize=0.2e-3)
# shape = shape_dict.torus_3d(r1=R-0.2*R, r2=0.2*R, meshsize=R/meshsize_sc)
shape = shape_dict.torus_inscribed_3d(r1=R, r2=0.1*R, meshsize=R/meshsize_sc)
# shape = shape_dict.box_3d([-R/2, -R/2, -R/2], R, R, R, meshsize=R/meshsize_sc)
# shape = shape_dict.box_extrude_3d([-R/2, -R/2, -R/2], R, R, R, meshsize=R/meshsize_sc)

mesh = genmesh(P_bdry=None, meshsize=None, msh_file=shape.msh_file, dimension=3)
trisurf = True
transparent = True
mesh.plot(plot_node_text=False, highlight_bdry_nodes=False, trisurf=trisurf, transparent=True, volume_proportional=False)

