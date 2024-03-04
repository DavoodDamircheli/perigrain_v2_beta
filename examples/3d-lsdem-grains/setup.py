import sys, os
sys.path.append(os.getcwd())

import shape_dict, material_dict
from shape_dict import Shape
from exp_dict import ShapeList, Wall3d, Contact, Experiment, plot3d_setup, GridDist
from arrangements import get_incenter_mesh_loc
from shape_params import Param
from genmesh import genmesh

L = 20e-3
meshsize = 1e-3
meshsize_sc = 10
# meshsize_sc = 5
R = 1e-3

ind = 68
msh_file = '/home/debdeep/ls-shapes/msh/mesh_' + str(ind) + '.msh'
shape = Shape(P=None, nonconvex_interceptor=None, msh_file=msh_file)
shape.plot()

# shape = shape_dict.torus_inscribed_3d(r1=R, r2=0.1*R, meshsize=R/meshsize_sc)
# shape = shape_dict.box_3d([-R/2, -R/2, -R/2], R, R, R, meshsize=R/meshsize_sc)
# shape = shape_dict.box_extrude_3d([-R/2, -R/2, -R/2], R, R, R, meshsize=R/meshsize_sc)

mesh = genmesh(P_bdry=None, meshsize=None, msh_file=shape.msh_file, dimension=3)
trisurf = False
transparent = False
mesh.plot(plot_node_text=False, highlight_bdry_nodes=True, trisurf=trisurf, transparent=transparent, volume_proportional=False)
# mesh.plot()

