import numpy as np
# import time
from random import seed
from random import random
import time

import sys
sys.path.insert(0, "../")
import sys, os
sys.path.append(os.getcwd())

import shape_dict, material_dict
from genmesh import genmesh
from exp_dict import ShapeList, Wall, Contact, Experiment, GridDist
from arrangements import get_incenter_mesh_loc
from shape_params import Param
from genmesh import genmesh

import h5py
import re

#######################################################################
# Mesh
rad = 1e-3
torus_thickness_ratio = 0.2
r2 = rad * torus_thickness_ratio/2
shape = shape_dict.torus_inscribed_3d(xyz=[0,0,0], r1=rad, r2=r2, meshsize=0.1e-3, meshdata_dir='meshdata', filename_suffix='00')
mesh = genmesh(P_bdry=None, meshsize=None, msh_file=shape.msh_file, dimension = 3)

vertices, faces = mesh.pos, mesh.bdry_edges
#######################################################################

from vispy import app, scene
from vispy.io import read_mesh, load_data_file
from vispy.scene.visuals import Mesh
from vispy.scene import transforms
from vispy.visuals.filters import ShadingFilter, WireframeFilter
import vispy.io as io


import argparse
parser = argparse.ArgumentParser()

# parser.add_argument('--filename', type=str, help='input filename', default='/home/debdeep/from_smic2/tc_00001.h5')
# parser.add_argument('--filename', type=str, help='input filename', default='/home/debdeep/from_smic2/tc_00500.h5')
parser.add_argument('--path', type=str, help='parent path', default='/home/debdeep/peri-wheel/examples_output/3d-chain-pattern/nofrac_')
parser.add_argument('--ind', type=int, help='filename index', default=1)

parser.add_argument('--shininess', default=100)
parser.add_argument('--wireframe-width', default=0.5)
args, _ = parser.parse_known_args()

def attach_headlight(view):
    light_dir = (0, 1, 0, 0)
    shading_filter.light_dir = light_dir[:3]
    initial_light_dir = view.camera.transform.imap(light_dir)

    @view.scene.transform.changed.connect
    def on_transform_change(event):
        transform = view.camera.transform
        shading_filter.light_dir = transform.map(initial_light_dir)[:3]

canvas = scene.SceneCanvas(keys='interactive', bgcolor='white')
view = canvas.central_widget.add_view()

view.camera = 'arcball'
# view.camera.depth_value = 1e3
view.camera.distance = 35e-3 # since particle is small, we need to zoom in a lot

wireframe_filter = WireframeFilter(width=args.wireframe_width)
wireframe_filter.enabled = False    # turn off the lines
shading_filter = ShadingFilter(shading='smooth', shininess=args.shininess)


fname = ('/tc_%05d.h5' % args.ind)
filename = args.path+fname
f = h5py.File(filename, "r")

# canvas.show()
# for i in range(90):
#         print('i=',i)
#         name = ('P_%05d' % i)
#
for name in f:
    if re.match(r'P_[0-9]+', name):
        pid = int(name[2:])
        print(pid)

        Pos = np.array(f[name+'/CurrPos'], dtype=float)

        # Create a colored `MeshVisual`.
        mesh = Mesh(Pos, faces, color=(.5, .7, .5, 1))
        view.add(mesh)
        mesh.attach(wireframe_filter)
        mesh.attach(shading_filter)


attach_headlight(view)
img = canvas.render()
vname = ('/vp_%05d.png' % args.ind)
io.write_png(args.path+vname,img)
