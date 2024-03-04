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

# shape = shape_dict.disk_w_hole_3d()
# shape = shape_dict.torus_3d()
# shape = shape_dict.ring_sharp_3d()

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


import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--setup', action='store_true', help='whether to show plot from setup.h5')
# parser.add_argument('--filename', type=str, help='input filename', default='/home/debdeep/from_smic2/tc_00001.h5')
# parser.add_argument('--filename', type=str, help='input filename', default='/home/debdeep/peri-wheel/examples_output/3d-chain-pattern/nofrac_/tc_00101.h5')
parser.add_argument('--path', type=str, help='parent path', default='/home/debdeep/peri-wheel/examples_output/3d-chain-pattern/nofrac_')
parser.add_argument('--ind', type=int, help='filename index', default=1)
parser.add_argument('--shininess', default=100)
parser.add_argument('--wireframe-width', default=0.5)
args, _ = parser.parse_known_args()

canvas = scene.SceneCanvas(keys='interactive', bgcolor='white')
view = canvas.central_widget.add_view()

view.camera = 'arcball'
# view.camera.depth_value = 1e3
view.camera.distance = 100e-3 # since particle is small, we need to zoom in a lot

wireframe_filter = WireframeFilter(width=args.wireframe_width)
shading_filter = ShadingFilter(shininess=args.shininess)


# another with shift
# mesh = Mesh(vertices+np.array([0,0,1e-3]), faces, color=(.5, .7, .5, 1))
# view.add(mesh)
# mesh.attach(wireframe_filter)
# mesh.attach(shading_filter)


if args.setup:
    filename = args.path+'/setup.h5'
else:
    fname = ('/tc_%05d.h5' % args.ind)
    filename = args.path+fname

f = h5py.File(filename, "r")

# # Create a colored `MeshVisual`.
# Pos = np.array(f['P_00000'+'/CurrPos'], dtype=float)
# print('vert', vertices, len(vertices))
# print('pos', Pos, len(Pos))
# # mesh = Mesh(vertices, faces, color=(.5, .7, .5, 0.9))
# mesh = Mesh(Pos, faces, color=(.5, .7, .5, 0.9))
# view.add(mesh)
# mesh.attach(wireframe_filter)
# mesh.attach(shading_filter)


for name in f:
    if re.match(r'P_[0-9]+', name):
        pid = int(name[2:])
        print(pid)

# for i in range(25):
#         print('i=',i)
#         name = ('P_%05d' % i)


        if args.setup:
            Pos = np.array(f[name+'/Pos'], dtype=float)
        else:
            Pos = np.array(f[name+'/CurrPos'], dtype=float)

        # Create a colored `MeshVisual`.
        mesh = Mesh(Pos, faces, color=(.5, .7, .5, 1))
        view.add(mesh)
        mesh.attach(wireframe_filter)
        mesh.attach(shading_filter)

def attach_headlight(view):
    light_dir = (0, 1, 0, 0)
    shading_filter.light_dir = light_dir[:3]
    initial_light_dir = view.camera.transform.imap(light_dir)

    @view.scene.transform.changed.connect
    def on_transform_change(event):
        transform = view.camera.transform
        shading_filter.light_dir = transform.map(initial_light_dir)[:3]


attach_headlight(view)

shading_states = (
    dict(shading=None),
    dict(shading='flat'),
    dict(shading='smooth'),
)
shading_state_index = shading_states.index(
    dict(shading=shading_filter.shading))

wireframe_states = (
    dict(wireframe_only=False, faces_only=False,),
    dict(wireframe_only=True, faces_only=False,),
    dict(wireframe_only=False, faces_only=True,),
)
wireframe_state_index = wireframe_states.index(dict(
    wireframe_only=wireframe_filter.wireframe_only,
    faces_only=wireframe_filter.faces_only,
))


def cycle_state(states, index):
    new_index = (index + 1) % len(states)
    return states[new_index], new_index


@canvas.events.key_press.connect
def on_key_press(event):
    global shading_state_index
    global wireframe_state_index
    if event.key == 's':
        state, shading_state_index = cycle_state(shading_states,
                                                 shading_state_index)
        for attr, value in state.items():
            setattr(shading_filter, attr, value)
        mesh.update()
    elif event.key == 'w':
        wireframe_filter.enabled = not wireframe_filter.enabled
        mesh.update()
    elif event.key == 'f':
        state, wireframe_state_index = cycle_state(wireframe_states,
                                                   wireframe_state_index)
        for attr, value in state.items():
            setattr(wireframe_filter, attr, value)
        mesh.update()


canvas.show()


if __name__ == "__main__":
    app.run()
