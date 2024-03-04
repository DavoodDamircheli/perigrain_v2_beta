import numpy as np
# import time
from random import seed
import random
import matplotlib.pyplot as plt


import sys, os
sys.path.append(os.getcwd())

import shape_dict, material_dict
# from genmesh import genmesh
from exp_dict import ShapeList, Wall, Contact, Experiment, plot_setup, plot_node_bond_distribution, GridDist
from arrangements import get_incenter_mesh_loc

from shape_params import Param

from getargs import readargs
args = readargs()

## comments: delta=30e-3 with meshsize 50e-3/4 blows up, but delta=40e-3 works until timestep 600, them blows up
# # delta = 45e-3
# # grain_meshsize = 50e-3/4
# # contact_radius = 50e-3/4
delta = args.delta
meshsize = args.grain_meshsize 
contact_radius = args.contact_radius

g_val = -1e1

# at least this much space to leave between particles. Will get halved while selecting max radius
min_particle_spacing = contact_radius*1.001

## wall half-length
# L = 3000e-3
# hh = 1200e-3
L = args.L
hh = args.hh
# particle generation spacing
P_meshsize = args.P_meshsize
if args.manual_P_trim:
    min_trim_rad = args.P_trim_val
else:
    min_trim_rad = P_meshsize / 1000 * 2.5


#### material
wheel_material = material_dict.peridem_deformable(delta, G_scale=1, rho_scale=args.wheel_rho_scale)
grain_material = material_dict.peridem_deformable(delta, G_scale=args.grain_G_scale, Gnot_scale=args.grain_Gnot_scale, rho_scale=args.grain_rho_scale)


######## Wheel geometry
## Wheel radius: standard is 14 inch
# # wheel_rad = hh/2
if args.fit_wheel:
    wheel_rad = hh/args.fit_wheel_scale
else:
    wheel_rad = args.wheel_rad
# wheel_inner_ratio = 0.8
wheel_inner_ratio = 0.9
# wheel_x_offset = 2*contact_radius
wheel_x_offset = 2*wheel_rad

## torque density
# torque_val = -1

# save wheel radius to file to load into the config file via script
with open(args.wheel_rad_file, 'w') as f:
    f.write(str(wheel_rad))


#######################################################################
## Wall 

wall_left   = -L
wall_right  = L
wall_top    = hh * 1.2
# wall_top    = hh * 2
wall_bottom = -hh

bulk_wall_left   = -L
bulk_wall_right  = L
bulk_wall_top    = 0
bulk_wall_bottom = -hh

# wall_bottom = -5e-3
wall = Wall(1, wall_left, wall_right, wall_top, wall_bottom)


seed(args.seed)

#######################################################################
## Particle arrangement

if args.arrangement == 'interior':
    # particle generation boundary
    c = min_particle_spacing/2
    P = np.array([
        [bulk_wall_left + c, bulk_wall_bottom + c ],
        [bulk_wall_right - c, bulk_wall_bottom + c],
        [bulk_wall_right - c, bulk_wall_top - c],
        [bulk_wall_left + c, bulk_wall_top - c]
        ])
    msh = get_incenter_mesh_loc(P, P_meshsize, modify_nodal_circles= True, gen_bdry = False )
    # reduce radius to avoid contact
    msh.incircle_rad -= min_particle_spacing/2

    # trim minimum radius
    msh.trim(min_rad=min_trim_rad)

    msh.info()
    # msh.plot(plot_edge = True)
    msh.plot_distribution(bins=20, filename=args.size_distribution_pngfile)
    incenter = msh.incenter
    incircle_rad = msh.incircle_rad

elif args.arrangement == 'grid':
    c = min_particle_spacing/2

    nx =  int( (bulk_wall_right - bulk_wall_left)/ (2 * (args.max_rad + contact_radius/2) ) )
    ny =  int( (bulk_wall_top - bulk_wall_bottom)/ (2 * (args.max_rad + contact_radius/2) ) )
    # nx =  int( (bulk_wall_right - bulk_wall_left)/ (2 * (args.max_rad) ) )
    # ny =  int( (bulk_wall_top - bulk_wall_bottom)/ (2 * (args.max_rad) ) )

    print('nx, ny', nx, ny)

    GD = GridDist(bulk_wall_left, bulk_wall_right, bulk_wall_top, bulk_wall_bottom, nx, ny)
    incenter = GD.gridloc()
    print('length', len(incenter))
    incircle_rad = np.zeros(len(incenter))

    # how many types of sizes in polydispersed mixture; use polydis_n = n*nx, n=1,2 to generate left-to-right arrangement and layering
    polydis_n = 5

    rad_1 = args.min_rad
    rad_2 = args.max_rad
    polydis_rads = np.linspace(rad_1, rad_2, num=polydis_n, endpoint=True)

    # order = 'gradual'
    # order = 'rand_unif'
    # order = 'polydispersed_rnd'
    # args.arr_grid_order = 'polydispersed_seq'
    order = args.arr_grid_order

    ## paste here

    for i,pos in enumerate(incenter):

        if order=='gradual_up':
            # theta(y) goes from 0 to 1 as y goes from wall_bottom to wall_top
            theta = (pos[1] - bulk_wall_bottom)/(bulk_wall_top - bulk_wall_bottom)
            # convex combination
            thisrad = (1- theta) * rad_1 + theta * rad_2

        elif order=='gradual_down':
            # theta(y) goes from 0 to 1 as y goes from bulk_wall_bottom to bulk_wall_top
            theta = (pos[1] - bulk_wall_bottom)/(bulk_wall_top - bulk_wall_bottom)
            # convex combination
            thisrad = (1- theta) * rad_2 + theta * rad_1

        elif order == 'rand_unif':
            thisrad = random.uniform(rad_1, rad_2)

        elif order == 'polydispersed_rnd':
            thisrad = random.choice(polydis_rads)

        elif order == 'polydispersed_seq':
            thisrad = polydis_rads[i % polydis_n]


        incircle_rad[i] = thisrad

        # draw circle
        # tt = np.linspace(0, 2*np.pi, 20, endpoint=False)
        # xx = thisrad * np.cos(tt)
        # yy = thisrad * np.sin(tt)
        # plt.plot(xx+pos[0], yy+pos[1], 'k', linewidth=0.5)

    print('[', np.min(incircle_rad), ',', np.max(incircle_rad), ']', np.mean(incircle_rad))
    # plt.axis('scaled')
    # plt.grid()
    # #plt.savefig(filename, dpi=300, bbox_inches='tight')
    # plt.show()
    # plt.close()

else:
    print('Wrong particle arrangement type')


print('total particles', len(incircle_rad))


#######################################################################
## Shapes

# Create a list of shapes
SL = ShapeList()

# add a wheel
# SL.append(shape=shape_dict.wheel_annulus(scaling=wheel_rad, inner_circle_ratio=wheel_inner_ratio, meshsize=args.wheel_meshsize, filename_suffix='wheel') , count=1, meshsize=args.wheel_meshsize, material=wheel_material)

# append each shape with own scaling
for i in range(len(incenter)):
    if args.shape == 'plus':
        shape = shape_dict.plus_inscribed(notch_dist=args.plus_notch_dist, scaling=incircle_rad[i])

    elif args.shape == 'pertdisk':
        shape = shape_dict.pertdisk_incirc(seed=i, steps=args.pertdisk_steps, scaling=incircle_rad[i], incirc_ratio=args.pertdisk_incirc_ratio)

    elif args.shape == 'n4':
        shape = shape_dict.rect(scaling=incircle_rad[i], ratio_vol=args.n4_ratio_vol)

    elif args.shape == 'circ':
        shape=shape_dict.small_disk(steps=args.circ_steps, scaling=incircle_rad[i])

    elif args.shape == 'ring':
        shape = shape_dict.wheel_annulus(scaling=incircle_rad[i], inner_circle_ratio=args.ring_inner_ratio, meshsize=args.grain_meshsize, filename_suffix=str(i))

    else:
        print('Wrong input shape specified')

    SL.append(shape=shape, count=1, meshsize=args.grain_meshsize, material=grain_material)

# generate the mesh for each shape
particles = SL.generate_mesh(dimension = 2, contact_radius = contact_radius, plot_mesh=False, plot_shape=False, shapes_in_parallel=True)

# for the wheel
# wheel = particles[0][0]
# wheel.breakable = 0
# wheel.shift([-L+wheel_rad+contact_radius+wheel_x_offset, bulk_wall_top+wheel_rad+contact_radius/2])
# wheel.acc += [0, g_val]
# wheel.extforce += [0, g_val * wheel.material.rho]

# torque
# wheel.torque_axis = 2
# wheel.torque_val = torque_val * wheel.material.rho

# keep the wheel fixed in the beginning, mobilize later externally via config file and scripts
# wheel.movable = 0
# wheel.stoppable = 0

## Apply transformation
for i in range(len(incircle_rad)):
    ## scaling is done while generating the mesh
    ## generate random rotation
    part = particles[i+0][0]
    part.rotate(0 + (random.random() * (2*np.pi - 0)) )
    part.shift(incenter[i])

    # Initial data
    part.vel += [0, -1]
    part.acc += [0, g_val]
    part.extforce += [0, g_val * part.material.rho]

# contact properties
normal_stiffness = 18 * material_dict.peridem(delta).bulk_modulus /( np.pi * np.power(delta,4));
damping_ratio = 0.8
friction_coefficient = 0.8
contact  = Contact(contact_radius, normal_stiffness, damping_ratio, friction_coefficient)


show_plot = args.plot
plot_setup(particles, dotsize=args.dotsize, contact_radius=contact_radius, delta=delta, show_plot=args.plot, wall=wall, linewidth=0.1, save_filename=args.pngfile)

# return Experiment(particles, wall, contact)
exp =  Experiment(particles, wall, contact)

#######################################################################

# save the data
print('saving experiment setup to', args.setup_file)
exp.save(args.setup_file)
