## Add the project root as a path
import sys, os
sys.path.append(os.getcwd())
from random import seed
import random

import exp_dict

import numpy as np
import matplotlib.pyplot as plt

import shape_dict

# wall info
nx = 50
ny = 50
rad = 1e-3

L = rad * nx
wall_left   = -L
wall_right  = L
wall_top    = L
wall_bottom = -L


GD = exp_dict.GridDist(wall_left, wall_right, wall_top, wall_bottom, nx, ny)
incenter = GD.gridloc()
incircle_rad = np.zeros(len(incenter))

# upper and lower bound on radius
rad_1 = rad
rad_2 = rad/2

# how many types of sizes in polydispersed mixture; use polydis_n = n*nx, n=1,2 to generate left-to-right arrangement and layering
polydis_n = 5

polydis_rads = np.linspace(rad_1, rad_2, num=polydis_n, endpoint=True)

order = 'gradual_up'
# order = 'gradual_down'
# order = 'rand_unif'
# order = 'polydispersed_rnd'
# order = 'polydispersed_seq'


for i,pos in enumerate(incenter):

    if order=='gradual_up':
        # theta(y) goes from 0 to 1 as y goes from wall_bottom to wall_top
        theta = (pos[1] - wall_bottom)/(wall_top - wall_bottom)
        # convex combination
        thisrad = (1- theta) * rad_1 + theta * rad_2

    elif order=='gradual_down':
        # theta(y) goes from 0 to 1 as y goes from wall_bottom to wall_top
        theta = (pos[1] - wall_bottom)/(wall_top - wall_bottom)
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
    tt = np.linspace(0, 2*np.pi, 20, endpoint=False)
    xx = thisrad * np.cos(tt)
    yy = thisrad * np.sin(tt)
    plt.plot(xx+pos[0], yy+pos[1], 'k', linewidth=0.5)

plt.axis('scaled')
plt.grid()
#plt.savefig(filename, dpi=300, bbox_inches='tight')
plt.show()
plt.close()
