import sys, os
sys.path.append(os.getcwd())

import numpy as np
import h5py

import re
# from multiprocessing import Pool
from pathos.multiprocessing import ProcessingPool as Pool
import time

import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
from matplotlib.collections import LineCollection
from sys import argv


import argparse
# Instantiate the parser
parser = argparse.ArgumentParser(description='Optional app description')
# Optional argument
parser.add_argument('--quantity', type=str, help='quantity to compute', default='kinetic_energy')
parser.add_argument('--directory', type=str, help='directory location of data', default='output/hdf5')
parser.add_argument('--channel_y', type=float, help='y-value to be considered the hopped opening', default=0.02)
# finish parsing
args = parser.parse_args()

#directory = 'examples_output/hopper/nofrac__x-5_y-5_r-20e-3_z-5_w-5'
#args.directory = 'examples_output/hopper/nofrac__z-8_w-8_K-0.00001_G-0.00001'
#args.quantity = 'kinetic_energy'
#quantity = 'mass_loss'


# plot info
#p = h5py.File('output/hdf5/plotinfo.h5', "r")
p = h5py.File(args.directory+'/plotinfo.h5', "r")
fc = int(p['f_l_counter'][0])
lc = int(p['f_l_counter'][1])
dt = float(p['dt'][0])
modulo = int(p['modulo'][0])

def total_neighbors(conn, N):
    """Compute the total number of neighbors from connectivity data
    :conn: connectivity matrix, mx2, m=total intact bonds
    :N: total number of nodes
    :returns: TODO
    """
    deg = np.zeros(N)
    for i in range(len(conn)):
        deg[conn[i][0]] += 1
        deg[conn[i][1]] += 1
    return deg

def mass_below(t):
    print(t)
    tc_ind = ('%05d' % t)

    #filename = 'output/hdf5/tc_'+tc_ind+'.h5'
    filename = args.directory+'/tc_'+tc_ind+'.h5'
    f = h5py.File(filename, "r")

    #orig_filename = 'data/hdf5/all.h5'
    orig_filename = args.directory+'/setup.h5'
    orig_f = h5py.File(orig_filename, "r")

    tot_vol = 0
    for name in f:
        if re.match(r'P_[0-9]+', name):
            pid = int(name[2:])
            Pos = np.array(f[name+'/CurrPos'])
            Vol = np.array(orig_f[name+'/Vol'])

            low_y_ind = np.squeeze(Pos[:,1] <= args.channel_y)
            part_out_vols = np.sum(Vol[low_y_ind])

            tot_vol += part_out_vols

    return tot_vol

def kinetic_energy_below(t):
    print(t)
    tc_ind = ('%05d' % t)

    #filename = 'output/hdf5/tc_'+tc_ind+'.h5'
    filename = args.directory+'/tc_'+tc_ind+'.h5'
    f = h5py.File(filename, "r")

    #orig_filename = 'data/hdf5/all.h5'
    orig_filename = args.directory+'/setup.h5'
    orig_f = h5py.File(orig_filename, "r")

    tot_vol = 0
    for name in f:
        if re.match(r'P_[0-9]+', name):
            pid = int(name[2:])
            Pos = np.array(f[name+'/CurrPos'])
            vel = np.array(f[name+'/vel'])
            Vol = np.array(orig_f[name+'/Vol'])

            low_y_ind = np.squeeze(Pos[:,1] <= args.channel_y)

            #part_out_vols = np.sum(Vol[low_y_ind])
            kinetic_energy = np.sum(0.5 * Vol[low_y_ind] * np.sum(vel[low_y_ind]**2, axis=0))

            #tot_vol += part_out_vols
            tot_vol += kinetic_energy

    return tot_vol

tt = range(fc, lc+1)

## serial j
# for t in range(fc, lc+1):
    # genplot(t)


## parallel
a_pool = Pool()
if args.quantity=='kinetic_energy':
    t_vols = a_pool.map(kinetic_energy_below, tt)
elif args.quantity=='mass_loss':
    t_vols = a_pool.map(mass_below, tt)
else:
    print('Wrong quantity provided')
a_pool.close()

out_png = args.directory+'/'+args.quantity+'.png'
plt.plot(tt, t_vols)
plt.title(args.quantity)
# plt.show()
print('Saving img to', out_png)
plt.savefig(out_png, dpi=300, bbox_inches='tight')
