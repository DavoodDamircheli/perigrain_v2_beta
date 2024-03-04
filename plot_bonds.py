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
parser.add_argument('--data_dir', type=str, help='input data directory', default='output/hdf5/')
parser.add_argument('--img_dir', type=str, help='output image directory', default='output/img/')
parser.add_argument('--setup_file', type=str, help='output setup file', default='data/hdf5/all.h5')
parser.add_argument('--all_dir', type=str, help='directory where all files are located')
parser.add_argument('--img_prefix', type=str, help='string prefix to go in the filename', default='img')

parser.add_argument('--fc', type=int, help='first counter')
parser.add_argument('--lc', type=int, help='last counter')
parser.add_argument('--dotsize', type=float, help='dotsize', default=1)
parser.add_argument('--vmin', type=float, help='lower limit of colorbar', default=None)
parser.add_argument('--vmax', type=float, help='upper limit of colorbar', default=None)
parser.add_argument('--colormap', type=str, help='quantity to plot in color', default='viridis')
parser.add_argument('--xlim', type=str, help='space-separated xlim, within quote. e.g. --xlim \'-1 1\'')
parser.add_argument('--ylim', type=str, help='same as ylim')
parser.add_argument('--nogrid', action='store_true', help='do not show grid')
parser.add_argument('--nowall', action='store_true', help='do not show wall')
parser.add_argument('--nocolorbar', action='store_true', help='do not show colorbar')
parser.add_argument('--seaborn', action='store_true', help='use seaborn package')

parser.add_argument('--quantity', type=str, help='quantity to plot in color', default='damage')
parser.add_argument('--serial', action='store_true', help='read timesteps in serial')
parser.add_argument('--plot_contact_rad', action='store_true', help='plot contact radius')
parser.add_argument('--contact_rad_ind', type=float, help='particle on which to draw contact radius', default=0)
parser.add_argument('--contact_rad_node', type=float, help='node on which to draw contact radius', default=0)
parser.add_argument('--adaptive_dotsize', action='store_true', help='scale dotsizes by volume element')
parser.add_argument('--colorize_only', nargs='+', help='only colorize particle indices, separate with space', required=False)


# for plotting contact forces
parser.add_argument('--plot_contact_force', action='store_true', help='plot nodes experiencing contact force')
parser.add_argument('--wheel_ind', type=int, help='index of the particle experiencing contact force', default=0)
parser.add_argument('--nbr_ind', type=int, help='index of neighboring particle with which particle[wheel_ind] is in contact', default=4)
# zoomed in near the wheel
parser.add_argument('--zoom_wheel', action='store_true', help='create plot near the wheel. Uses --wheel_ind.')

# finish parsing
args = parser.parse_args()

data_dir = args.data_dir
img_dir = args.img_dir
setup_filename = args.setup_file
dotsize = args.dotsize

if args.all_dir:
    data_dir       = args.all_dir
    img_dir        = args.all_dir
    setup_filename = args.all_dir+'/setup.h5'

# plot_damage = 1

# plot info
plotinfo_file = data_dir+'/plotinfo.h5'
p = h5py.File(plotinfo_file, "r")
fc = int(p['f_l_counter'][0])
lc = int(p['f_l_counter'][1])
dt = float(p['dt'][0])
modulo = int(p['modulo'][0])

if args.fc:
    fc = args.fc
if args.lc:
    lc = args.lc

if args.seaborn:
    import seaborn as sns
    sns.set()

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


# ref_connectivity 
f = h5py.File(setup_filename, "r")
ref_n = {}
for name in f:
    if re.match(r'P_[0-9]+', name):
        # pid = int(name[2:])
        ref_connectivity = np.array(f[name+'/Connectivity']).astype(int)
        N = len(f[name+'/Pos'])
        orig_tot_nbrs = total_neighbors(ref_connectivity, N)
        ref_n[name] = orig_tot_nbrs

        if args.zoom_wheel:
            pid = int(name[2:])
            if pid == args.wheel_ind:
                Pos = np.array(f[name+'/Pos'])
                wheel_centroid = np.mean(Pos, axis=0, keepdims=True)
                wheel_rad = np.max(np.sqrt(np.sum((Pos - wheel_centroid)**2, axis=1, keepdims=False)))

contact_radius = np.array(f['pairwise/contact_radius']).astype(float)[0][0]

if args.adaptive_dotsize:
    N = (f['total_particles_univ']).astype(int)[0][0]
    slot = np.zeros(N)
    for name in f:
        if re.match(r'P_[0-9]+', name):
            pid = int(name[2:])
            vol = np.array(f[name+'/Vol'])
            slot[pid] = np.min(vol)
    min_vol = np.min(slot)
    print(min_vol)

    scaled_ds = [[]] * N
    for name in f:
        if re.match(r'P_[0-9]+', name):
            pid = int(name[2:])
            vol = np.array(f[name+'/Vol'])
            scaled_ds[pid] = np.sqrt(vol/min_vol) * float(dotsize)



def genplot(t):
    print(t)
    tc_ind = ('%05d' % t)
    filename = data_dir+'/tc_'+tc_ind+'.h5'
    wall_filename = data_dir+'/wall_'+tc_ind+'.h5'
    # out_png = img_dir+'/img_'+tc_ind+'.png'
    out_png = img_dir+'/'+args.img_prefix+'_'+tc_ind+'.png'
    f = h5py.File(filename, "r")

    wheel_centroid = []

    for name in f:
        if re.match(r'P_[0-9]+', name):
            pid = int(name[2:])
            Pos = np.array(f[name+'/CurrPos'])

            if args.zoom_wheel:
                if pid == args.wheel_ind:
                    wheel_centroid = np.mean(Pos, axis=0, keepdims=False)
                
            c = 'blue'
            if (args.quantity == 'damage'):
                # print('plotting damage')
                # damage
                P_conn = np.array(f[name+'/Connectivity']).astype(int)
                N = len(Pos)
                now_nbrs = total_neighbors(P_conn, N)
                orig_nbrs = np.array(ref_n[name])
                damage = (orig_nbrs - now_nbrs)/ orig_nbrs
                c = damage
                # print('damage', damage)
                vmin=0
                vmax=1
            if (args.quantity == 'force'):
                # print('plotting damage')
                # damage
                force = np.array(f[name+'/force'])
                f_norm = np.sqrt(np.sum(force**2, axis=1))
                c = f_norm
                # print('damage', damage)
                vmin=args.vmin
                vmax=args.vmax

            if args.colorize_only:
                indx = [int(ii) for ii in args.colorize_only]
                if pid in indx:
                    pass
                else:
                    c = 'blue'
                    vmin = 0
                    vmax = 1

 
            P_conn = np.array(f[name+'/Connectivity']).astype(int)

            E0 = P_conn[:,0]
            E1 = P_conn[:,1]

            P0 = Pos[E0]
            P1 = Pos[E1]

            ls = [ [p0, p1] for p0,p1 in zip(P0, P1)]
            # print('lc', lc)
            
            
            # a = [wi[0],  wi[3]]
            # b = [wi[1], wi[3]]
            # c = [wi[1], wi[2]]
            # d = [wi[0],  wi[2]]
            # ls = [ [a, b], [b,c], [c,d], [d,a] ]
            lc = LineCollection(ls, linewidths=0.1, colors='b')
            plt.gca().add_collection(lc)

            # plt.scatter(Pos[:,0], Pos[:,1], c=c, s=dotsize, marker='.', linewidth=0, cmap='viridis', vmin=vmin, vmax=vmax)

            if args.adaptive_dotsize:
                dval = scaled_ds[pid]
            else:
                dval = dotsize
            plt.scatter(Pos[:,0], Pos[:,1], c=c, s=dval, marker='.', linewidth=0, cmap=args.colormap, vmin=vmin, vmax=vmax)

            if args.plot_contact_force:
                if pid == args.wheel_ind:
                    cnode_list = []
                    nbr_name = ('P_%05d' % args.nbr_ind)
                    nbr_currpos = np.array(f[nbr_name+'/CurrPos'])
                    for node in range(len(Pos)):
                        for i in range(len(nbr_currpos)):
                            dist = np.sqrt(np.sum((nbr_currpos[i] - Pos[node])**2))
                            if dist <= contact_radius:
                                cnode_list.append(node)
                                break
                    # plt.scatter(Pos[cnode_list,0], Pos[cnode_list,1], color='red')
                    plt.scatter(Pos[cnode_list,0], Pos[cnode_list,1], color='red', s=dval, marker='.', linewidth=0, cmap=args.colormap)


            if args.plot_contact_rad:
                if pid == args.contact_rad_ind:
                    node = args.contact_rad_node
                    x0 = Pos[node][0]
                    y0 = Pos[node][1]
                    steps=20
                    tt = np.linspace(0, 2*np.pi, steps, endpoint=True)
                    xx = contact_radius * np.cos(tt)
                    yy = contact_radius * np.sin(tt)
                    plt.plot( xx + x0, yy + y0)

    # wall
    if args.nowall:
        pass
    else:
        w = h5py.File(wall_filename, "r")
        wi = np.array(w['wall_info'])[:,0]
        a = [wi[0],  wi[3]]
        b = [wi[1], wi[3]]
        c = [wi[1], wi[2]]
        d = [wi[0],  wi[2]]
        ls = [ [a, b], [b,c], [c,d], [d,a] ]
        lc = LineCollection(ls, linewidths=1, colors='b')
        plt.gca().add_collection(lc)

    if not args.nocolorbar:
        plt.colorbar()
    if not args.nogrid:
        plt.grid()
    if args.xlim:
        print('xlim', args.xlim)
        # xx = args.xlim.split(',')
        xx = args.xlim.split(' ')
        print(xx[0], xx[1])
        plt.xlim(float(xx[0]), float(xx[1]))
    if args.ylim:
        # yy = args.ylim.split(',')
        yy = args.ylim.split(' ')
        plt.ylim(float(yy[0]), float(yy[1]))

    plt.axis('scaled')
    plt.savefig(out_png, dpi=300, bbox_inches='tight')

    if args.zoom_wheel:
        # find wheel centroid
        # print(wheel_centroid)
        # print('wheel_rad', wheel_rad)
        # specify boundary
        plt.xlim(wheel_centroid[0] - 1.5* wheel_rad, wheel_centroid[0] + 1.5 *wheel_rad)
        plt.ylim(wheel_centroid[1] - 3*wheel_rad, wheel_centroid[1] + 1.5 * wheel_rad)

        # scaled axes (note: if xlim/ylim specified, axis('scaled') won't work
        plt.gca().set_aspect('equal', adjustable='box')

        zoom_out_png = img_dir+'/zoom_'+args.img_prefix+'_'+tc_ind+'.png'
        plt.savefig(zoom_out_png, dpi=300, bbox_inches='tight')

    plt.close()

if args.serial:
    ## serial j
    for t in range(fc, lc+1):
        genplot(t)
else:
    ## parallel
    a_pool = Pool()
    a_pool.map(genplot, range(fc, lc+1))
    a_pool.close()
