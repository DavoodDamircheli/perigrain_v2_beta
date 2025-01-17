import numpy as np
import h5py
import re
from pathos.multiprocessing import ProcessingPool as Pool
import time
import matplotlib.pyplot as plt
import matplotlib
#matplotlib.use('Agg')
matplotlib.use('TkAgg')
from matplotlib.collections import LineCollection
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection
from matplotlib.colors import LinearSegmentedColormap
import pdb
import argparse
import seaborn as sns
# Instantiate the parser
parser = argparse.ArgumentParser(description='Optional app description')

# Optional arguments
parser.add_argument('--data_dir', type=str, help='input data directory', default='output/hdf5/')
parser.add_argument('--img_dir', type=str, help='output image directory', default='output/img/')
parser.add_argument('--setup_file', type=str, help='output setup file', default='data/hdf5/all.h5')
parser.add_argument('--all_dir', type=str, help='directory where all files are located')
parser.add_argument('--img_prefix', type=str, help='string prefix to go in the filename', default='img')

parser.add_argument('--fc', type=int, help='first counter')
parser.add_argument('--lc', type=int, help='last counter')
parser.add_argument('--dotsize', type=float, help='dotsize', default=1)
parser.add_argument('--alpha', type=float, help='alpha', default=0.5)
parser.add_argument('--scale_factor', type=float, help='scale_factor', default=1)
parser.add_argument('--cross_section', type=float, help='scale_factor', default=0)
parser.add_argument('--filter_yellow', type=float, help='scale_factor', default=0)

parser.add_argument('--kalthoff', type=float, help='scale_factor', default=0)
parser.add_argument('--colormap', type=str, help='quantity to plot in color- including coolwarm or Spectral ', default='viridis')
parser.add_argument('--xlim', type=str, help='space-separated xlim, within quote. e.g. --xlim \'-1 1\'')
parser.add_argument('--ylim', type=str, help='same as ylim')
parser.add_argument('--zlim', type=str, help='only for 3d')
parser.add_argument('--nogrid', action='store_true', help='do not show grid')
parser.add_argument('--nowall', action='store_true', help='do not show wall')
parser.add_argument('--nocolorbar', action='store_true', help='do not show colorbar')
parser.add_argument('--seaborn', action='store_true', help='use seaborn package')

parser.add_argument('--view_az', type=float, help='camera azimuthal', default=30)
parser.add_argument('--view_angle', type=float, help='camera horizontal angle', default=45)
parser.add_argument('--motion', action='store_true', help='camera motion', default=1)
parser.add_argument('--motion_az_stepsize', type=float, help='camera azimuthal step size', default=0)
parser.add_argument('--motion_angle_stepsize', type=float, help='camera horizontal angle step size', default=0)

parser.add_argument('--quantity', type=str, help='quantity to plot in color', default='damage')
parser.add_argument('--serial', action='store_true', help='read timesteps in serial')
parser.add_argument('--plot_contact_rad', action='store_true', help='plot contact radius')
parser.add_argument('--contact_rad_ind', type=float, help='particle on which to draw contact radius', default=0)
parser.add_argument('--contact_rad_node', type=float, help='node on which to draw contact radius', default=0)
parser.add_argument('--adaptive_dotsize', action='store_true', help='scale dotsizes by volume element')
parser.add_argument('--colorize_only', nargs='+', help='only colorize particle indices', required=False)

# for plotting contact forces
parser.add_argument('--plot_contact_force', action='store_true', help='plot nodes experiencing contact force')
parser.add_argument('--wheel_ind', type=int, help='index of the particle experiencing contact force', default=0)
parser.add_argument('--nbr_ind', type=int, help='index of neighboring particle with which particle[wheel_ind] is in contact', default=4)

# finish parsing
args = parser.parse_args()

data_dir = args.data_dir
img_dir = args.img_dir
setup_filename = args.setup_file
dotsize = args.dotsize

if args.all_dir:
    data_dir = args.all_dir
    img_dir = args.all_dir
    setup_filename = args.all_dir + '/setup.h5'

# plot info
plotinfo_file = data_dir + '/plotinfo.h5'
p = h5py.File(plotinfo_file, "r")
fc = int(p['f_l_counter'][0, 0])
lc = int(p['f_l_counter'][1, 0])
dt = float(p['dt'][0, 0])
modulo = int(p['modulo'][0, 0])

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
        ref_connectivity = np.array(f[name + '/Connectivity']).astype(int)
        N = len(f[name + '/Pos'])
        orig_tot_nbrs = total_neighbors(ref_connectivity, N)
        ref_n[name] = orig_tot_nbrs
contact_radius = np.array(f['pairwise/contact_radius']).astype(float)[0][0]



#--------------------------------------------------            
#--------------------Heart-------------------------
#-------------------------------------------------- 



def genplot(t):
    print(t)
    tc_ind = ('%05d' % t)
    filename = data_dir + '/tc_' + tc_ind + '.h5'
    wall_filename = data_dir + '/wall_' + tc_ind + '.h5'
    out_png = img_dir + '/' + args.img_prefix + '_' + tc_ind + '.png'

    # Use a clean style
    plt.style.use('seaborn-v0_8-white')

    # Set the color palette
    sns.set_palette("colorblind")

    fig = plt.figure()
    
    ax = fig.add_subplot(111, projection='3d')
 
    f = h5py.File(filename, "r")
    for name in f:
        if re.match(r'P_[0-9]+', name):
            pid = int(name[2:])
            Pos = np.array(f[name + '/CurrPos'])
            c = 'blue'
            vmin = None
            vmax = None
            
           

            if len(Pos) == 0:
                continue
           #----------------------------cross section ------------------------------ 
           # Initialize P_conn
            P_conn = np.array(f[name + '/Connectivity']).astype(int) 
            
            y_cut=5
            # Apply a cross-sectional filter
            if y_cut is not None:
                
                mask = Pos[:, 1] <= y_cut  # Keep points where y <= y_cut
                Pos = Pos[mask]
                
                # Map old indices to new indices
                valid_indices = np.where(mask)[0]
                index_map = {old: new for new, old in enumerate(valid_indices)}
                
                # Filter P_conn to exclude connections with masked-out nodes
                filtered_conn = []
                for conn in P_conn:
                    if conn[0] in index_map and conn[1] in index_map:
                        filtered_conn.append([index_map[conn[0]], index_map[conn[1]]])
                P_conn = np.array(filtered_conn)
                N = len(Pos)  # Update N to reflect the new number of nodes 

            
                # Map orig_nbrs to the filtered indices
                orig_nbrs_full = np.array(ref_n[name])  # Original neighbors for all nodes
                orig_nbrs = np.zeros(N)  # Initialize with zeros for the filtered nodes
                for old_index, new_index in index_map.items():
                    orig_nbrs[new_index] = orig_nbrs_full[old_index]

                # Calculate now_nbrs from the filtered P_conn
                now_nbrs = total_neighbors(P_conn, N)
            #-------The proportion of lost neighbors for each node in a network
                damage = (orig_nbrs - now_nbrs) / orig_nbrs
                c = damage
              
           #----------------------------cross section ------------------------------ 

            if (args.quantity == 'force'):
                force = np.array(f[name + '/force'])
                f_norm = np.sqrt(np.sum(force ** 2, axis=1))
            if args.adaptive_dotsize:
                dval = scaled_ds[pid]
            else:
                dval = dotsize 



            filter_yellow=1  
            if filter_yellow:
                print("I am filtering fire ants") 
                if isinstance(c, np.ndarray):
                    non_yellow_indices = (c < 0.9) | (c > 1.1)  # Adjust the range as needed
                    Pos = Pos[non_yellow_indices]
                    c = c[non_yellow_indices]
            scale_factor=args.scale_factor
            forcolorbar = ax.scatter(Pos[:, 0], Pos[:, 1], Pos[:, 2], c=c, edgecolors='black', s=dval *scale_factor, marker='.', linewidth=0,
                                     cmap=args.colormap, vmin=vmin, vmax=vmax, alpha=args.alpha)
    # -------------------------contact-----------------------------------
            if args.plot_contact_force:
                if pid == args.wheel_ind:
                    cnode_list = []
                    nbr_name = ('P_%05d' % args.nbr_ind)
                    print(nbr_name + '/CurrPos')  
                    nbr_currpos = np.array(f[nbr_name + '/CurrPos'])
                    print(nbr_currpos) 
                    for node in range(len(Pos)):
                        for i in range(len(nbr_currpos)):
                            dist = np.sqrt(np.sum((nbr_currpos[i] - Pos[node]) ** 2))
                            if dist <= contact_radius:
                                cnode_list.append(node)
                                break
                    plt.scatter(Pos[cnode_list, 0], Pos[cnode_list, 1], Pos[cnode_list, 2], color='red', s=dval,
                                marker='.', linewidth=0, cmap=args.colormap)

            if args.plot_contact_rad:
                if pid == args.contact_rad_ind:
                    node = args.contact_rad_node
                    x0 = Pos[node][0]
                    y0 = Pos[node][1]
                    steps = 20
                    tt = np.linspace(0, 2 * np.pi, steps, endpoint=True)
                    xx = contact_radius * np.cos(tt)
                    yy = contact_radius * np.sin(tt)
                    plt.plot(xx + x0, yy + y0)

    # -------------------------wall-----------------------------------
    if args.nowall:
        pass
    else:
        wall_alpha = 0.5
        wall_linewidth = .1

        w = h5py.File(wall_filename, "r")
        ss = np.array(w['wall_info'])[:,0]
        x_min = ss[0]
        y_min = ss[1]
        z_min = ss[2]
        x_max = ss[3]
        y_max = ss[4]
        z_max = ss[5]
        
        kalthoff=args.kalthoff
        print(kalthoff)
        if 0==2:
            x_min = ss[0] + 0.1
            y_min = ss[1] + 0.1
            z_min = ss[2] + 0.1
            x_max = ss[3] - 0.1
            y_max = ss[4]
            z_max = ss[5] - 0.1


        Z = np.array([
            [x_min, y_min, z_min],
            [x_max, y_min, z_min],
            [x_max, y_max, z_min],
            [x_min, y_max, z_min],
            [x_min, y_min, z_max],
            [x_max, y_min, z_max],
            [x_max, y_max, z_max],
            [x_min, y_max, z_max]
            ])

        # list of faces
        faces = [
             [Z[0],Z[1],Z[2],Z[3]],
             [Z[4],Z[5],Z[6],Z[7]],
             [Z[0],Z[1],Z[5],Z[4]],
             [Z[2],Z[3],Z[7],Z[6]],
             [Z[1],Z[2],Z[6],Z[5]],
             [Z[4],Z[7],Z[3],Z[0]]
             ]

        ls = np.array([
            [Z[0], Z[1]],
            [Z[1], Z[2]],
            [Z[2], Z[3]],
            [Z[3], Z[0]],

            [Z[4], Z[5]],
            [Z[5], Z[6]],
            [Z[6], Z[7]],
            [Z[7], Z[4]],

            [Z[0], Z[4]],
            [Z[1], Z[5]],
            [Z[2], Z[6]],
            [Z[3], Z[7]]
            ])

        wall_alpha = 0.2
        lc = Line3DCollection(ls, linewidths=wall_linewidth, colors='gray', alpha=wall_alpha)

        # Define the colors for specific faces
        # Only the faces listed here will be colored and plotted
        face_colors = {
            0: '#d9d9d9',  # Back face
            2: '#bdd7e7',  # Back face
            5: '#bdd7e7',  # Right face
            1: '#decbe4',  # Tope moving face
            # Note: Faces 1 (Top) and 4 (Left) are omitted and won't be plotted
        }

        wall_alpha = 0.3  # Transparency for the walls

        # Loop through each face and add it to the plot if it has a color
        for i, face in enumerate(faces):
            if i not in face_colors:
                continue  # Skip this face if not in the face_colors dictionary

            face_collection = Poly3DCollection([face], alpha=wall_alpha)
            face_collection.set_facecolor(face_colors[i])  # Use the specified color
            face_collection.set_edgecolor('k')  # Black edges for visibility
            
            ax.add_collection3d(face_collection)
        ax.add_collection(lc)

    if  args.nocolorbar:
        fig.colorbar(forcolorbar)
    if not args.nogrid:
        plt.grid()
    if args.xlim:
        xx = args.xlim.split(' ')
        ax.set_xlim3d([float(xx[0]), float(xx[1])])
    if args.ylim:
        yy = args.ylim.split(' ')
        ax.set_ylim3d([float(yy[0]), float(yy[1])])
    if args.zlim:
        zz = args.zlim.split(' ')
        ax.set_zlim3d([float(zz[0]), float(zz[1])])

    # Set plot limits and aspect ratio
    mx = np.amax(np.abs(ss))
    XYZlim = [-mx, mx]
    XYlim = [-0.27,0.27 ]
    
    Xlim = [-0.02,0.27 ]
    ax.set_xlim3d(XYZlim)
    ax.set_ylim3d(XYZlim) 
    ax.set_zlim3d(XYZlim) 
    ax.set_box_aspect((1, 1, 1))
    if kalthoff:
        XYZlim = [-mx, mx]
        # Adjust limits to zoom in on the object
        ax.set_xlim3d([-0.09, 0.2])  # Adjust these values to zoom in on your object
        ax.set_ylim3d([-0.06, 0.03])
        ax.set_zlim3d([-0.05, 0.2])
   
        ax.set_box_aspect((1, 1, 1))
    if args.motion:
        view_az = args.view_az + args.motion_az_stepsize * (t - fc)
        view_angle = args.view_angle + args.motion_angle_stepsize * (t - fc)
        ax.view_init(view_az, view_angle)
    else:
        ax.view_init(args.view_az, args.view_angle)

    #--------------------plot setup---------------
    # Thoroughly remove grid
    ax.grid(False)       # Turn off grid at the axis level
    sns.despine(ax=ax)   # Remove spines while keeping the box
    ax.set_box_aspect([1, 1, 1])  # Adjust the aspect ratio if necessary
    #plt.show() 
    plt.savefig(out_png, dpi=300, bbox_inches='tight')
    plt.close()



#--------------------------------------------------            
#--------------------Heart-------------------------
#-------------------------------------------------- 


# if args.serial:
#     ## serial j
#     for t in range(fc, lc + 1):
#         genplot(t)
# else:
#     ## parallel
#     a_pool = Pool()
#     a_pool.map(genplot, range(fc, lc + 1))
#     a_pool.close()
#
# genplot(491)
# genplot(491)
# genplot(931)
#

# genplot(81)
# genplot(181)
# genplot(381)
# genplot(810)
#
# genplot(1)
# genplot(51)
# genplot(212)
# genplot(551)
#
genplot(1)



