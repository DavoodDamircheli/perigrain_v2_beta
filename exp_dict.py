import numpy as np
import copy
import h5py
import matplotlib.pyplot as plt
# to plot a collection of lines
from matplotlib.collections import LineCollection
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection
# from gekko import GEKKO

#from multiprocessing import Pool
from pathos.multiprocessing import ProcessingPool as Pool
from itertools import combinations
# from functools import partial
import time
from random import seed
from random import random

import shape_dict
import material_dict
from genmesh import genmesh
from arrangements import get_incenter_mesh_loc
import pdb
#pdb.set_trace()
#######################################################################
class Material(object):
    """docstring for Material"""
    def __init__(self, delta, rho, snot, cnot, bulk_modulus, E = None, nu = None, Gnot = None, shear_modulus = None, name=None):
        super(Material, self).__init__()
        self.delta  = delta
        self.rho    = rho     
        self.snot   = snot
        self.cnot   = cnot
        self.bulk_modulus = bulk_modulus

        self.E = E
        self.nu = nu
        self.Gnot = Gnot
        self.shear_modulus = shear_modulus

        self.name = name

def peri_deformable(delta, rho_scale=1, K_scale=1, G_scale=1, Gnot_scale=1):
    """ Smaller fracture toughness
    """
    bulk_modulus = 2.0e9 * K_scale
    shear_modulus = 1.33e+09 * G_scale
    rho=1200.0	* rho_scale
    Gnot = 135.0 * Gnot_scale

    nu = 1/3
    E = 9 * bulk_modulus * shear_modulus / ( 9 * bulk_modulus + shear_modulus);

    cnot = 24 * E /( (1 - nu) * np.pi * (delta**3) );
    snot = np.sqrt(4 * np.pi * Gnot /(9*E*delta));

    return Material(delta, rho, snot, cnot, bulk_modulus, E = E, nu = nu, Gnot = Gnot, shear_modulus = shear_modulus )

#######################################################################


def single_bond_ext(ij_pair, pos, delta, remove_ncvx_bonds, interceptor):
    i = ij_pair[0]
    j = ij_pair[1]

    p_i = pos[i]
    p_j = pos[j]
    d = np.sqrt(np.sum(np.square(p_j - p_i))) #norm
    if (d <= delta):
        ##### remove nonconvex bonds
        if remove_ncvx_bonds:
            intersects = False
            for k in range(len(interceptor)):
                # print(k)
                A1 = interceptor[k][0:2]
                A2 = interceptor[k][2:4]
                if (lines_intersect(A1, A2, p_i, p_j)):
                    intersects = True
                    break
            if not intersects:
                return [i,j]
            # else:
                # return None
        else:
            return [i,j]
    # else:
        # return None

import numpy as np


#-----------------intersection of line segment with 3d wall-------------
def wall_line_intrsect(p1, p2, p6, p7, z_extrusion):
    """
    Check if the line segment from p6 to p7 intersects with a wall created by extruding the line segment p1-p2 in the Z direction.
    The wall is perpendicular to the X-axis and parallel to the Y and Z axes.
    """
    # X-coordinate of the wall
    wall_x = p1[0]

    # X-coordinates of the line segment
    x6, x7 = p6[0], p7[0]
    # Check if the line segment crosses the X-coordinate of the wall
    if min(x6, x7) <= wall_x <= max(x6, x7):
        # Calculate the parameter t for the intersection point
        t = (wall_x - x6) / (x7 - x6)

        # Intersection point in 3D
        intersection_3d = np.array(p6) + t * (np.array(p7) - np.array(p6))

        # Check if the intersection point's Y and Z coordinates are within the wall's bounds
        y_bounds = sorted([p1[1], p2[1]])
        z_bounds = sorted([p1[2], p1[2] + z_extrusion])
        if y_bounds[0] <= intersection_3d[1] <= y_bounds[1] and z_bounds[0] <= intersection_3d[2] <= z_bounds[1]:
           return True
        
    return False 




#--------------end intersection of line segment with 3d wall------------



def line_intersection_3d(A1, A2, A3, A4):
    # Convert points to NumPy arrays
    A1, A2, A3, A4 = np.array(A1), np.array(A2), np.array(A3), np.array(A4)

    # Direction vectors for the lines
    d1 = A2 - A1
    d2 = A4 - A3

    # The vector from A1 to A3
    r = A3 - A1

    # Cross products
    cross_d1_d2 = np.cross(d1, d2)
    cross_r_d2 = np.cross(r, d2)

    # Check if lines are parallel (cross product of direction vectors is zero)
    if np.linalg.norm(cross_d1_d2) == 0:
        # Check if lines are collinear
        if np.linalg.norm(cross_r_d2) == 0:
            # Lines are collinear - they may intersect or overlap
            return True
        else:
            # Lines are parallel but not collinear - they do not intersect
            return False

    # The scalar projection of r onto the cross product of d1 and d2
    t = np.dot(cross_r_d2, cross_d1_d2) / np.linalg.norm(cross_d1_d2)**2

    # Check if the intersection point is within the line segments
    intersection_point = A1 + t * d1
    within_segment = all(np.min([A1[i], A2[i]]) <= intersection_point[i] <= np.max([A1[i], A2[i]]) for i in range(3))

    return within_segment



def lines_intersect(A1, A2, A3, A4):
    """ Whether the lines A1-A2 and A3-A4 intersect
    :returns: TODO
    : Ai = [Ai_x, Ai_y] is an np.array
    """

    A = np.array([A1, A2, A3, A4])
    x= A[:,0]
    y= A[:,1]

    M = np.array([ [x[0] - x[1], x[3] - x[2]], [y[0] - y[1], y[3] - y[2]] ])
    b = np.array( [ x[3] - x[1], y[3] - y[1] ] )

    if (np.linalg.det(M) == 0):
        return False
    else:
        sols = np.linalg.solve(M, b)
        if ((sols[0] < 1) and (sols[0] > 0) and (sols[1] < 1) and (sols[1] > 0)):
            return True
        else:
            return False

def get_nbd_vertex_elem(T, total_vertices):
    """ From the element (facet) array, compute the facet-neighborhood array of the vertices

    :T: vertex indices of all elements
    :returns: neighboring facet indices of all vertices

    """

    full_list = []

    for i in range(total_vertices):
        # vertex i won't appear more than once in a single row
        row_list = np.where(T == i)[0]
        full_list.append(row_list.tolist())

    return full_list
    

class Particle(object):
    """ A discretized particle with nodes 
    Convention: All member variables are made accessible in the first level to make it compatible with other functions, i.e. self.pos instead of self.mesh.pos
    """
    # def __init__(self, shape, meshsize, material):
    def __init__(self, mesh, shape, material, nbdarr_in_parallel=True):
        self.shape = shape


        # generate mesh
        # mesh = genmesh(shape.P, meshsize)
        # self.mesh = mesh
        #redundant assignment.  Should be self.mesh.pos
        #Keeping it for compatibility with other functions.

        self.pos = mesh.pos
        self.vol = mesh.vol


        # boundary nodes
        self.edge_nodes = mesh.bdry_nodes

        # default nonlocal boundary nodes: all
        self.nonlocal_bdry_nodes = range(len(self.pos))

        # print('Computing boundary nodes: ')
        # print('With hard-coded contact radius: ')
        # c_R = 1e-3/5
        # temp = []
        # for p in range(len(self.edge_nodes)):
            # e_node = self.edge_nodes[p]
            # pos_node = self.pos[e_node]
            # for q in range(len(self.pos)):
                # pos_other = self.pos[q]
                # d = np.sqrt(np.sum(np.square(pos_node - pos_other))) #norm
                # # if (d <= self.contact.contact_radius):
                # if (d <= c_R*1.1):
                    # temp.append(q)
        # self.nonlocal_bdry_nodes = list(set(temp))

        self.bdry_edges = mesh.bdry_edges

        # dimension
        total_nodes = len(self.pos)
        self.dim = len(self.pos[0])
        print('dimension = ', self.dim)

        # default initial data
        self.disp = np.zeros((total_nodes, self.dim));
        self.vel = np.zeros((total_nodes, self.dim));
        self.acc = np.zeros((total_nodes, self.dim));
        self.extforce = np.zeros((total_nodes, self.dim));

        self.clamped_nodes = []

        self.material = material

        # torque info
        self.torque_axis = 2
        self.torque_val = 0

        # extra properties
        self.movable = 1
        self.breakable = 1
        self.stoppable = 1

        # neighborhood array generation
        # This is connectivity: a list of edges that are connected
        # print('Computing neighborhood connectivity.')
        self.NArr = []
        nci = self.shape.nonconvex_interceptor
        #print('nci', nci)
        self.remove_ncvx_bonds = False
        self.interceptor = []
        self.remove_ncvx_bonds_3d = False 
        if (nci and self.dim ==2):
            self.remove_ncvx_bonds = True
            self.interceptor = nci.all()
            # print('(Will remove nonconvex bonds)')
        if (nci and self.dim ==3):
            self.remove_ncvx_bonds_3d = True
            self.interceptor = nci.all()
        print("-----------------------this is nonconcex-intercetor nci----------")
        #print(nci.all())
        #----------------this is just for kalthoff-winfler-------------V
        #----------------this is just for kalthoff-winfler-------------V
        #----------------this is just for kalthoff-winfler-------------V
        #----------------this is just for kalthoff-winfler-------------V
        #----------------this is just for kalthoff-winfler-------------V
        #----------------this is just for kalthoff-winfler-------------V
        r = self.interceptor
        self.wallsCoordinate = r 
        #print("-----------------------this is nonconcex-intercetor----------")
        '''
        #print(r)
        # Separate 'result' back into the original lists
        Al1raw = [tuple(r[i:i+3]) for i in range(0, len(r), 12)]
        Al2raw = [tuple(r[i:i+3]) for i in range(3, len(r), 12)]
        Ar1raw = [tuple(r[i:i+3]) for i in range(6, len(r), 12)]
        Ar2raw = [tuple(r[i:i+3]) for i in range(9, len(r), 12)]
        #print("this is AL")
        # remove all None
        self.Al1 = [i for i in Al1raw if i !=( None,None,None)]
        self.Al2 = [i for i in Al2raw if i !=( None,None,None)]
        self.Ar1 = [i for i in Ar1raw if i !=( None,None,None)]
        self.Ar2 = [i for i in Ar2raw if i !=( None,None,None)]
        

        #plot_list_3d(Ar1raw)
        #plot_list_3d(Ar1raw)
        # start = time.time()
        #nbdarr_in_parallel=False
        debug = False 
        if debug:
            fig =plt.figure()
            ax = fig.add_subplot(111,projection='3d')
            x0,y0,z0 = zip(*self.Al1) 
            x1,y1,z1 = zip(*self.Al2) 
            
            removed_nodes =  np.array([point for point in self.pos if not any(np.array_equal(point, p) for p in self.Al1)])
            removed_nodes =  np.array([point for point in removed_nodes if not any(np.array_equal(point, p) for p in self.Al2)])

            x5,y5,z5 = zip(*removed_nodes)


            ax.scatter(x0,y0,z0,color='red',marker = '*',s=2)
            ax.scatter(x1,y1,z1,color='red',marker = '*',s=2)
            

            ax.scatter(x5,y5,z5,color='blue',marker = '*',s=2)
            plt.show()
        '''
        #nbdarr_in_parallel=False
        if nbdarr_in_parallel:
            print('nbdarr in parallel')
            ## parallel attempt
            a_pool = Pool()
            # all_bonds = a_pool.map(oneparam_f_bond, combinations(range(total_nodes),2)) 
            if (self.dim==3):
                #all_bonds = a_pool.map(self.kalthoff_single_bond_3d, combinations(range(total_nodes),2)) 
                all_bonds = a_pool.map(self.single_bond_3d_alpha, combinations(range(total_nodes),2)) 
                #all_bonds = a_pool.map(self.single_bond_3d, combinations(range(total_nodes),2)) 
            else:
                all_bonds = a_pool.map(self.single_bond, combinations(range(total_nodes),2)) 
                #all_bonds = a_pool.map(self.single_bond, combinations(range(total_nodes),2)) 
            # print(all_bonds)
        else:
            ## Serial version
            print('nbdarr in serial')
            # for ij in list(combinations(range(total_nodes), 2)):
            # print(out_bonds(self, ij))

            all_bonds = []
            for ij in list(combinations(range(total_nodes), 2)):
                if(self.dim==3):
                    #print(self.kalthoff_single_bond_3d(ij))
                    #all_bonds.append( self.kalthoff_single_bond_3d(ij))
                    
                    all_bonds.append( self.single_bond_3d_alpha(ij))
                    #all_bonds.append( self.single_bond_3d(ij))
                    #print(all_bonds)
                else:
                    all_bonds.append( self.single_bond(ij))
                    #all_bonds.append( self.single_bond(ij))

        # remove all None
        self.NArr = np.array([i for i in all_bonds if i is not None])
        print("\n  first ten element of bonds----> \n")
        print(all_bonds[:10])
        print("\n total number of all bonds in this particle are  ::\n")
        print('_'+str(len(all_bonds)), end=' ', flush=True )
  
        # print the total number of bonds in a particle
        print("\ntotal number of bonds after removing outsider bonds  ::\n")
        print('_'+str(len(self.NArr)), end=' ', flush=True)
        #print(all_bonds.count())
        # print('Done generating nbdarr')
        
        
            # ## Serial version

                # for j in range(i+1, total_nodes):
                    # # if (j != i):
                        # p_j = self.pos[j]
                        # d = np.sqrt(np.sum(np.square(p_j - p_i))) #norm
                        # if (d <= self.material.delta):
                            # ##### remove nonconvex bonds
                            # if remove_ncvx_bonds:
                                # intersects = False
                                # for k in range(len(interceptor)):
                                    # # print(k)
                                    # A1 = interceptor[k][0:2]
                                    # A2 = interceptor[k][2:4]
                                    # if (lines_intersect(A1, A2, p_i, p_j)):
                                        # intersects = True
                                        # break

                                # if not intersects:
                                    # self.NArr.append([i, j])
            # self.NArr  = np.array(self.NArr)

            # ## remove nonconvex bonds
            # if self.dim==2:
                # # print(self.shape.nonconvex_interceptor)
                # if self.shape.nonconvex_interceptor is not None:
                    # interceptor = self.shape.nonconvex_interceptor.all()
                    # print('Removing nonconvex bonds')
                    # to_delete = []
                    # for j in range(len(interceptor)):
                        # A1 = interceptor[j][0:2]
                        # A2 = interceptor[j][2:4]
                        # for i in range(len(self.NArr)):
                            # edge = self.NArr[i]
                            # A3 = self.pos[edge[0]]
                            # A4 = self.pos[edge[1]]

                            # if lines_intersect(A1, A2, A3, A4):
                                # # print('non-convex bond: ', edge)
                                # to_delete.append(i)
                                # # print(len(self.NArr))
                
                    # # delete at once. delete is not in-place
                    # self.NArr = np.delete(self.NArr, to_delete, 0)
                    # pass
            # else:
                # pass
                # # 3d code for nonconvex bond removal goes here


        # print(self.NArr)
        # print('time taken ', time.time() - start)

    def test_single_bond(self, ij_pair):
        # print('ij', ij_pair, end='\n', flush=True)
        i = ij_pair[0]
        j = ij_pair[1]

        p_i = self.pos[i]
        p_j = self.pos[j]
        d = np.sqrt(np.sum(np.square(p_j - p_i))) #norm
        if (d <= self.material.delta):
                return [i,j]
        # else:
            # return None


    def single_bond(self, ij_pair):
        # print('ij', ij_pair, end='\n', flush=True)
        i = ij_pair[0]
        j = ij_pair[1]

        p_i = self.pos[i]
        p_j = self.pos[j]
        d = np.sqrt(np.sum(np.square(p_j - p_i))) #norm
        # material.delta is peridynamic horizon
        if (d <= self.material.delta):
            ##### remove nonconvex bonds
            if self.remove_ncvx_bonds:
                intersects = False
                for k in range(len(self.interceptor)):
                    # print(k)
                    A1 = self.interceptor[k][0:2]
                    A2 = self.interceptor[k][2:4]
                    if (lines_intersect(A1, A2, p_i, p_j)):
                        intersects = True
                        break
                if not intersects:
                    return [i,j]
                # else:
                    # return None
            else:
                return [i,j]



      #--------------eliminiating bonds in 3d--------------
    def kalthoff_single_bond_3d(self, ij_pair):
        # print('ij', ij_pair, end='\n', flush=True)
        i = ij_pair[0]
        j = ij_pair[1]

        p_i = self.pos[i]
        p_j = self.pos[j]
        d = np.sqrt(np.sum(np.square(p_j - p_i))) #norm
        r = self.wallsCoordinate

                #print("-----------------------this is nonconcex-intercetor----------")
        w1A = r[0]
        w1B = r[1]  


        w2A = r[2]
        w2B = r[3]
        
        z = r[4][2]
        # material.delta is peridynamic horizon
        #print(d<=self.material.delta)
        if (d <= self.material.delta):
            ##### remove nonconvex bonds
            if self.remove_ncvx_bonds_3d: 

                #print("remove is on") 
                intersects = True 
                    
                    #-- check the first notches
                #print(wall_line_intrsect(w1A,w1B,p_i,p_j, z))
                print(wall_line_intrsect(w2A,w2B,p_i,p_j, z))
                if wall_line_intrsect(w1A,w1B,p_i,p_j, z):
                   intersects = False 
                   #print('detected')
                    #-- check the second notches
                    #print(line_intersection_3d(B1, B2, p_i, p_j))
                if  wall_line_intrsect(w2A,w2B,p_i,p_j, z):
                    intersects = False
                     
                    #print('detected')

                if  intersects:
                    return [i,j]
                else: 
                    return None
            else:
                return [i,j]




     #--------------eliminiating bonds in 3d alpha--------------
    def single_bond_3d_alpha(self, ij_pair):
        # print('ij', ij_pair, end='\n', flush=True)
        i = ij_pair[0]
        j = ij_pair[1]

        p_i = self.pos[i]
        p_j = self.pos[j]
        d = np.sqrt(np.sum(np.square(p_j - p_i))) #norm
        
        if self.remove_ncvx_bonds_3d: 

            r = self.interceptor
                    #print("-----------------------this is nonconcex-intercetor----------")
            #print(r)
            w1A = r[0]
            w1B = r[1]  


            w2A = r[2]
            w2B = r[3]
            
            z = r[4][2]

        #print(d<=self.material.delta)
        if (d <= self.material.delta):
            ##### remove nonconvex bonds
            if self.remove_ncvx_bonds_3d: 

                #print("remove is on") 
                intersects = False
                for k in range(1):
                    #print(k)
                    #print("this is A1")
                    #print(A1)
                    #plot_single_bonds(A1,A2,p_i,p_j) 
                    #print('line_intersection_3d(A1, A2, p_i, p_j)-------->``')
                    #print(line_intersection_3d(A1, A2, p_i, p_j))
                    
                    if (wall_line_intrsect(w1A,w1B,p_i,p_j, z)):
                        #print('detected w1')
                    #if (line_intersection_3d(A1, A2, p_i, p_j)):
                        #print(wall_line_intrsect(w1A,w1B,p_i,p_j, z))
                        intersects = True
                        break
                    #print('line_intersection_3d(B1, B2, p_i, p_j)-------->``')
                    #print(line_intersection_3d(B1, B2, p_i, p_j))
                    

                    if  (wall_line_intrsect(w2A,w2B,p_i,p_j, z)):
                    #if (line_intersection_3d(B1, B2, p_i, p_j)):
                        #print('detected w2')
                        #print("shit")
                        #print(wall_line_intrsect(w2A,w2B,p_i,p_j, z))
                        intersects = True
                        break

                if not intersects:
                    #print("flaseeeee")
                    return [i,j]
                # else:
                    # return None
            else:
                return [i,j]
        # else:
            # return None




     #--------------eliminiating bonds in 3d--------------
    def single_bond_3d(self, ij_pair):
        # print('ij', ij_pair, end='\n', flush=True)
        i = ij_pair[0]
        j = ij_pair[1]

        p_i = self.pos[i]
        p_j = self.pos[j]
        d = np.sqrt(np.sum(np.square(p_j - p_i))) #norm
        r = self.interceptor
                #print("-----------------------this is nonconcex-intercetor----------")
        Al1=self.Al1
        Al2=self.Al2
        Ar1=self.Ar1
        Ar2=self.Ar2

        # material.delta is peridynamic horizon
        tnodes = min(len(Al1),len(Al2),len(Ar1),len(Ar2))
        #print(d)
        #print("this is delta")
        #print(self.material.delta)
        #print(d<=self.material.delta)
        if (d <= self.material.delta):
            ##### remove nonconvex bonds
            if self.remove_ncvx_bonds_3d: 

                #print("remove is on") 
                intersects = False
                for k in range(tnodes):
                    #print(k)
                    A1 = Al1[k][0:3]
                    A2 = Al2[k][0:3]
                    #print("this is A1")
                    #print(A1)
                    #plot_single_bonds(A1,A2,p_i,p_j) 
                    #print('line_intersection_3d(A1, A2, p_i, p_j)-------->``')
                    #print(line_intersection_3d(A1, A2, p_i, p_j))
                    if (line_intersection_3d(A1, A2, p_i, p_j)):
                        intersects = True
                        break
                    #print("this is problematic K")
                    #print(k)
                    B1 = Ar1[k][0:3]
                    B2 = Ar2[k][0:3]
                    #plot_single_bonds(B1,B2,p_i,p_j) 
                    #print('line_intersection_3d(B1, B2, p_i, p_j)-------->``')
                    #print(line_intersection_3d(B1, B2, p_i, p_j))
                    if (line_intersection_3d(B1, B2, p_i, p_j)):
                        intersects = True
                        break

                if not intersects:
                    return [i,j]
                # else:
                    # return None
            else:
                return [i,j]
        # else:
            # return None

        # one parameter function to be used by pool
        # def oneparam_f_bond(ij):
        # return self.single_bond(self,ij)
        ## Using an external function to call for parallelization, otherwise, I get error 
        # return single_bond_ext(ij, self.pos, self.material.delta, remove_ncvx_bonds, interceptor)
        # print('Will generate nbdarr')
                    
    def gen_nonlocal_boundry_nodes(self, contact_radius, scaled = 1):
        """ Generate nonlocal boundary nodes
        :contact_radius: 
        """
        # print('Computing nonlocal boundary nodes: ', end = ' ', flush = True)
        if contact_radius is not None:
            temp = []
            for p in range(len(self.edge_nodes)):
                e_node = self.edge_nodes[p]
                pos_node = self.pos[e_node]
                for q in range(len(self.pos)):
                    pos_other = self.pos[q]
                    d = np.sqrt(np.sum(np.square(pos_node - pos_other))) #norm
                    # if (d <= self.contact.contact_radius):
                    if (d <= contact_radius*scaled):
                        temp.append(q)
            self.nonlocal_bdry_nodes = list(set(temp))
        # print(len(self.nonlocal_bdry_nodes))


    def rotate(self, angle):
        """ Rotate by angle in the clockwise direction
        Only for 2D particles
        """
        if self.dim!=2:
            print('Wrong dimension for rotation. Need 2. Use rotate3d for 3D particle.')
        Rot = np.array([[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)] ])
        self.pos = self.pos @ Rot   # @ is matrix multiplication, so is a.dot(b) 

    def rotate3d(self, axis, angle):
        """ Rotate by angle in the counterclockwise direction about an axis
        Only for 3d particles
        """
        if axis=='x':
            Rot = np.array([[1,0,0], [0, np.cos(angle), -np.sin(angle)], [0, np.sin(angle), np.cos(angle)] ])
        if axis=='y':
            Rot = np.array([[np.cos(angle),0, -np.sin(angle)], [0, 1, 0], [np.sin(angle), 0, np.cos(angle)] ])
        if axis=='z':
            Rot = np.array([[np.cos(angle), -np.sin(angle), 0], [np.sin(angle), np.cos(angle), 0], [0, 0, 1] ])

        self.pos = self.pos @ Rot   # @ is matrix multiplication, so is a.dot(b) 

    def trasnform_bymat(self, matrix):
        """ Transform using a matrix
        """
        # Rot = np.array([[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)] ])
        self.pos = self.pos @ matrix   # @ is matrix multiplication, so is a.dot(b) 
        print('Caution: did not transform the volume to maintain the density.')

    def shift(self, shift):
        a = self.pos
        b = shift
        #print("this is pos",a)
        #print('this is shift',shift)
        self.pos += shift

    def scale(self, scale):
        self.pos *= scale
        # update the volume
        self.vol *= (scale**self.dim)

        # delta and cnot are modified since the neighbors are not being recomputed

class GridDist(object):
    """ Particle location distribution on a grid
    """
    def __init__(self, left, right, top, bottom, nx, ny):
        self.left    = left   
        self.right   = right  
        self.top     = top    
        self.bottom  = bottom 

        self.nx = nx
        self.ny = ny

    def maximalrad(self):
        lx = (self.right - self.left)/self.nx/2
        ly = (self.top - self.bottom)/self.ny/2
        print(lx,ly)
        return min(lx, ly)

    def gridloc(self):
        """TODO: Docstring for gridloc.
        :arg1: TODO
        :returns: TODO
        """
        lx = (self.right - self.left)/self.nx
        ly = (self.top - self.bottom)/self.ny
        xx = np.linspace(self.left + lx/2, self.right-lx/2, num=self.nx, endpoint=True)
        yy = np.linspace(self.bottom + ly/2, self.top-ly/2, num=self.ny, endpoint=True)
        xv, yv = np.meshgrid(xx, yy)
        loc_list = np.array([xv.flatten(), yv.flatten()]).transpose()
        return np.array(loc_list)


class ShapeList(object):
    """Contains a list of instances of Shape class, and other info """
    def __init__(self):
        self.shape_list = []
        self.count_list = []
        self.meshsize_list = []
        self.material_list = []

        # self.msh_file_list = []

    def append(self, shape, count, meshsize, material, plot_shape = True):
    # def append(self, shape, count, meshsize, material, msh_file= None):
        """TODO: Docstring for append.

        :shape: TODO
        :count: TODO
        :meshsize: TODO
        :material: TODO
        :returns: TODO

        """
        self.shape_list.append(shape)
        self.count_list.append(count)
        self.meshsize_list.append(meshsize)
        self.material_list.append(material)



    def generate_mesh(self, dimension = 2, contact_radius = None, plot_mesh = True, plot_node_text=False, plot_shape = True, shapes_in_parallel=False, print_nnodes=True):
        """Returns a rank-2 array of particle meshes
        : shapes_in_parallel: if set to true, it computes the particle properties in parallel for each shape. Otherwise, the NbdArr and boundary node computation happens in parallel for each node. 
        Outermost parallelization is most preferable. So, for multiple shapes/particles with small number of nodes requires this option to be true. On the other hands, for single particle with many nodes, setting this to false works best.
        :returns: particles[sh][i] for shape sh count i
        """

        # if the shapes are computed in parallel, do not compute NbdArr etc in parallel
        # outermost parallelism
        #nbdarr_in_parallel  = not shapes_in_parallel

        nbdarr_in_parallel  = not shapes_in_parallel
        print('Shape:', end=' ', flush=True)

        # Helper function
        def gen_particle(sh):
            # print(sh, end=' ', flush=True)
            shape = self.shape_list[sh]
            if shape.msh_file is None:
                mesh = genmesh(P_bdry=shape.P, pygmsh_geom=shape.pygmsh_geom, meshsize=self.meshsize_list[sh], dimension = dimension)

                # print(shape.P)
                if (plot_shape and (shape.P is not None)):
                    shape.plot(bdry_arrow=True, extended_bdry=True, angle_bisector=True)
            else:
                # when the .msh is specified
                # mesh = genmesh(P_bdry=None, meshsize=None, msh_file=shape.msh_file, dimension = dimension, do_plot = plot_mesh)
                mesh = genmesh(P_bdry=None, meshsize=None, msh_file=shape.msh_file, dimension = dimension)
                if plot_mesh:
                    mesh.plot(dotsize=22, plot_node_text=plot_node_text, highlight_bdry_nodes=True)



            # print('total mesh volume: ', np.sum(mesh.vol))
            
            PP = Particle(mesh=mesh, shape=self.shape_list[sh], material=self.material_list[sh], nbdarr_in_parallel=nbdarr_in_parallel)
            
            
            nnodes = len(PP.pos)
            print("\n The number of pos in shape "+str(sh)+" are:  "+'('+str(nnodes)+')', end=' ', flush=True)

            # print('Done generating particles')

            # generate nonlocal boundary nodes
            PP.gen_nonlocal_boundry_nodes(contact_radius)


            particles_sh = []
            for i in range(self.count_list[sh]):
                # shallow copy (=) refers to the same object, instead of creating a new one
                this_PP = copy.deepcopy(PP)
                particles_sh.append(this_PP)

            # particles.append(particles_sh)
            return particles_sh


        start_sh = time.time()
        if shapes_in_parallel:
            print('Shapes in parallel')
            # parallel attempt
            shape_pool = Pool()
            particles = shape_pool.map(gen_particle, range(len(self.shape_list))) 
        else:
            print('Shapes in serial')
            # serial version
            particles = []
            print("shapeList is",len(self.shape_list)) 
            for sh in range(len(self.shape_list)):
                particles.append(gen_particle(sh))

        print('\n')
        print('\n time taken to generate all shapes\n', time.time() - start_sh)
        
        return particles

        
# def generate_particles(shapes, shape_count, meshsize, material_list):
    # """ Generate Particle related to each shape
    # :shapes: a list of instances of the shape class generated from shape_dict.py
    # :shape_count: how many members of each shape to create
    # :returns: 2d array of particles, particles[shape][count]
    # """
    # particles = []
    # for sh in range(len(shapes)):
        # print('Shape', sh)
        # PP = Particle(shapes[sh], meshsize, material_list[sh])

        # particles_sh = []
        # for i in range(shape_count[sh]):
            # # shallow copy (=) refers to the same object, instead of creating a new one
            # this_PP = copy.deepcopy(PP)
            # particles_sh.append(this_PP)

        # particles.append(particles_sh)
    
    # return particles


class Wall(object):
    """ Geometric wall """
    def __init__(self, allow = 0, left = None, right = None, top = None, bottom = None ):
        self.allow  = allow
        self.left        = left
        self.right       = right
        self.top         = top 
        self.bottom      = bottom

        self.reaction = None

    def get_lrtp(self):
        """ Returns the wall dimensions in an array
        :returns: TODO

        """
        return np.array([self.left, self.right, self.top, self.bottom])

    def set_size(self, ss):
        """ Reads the wall size from an array
        """
        if (len(ss) !=4):
            print('Error: input array dimension is not 4.')
        self.left        = ss[0]
        self.right       = ss[1]
        self.top         = ss[2]
        self.bottom      = ss[3]

    def get_lines(self):
        """ List of lines that are wall boundaries
        To feed to the collection
        """
        a = [self.left, self.bottom]
        b = [self.right, self.bottom]
        c = [self.right, self.top]
        d = [self.left, self.top]

        return [ [a, b], [b,c], [c,d], [d,a] ]

    def print(self):
        print('allow =  ', self.allow)
        print('left: ', self.left)
        print('right: ', self.right)
        print('top: ', self.top)
        print('bottom: ', self.bottom)

    def get_h(self):
        """get the horizontal length of the wall
        """
        return self.right - self.left
    def get_v(self):
        """get the vertical height of the wall
        """
        return self.top - self.bottom
    def wall_vol(self):
        """returns the wall volume
        """
        return self.get_v() * self.get_h()

class Wall3d(object):
    """ Geometric wall """
    def __init__(self, allow = 0, x_min = None, y_min = None, z_min = None, x_max = None, y_max = None, z_max = None):
        self.allow  = allow
        self.x_min = x_min
        self.y_min = y_min
        self.z_min = z_min
        self.x_max = x_max
        self.y_max = y_max
        self.z_max = z_max

    def get_lrtp(self):
        """ Returns the wall dimensions in an array
        :returns: TODO

        """
        return np.array([self.x_min, self.y_min, self.z_min, self.x_max, self.y_max, self.z_max])

    def set_size(self, ss):
        """ Reads the wall size from an array
        """
        if (len(ss) !=6):
            print('Error: input array dimension is not 4.')

        self.x_min = ss[0] 
        self.y_min = ss[1]
        self.z_min = ss[2]
        self.x_max = ss[3]
        self.y_max = ss[4]
        self.z_max = ss[5]





    def get_faces(self):

        # list of vertices
        Z = np.array([
            [self.x_min, self.y_min, self.z_min],
            [self.x_max, self.y_min, self.z_min],
            [self.x_max, self.y_max, self.z_min],
            [self.x_min, self.y_max, self.z_min],
            [self.x_min, self.y_min, self.z_max],
            [self.x_max, self.y_min, self.z_max],
            [self.x_max, self.y_max, self.z_max],
            [self.x_min, self.y_max, self.z_max]
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

        return faces
    
    def get_lines(self):
        # list of vertices
        Z = np.array([
            [self.x_min, self.y_min, self.z_min],
            [self.x_max, self.y_min, self.z_min],
            [self.x_max, self.y_max, self.z_min],
            [self.x_min, self.y_max, self.z_min],
            [self.x_min, self.y_min, self.z_max],
            [self.x_max, self.y_min, self.z_max],
            [self.x_max, self.y_max, self.z_max],
            [self.x_min, self.y_max, self.z_max]
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

        lines = np.array([
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

        return lines

    def print(self):
        print('allow =  ', self.allow)
        print('x_min: ', self.x_min)
        print('y_min: ', self.y_min)
        print('z_min: ', self.z_min)
        print('x_max: ', self.x_max)
        print('y_max: ', self.y_max)
        print('z_max: ', self.z_max)

class Contact(object):
    """Pairwise contact properties"""
    def __init__(self, contact_radius, normal_stiffness, damping_ratio, friction_coefficient):
        self.contact_radius         = contact_radius
        self.normal_stiffness        = normal_stiffness
        self.damping_ratio           = damping_ratio
        self.friction_coefficient    = friction_coefficient 

    def print(self):
        """print info

        :f: TODO
        :returns: TODO

        """
        print('contact_radius: ', self.contact_radius)
        print('normal_stiffness: ', self.normal_stiffness)
        print('damping_ratio: ', self.damping_ratio)
        print('friction_coefficient: ', self.friction_coefficient)

#---------------------Plot Bonds---------------------
def plot_single_bonds(A1,A2,B1,B2):
            
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Plotting lines
    #?????????????????????????????????/in qalat nist???????????????? start and end
    ax.plot([A1[0], A2[0]], [A1[1], A2[1]], [A1[2], A2[2]])
    ax.plot([B1[0], B2[0]], [B1[1], B2[1]], [B1[2], B2[2]])
    # Setting labels
    plt.title('bonds after removing nonconvex part')
    ax.set_xlabel('X axis')
    ax.set_ylabel('Y axis')
    ax.set_zlabel('Z axis')

    # Show plot
    plt.show()


def plot_bonds(particles):
    

    totalbonds = []
    total_convex_bond = []
    for sh in range(len(particles)):
        shape = particles[sh]
        for count in range(len(shape)):
            
            part = particles[sh][count]
            #nnodes_vec.append(len(part.pos))
            total_nodes = len(part.pos)
            #total_convex_bond =[(part.pos[a],part.pos[b]) for a,b in list(combinations(range(total_nodes), 2))] 
            total_convex_bond =[(part.pos[a],part.pos[b]) for a,b in part.NArr] 
            # Creating a 3D plot
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')

            # Plotting lines
            #?????????????????????????????????/in qalat nist???????????????? start and end
            for start, end in total_convex_bond:  
                ax.plot([start[0], end[0]], [start[1], end[1]], [start[2], end[2]])
            # Setting labels
            plt.title('bonds after removing nonconvex part')
            ax.set_xlabel('X axis')
            ax.set_ylabel('Y axis')
            ax.set_zlabel('Z axis')
            # save the image
            plt.savefig('bonds.png', dpi=300, bbox_inches='tight')
 
            # Show plot
            plt.show()
def plot_list_3d(A):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    print("no bond outside!!!")
    x2,y2,z2= zip(*A) 

    #Plotting the points from list1 in red and list2 in blue

    ax.scatter(x2, y2, z2, color='blue', label='concexNodes')
   
    plt.title('plot nodes')
    ax.set_xlabel('X axis')
    ax.set_ylabel('Y axis')
    ax.set_zlabel('Z axis')
    # save the image

    # Show plot
    plt.show() 

def plot_bonds_nonConvex(particles):
    

    totalbonds = []
    total_convex_bond = []
    for sh in range(len(particles)):
        shape = particles[sh]
        for count in range(len(shape)):
            
            part = particles[sh][count]
 
            #------------------------
            total_convex_node_Idx  = list({element for sublist in part.NArr for element in sublist})
             
            total_nodes = part.pos
            convx_nodes =[part.pos[a] for a in total_convex_node_Idx] 
            removed_nodes =  np.array([point for point in total_nodes if not any(np.array_equal(point, p) for p in convx_nodes)])
 
            total_Nonconvex_bond =[(part.pos[a],part.pos[b]) for a,b in removed_nodes] 
            #----------------------
            # Creating a 3D plot

            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')

            # Plotting lines
            #?????????????????????????????????/in qalat nist???????????????? start and end
            for start, end in total_Nonconvex_bond:  
                ax.plot([start[0], end[0]], [start[1], end[1]], [start[2], end[2]])
            # Setting labels
            plt.title('bonds after removing nonconvex part')
            ax.set_xlabel('X axis')
            ax.set_ylabel('Y axis')
            ax.set_zlabel('Z axis')
            # save the image
            plt.savefig('bondsNon.png', dpi=300, bbox_inches='tight')
 
            # Show plot
            plt.show()


#--------------------End Plot Bonds-------------------
def plot_nodes(particles):
    

    totalbonds = []
    total_convex_bond = []
    #for sh in range(len(particles)):
    
    for sh in [0]: 
        shape = particles[sh]
        for count in range(len(shape)):
            
            part = particles[sh][count]
            #nnodes_vec.append(len(part.pos))
            #x,y,z = zip(*part.pos) 
            total_convex_node_Idx  = list({element for sublist in part.NArr for element in sublist})
             
            total_nodes = part.pos
            convx_nodes =[part.pos[a] for a in total_convex_node_Idx] 
            #print(part.pos)
            #print(par.NArr)
            #removed_nodes = [point for point in total_nodes if point not in convx_nodes]   
            removed_nodes =  np.array([point for point in total_nodes if not any(np.array_equal(point, p) for p in convx_nodes)])
            #print(removed_nodes)
            #print(total_convex_node_Idx)
            #print(total_nodes)
                        # Creating a 3D plot
            print("\n number of removed\n"+ str(len(removed_nodes)))
            print("\n total nodes are \n"+ str(len(total_nodes)))
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            print("\n size removed nodes are \n ")
            print(len(removed_nodes))
            if len(removed_nodes)<1:
                print("\n no bond outside!!!\n")
                x2,y2,z2= zip(*convx_nodes) 
            else:
                print("\n ha we removed some bonds\n")
                x2,y2,z2= zip(*convx_nodes) 
                x1,y1,z1= zip(*removed_nodes) 
                ax.scatter(x1, y1, z1, color='green', label='totalNods',s=2)

            #Plotting the points from list1 in red and list2 in blue

            ax.scatter(x2, y2, z2, color='red', label='concexNodes',s=2)
           
            plt.title('plot nodes')
            ax.set_xlabel('X axis')
            ax.set_ylabel('Y axis')
            ax.set_zlabel('Z axis')
            # save the image
            plt.savefig('nodes.png', dpi=300, bbox_inches='tight')
 
            # Show plot
            plt.show() 
def plot_node_bond_distribution(particles, nnodes_filename='nnodes_distribution.png', nbonds_filename='nbonds_distribution.png', bins=None):
    nnodes_vec = []
    nbonds_vec = []
    for sh in range(len(particles)):
        shape = particles[sh]
        for count in range(len(shape)):
            part = particles[sh][count]
            nnodes_vec.append(len(part.pos))
            nbonds_vec.append(len(part.NArr))
    print(len(part.NArr))
    plt.hist(nnodes_vec, bins=bins)
    plt.title('Number of nodes: distribution')
    plt.xlabel('nnodes')
    plt.ylabel('frequency')
    # save the image
    plt.savefig(nnodes_filename, dpi=300, bbox_inches='tight')
    plt.close()
    
    plt.hist(nbonds_vec, bins=bins)
    plt.title('Number of bonds: distribution')
    plt.xlabel('nbonds')
    plt.ylabel('frequency')
    # save the image
    plt.savefig(nbonds_filename, dpi=300, bbox_inches='tight')
    plt.close()

def plot_setup(particles, dotsize=0.5, contact_radius=None, delta=None, linewidth=1, wall=None, show_plot=True, adaptive_dotsize=False, show_particle_index=False, save_filename='setup.png'):
    """plots the particle setup
    :particles: array of particles
    """

    dim = len(particles[0][0].pos[0])

    min_vol = 1
    if adaptive_dotsize:
        slot_shape = np.zeros(len(particles))
        for sh in range(len(particles)):
            shape = particles[sh]
            slot_count = np.zeros(len(shape))
            for count in range(len(shape)):
                part = particles[sh][count]
                slot_count[count] = np.min(part.vol)
            slot_shape[sh] = np.min(slot_count)
        min_vol = np.min(slot_shape)

    # print('min vol', min_vol)

    for sh in range(len(particles)):
        shape = particles[sh]
        for count in range(len(shape)):
            part = particles[sh][count]
            P = part.pos
            if adaptive_dotsize:
                dval = np.sqrt(part.vol/min_vol) * float(dotsize)
            else:
                dval = dotsize

            plt.scatter(P[:,0], P[:,1], s=dval, linewidths=0)

            # show clamped node
            cc = part.clamped_nodes
            for j in range(len(cc)):
                c = cc[j]
                plt.scatter(P[c,0], P[c,1], c='r', s=dval, linewidths=0)

            if show_particle_index:
                # # text
                plt.annotate(str(count), np.mean(P, axis = 0))

    if contact_radius:
        part = particles[0][0]
        node = 0
        x0 = part.pos[node][0]
        y0 = part.pos[node][1]
        tt = np.linspace(0, 2*np.pi, endpoint=True)
        xx = contact_radius * np.cos(tt)
        yy = contact_radius * np.sin(tt)
        plt.plot( xx + x0, yy + y0, linewidth=linewidth)

    if delta:
        part = particles[0][0]
        node = 0
        x0 = part.pos[node][0]
        y0 = part.pos[node][1]
        tt = np.linspace(0, 2*np.pi, endpoint=True)
        xx = delta * np.cos(tt)
        yy = delta * np.sin(tt)
        plt.plot( xx + x0, yy + y0, 'r', linewidth=linewidth)

    if wall:
        wi = wall.get_lrtp()
        a = [wi[0],  wi[3]]
        b = [wi[1], wi[3]]
        c = [wi[1], wi[2]]
        d = [wi[0],  wi[2]]
        ls = [ [a, b], [b,c], [c,d], [d,a] ]
        lc = LineCollection(ls, linewidths=1, colors='b')
        plt.gca().add_collection(lc)

    plt.axis('scaled')
                
    plt.grid()
    print('Saving experiment setup png to', save_filename)
    plt.savefig(save_filename, dpi=300, bbox_inches='tight')
    if show_plot:
        plt.show()

def plot3d_setup(particles, dotsize=0.5, contact_radius=None, delta=None, wall=None, show_plot=True, adaptive_dotsize=False, show_particle_index=False, save_filename='setup.png'):
    """plots the particle setup
    :particles: array of particles
    """

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    min_vol = 1
    if adaptive_dotsize:
        slot_shape = np.zeros(len(particles))
        for sh in range(len(particles)):
            shape = particles[sh]
            slot_count = np.zeros(len(shape))
            for count in range(len(shape)):
                part = particles[sh][count]
                slot_count[count] = np.min(part.vol)
            slot_shape[sh] = np.min(slot_count)
        min_vol = np.min(slot_shape)

    for sh in range(len(particles)):
        shape = particles[sh]
        for count in range(len(shape)):
            part = particles[sh][count]
            P = part.pos
            if adaptive_dotsize:
                dval = np.sqrt(part.vol/min_vol) * float(dotsize)
            else:
                dval = dotsize
            # check the color of the this part????????????????????????????
            ax.scatter(P[:,0], P[:,1], P[:,2], s=dval, marker = '.', linewidth=0, cmap='viridis')

            # check the color of the this part????????????????????????????
            # show clamped node
            cc = part.clamped_nodes
            for j in range(len(cc)):
                c = cc[j]
                ax.scatter(P[c,0], P[c,1], P[c,2], c='r', s=dval, linewidths=0)

            if show_particle_index:
                meanP = np.mean(P, axis=0)
                ax.text(meanP[0], meanP[1], meanP[2], str(count))

    if contact_radius:
        part = particles[0][0]
        node = 0
        x0 = part.pos[node][0]
        y0 = part.pos[node][1]
        z0 = part.pos[node][2]
        tt = np.linspace(0, 2*np.pi, endpoint=True)
        xx = contact_radius * np.cos(tt)
        yy = contact_radius * np.sin(tt)
        plt.plot( xx + x0, yy + y0, z0)

    if delta:
        part = particles[0][0]
        node = 0
        x0 = part.pos[node][0]
        y0 = part.pos[node][1]
        z0 = part.pos[node][2]
        tt = np.linspace(0, 2*np.pi, endpoint=True)
        xx = delta * np.cos(tt)
        yy = delta * np.sin(tt)
        plt.plot( xx + x0, yy + y0, z0, 'r')

    if wall:
        wall_alpha = 0.8
        wall_linewidth = 1
        ls = wall.get_lines()
        lc = Line3DCollection(ls, linewidths=wall_linewidth, colors='k', alpha=wall_alpha)
        ax.add_collection(lc)

    # mx = np.amax(np.abs(P))
    mx = np.amax(np.abs(wall.get_lrtp()))
    XYZlim = [-mx, mx]
    ax.set_xlim3d(XYZlim)
    ax.set_ylim3d(XYZlim)
    ax.set_zlim3d(XYZlim)
    ax.set_box_aspect((1, 1, 1))
                
    plt.grid()
    plt.savefig(save_filename, dpi=300, bbox_inches='tight')
    if show_plot:
        plt.show()



class Experiment(object):

    """Experiment setup"""

    def __init__(self, particles, wall, contact):
        """Collection of all

        :particles: TODO
        :wall: TODO
        :contact: TODO

        """
        self.particles = particles
        self.wall = wall
        self.contact = contact

        
        ## compute the boundary nodes, for each particle count
        # print('Computing boundary nodes: ')
        # for sh in range(len(self.particles)):
            # count = len(self.particles[sh])
            # for i in range(count):
                # particle = self.particles[sh][i]
                # temp = []
                # for p in range(len(particle.edge_nodes)):
                    # e_node = particle.edge_nodes[p]
                    # pos_node = particle.pos[e_node]
                    # for q in range(len(particle.pos)):
                        # pos_other = particle.pos[q]
                        # d = np.sqrt(np.sum(np.square(pos_node - pos_other))) #norm
                        # # if (d <= self.contact.contact_radius):
                        # if (d <= self.contact.contact_radius*1.1):
                            # temp.append(q)
                # particle.nonlocal_bdry_nodes = list(set(temp))
                ## particle.nonlocal_bdry_nodes = list(set(range(len(particle.pos))))

        
                # print('edge nodes: ', len(particle.edge_nodes))
                # print('bdry nodes: ', len(particle.nonlocal_bdry_nodes))
                # print(i, [len(particle.edge_nodes), len(particle.nonlocal_bdry_nodes)], 
                # print(i, end = ' ', flush=True)

        # print('')



    def save(self, filename):
        """Save experiment setup to file

        :filename: filename with path
        :returns: TODO

        """
        print('Saving to', filename,  'Univ particle: ')
        with h5py.File(filename, "w") as f:
            # universal index starts at 1 to be compatible with matlab
            # j = 1
            # universal index starts at 0 to be compatible with C++
            j = 0

            for sh in range(len(self.particles)):
                count = len(self.particles[sh])
                for i in range(count):
                    particle = self.particles[sh][i]
                    p_ind = ('P_%05d' % j)
                    f.create_dataset(p_ind + '/Pos', data=particle.pos)
                    f.create_dataset(p_ind + '/Vol', data=particle.vol)
                    f.create_dataset(p_ind + '/Connectivity', data=particle.NArr)
                    # needs to be a column matrix
                    f.create_dataset(p_ind + '/bdry_nodes', data=np.array([particle.nonlocal_bdry_nodes]).transpose())
                    f.create_dataset(p_ind + '/bdry_edges', data=np.array(particle.bdry_edges))

                    f.create_dataset(p_ind + '/clamped_nodes', data=np.array(particle.clamped_nodes))

                    # Initial conditions
                    f.create_dataset(p_ind + '/disp', data=particle.disp)
                    f.create_dataset(p_ind + '/vel', data=particle.vel)
                    f.create_dataset(p_ind + '/acc', data=particle.acc)
                    f.create_dataset(p_ind + '/extforce', data=particle.extforce)

                    # Material properties
                    f.create_dataset(p_ind + '/delta', data=[[particle.material.delta]])
                    f.create_dataset(p_ind + '/rho', data= [[particle.material.rho]])
                    f.create_dataset(p_ind + '/cnot', data= [[particle.material.cnot]])
                    f.create_dataset(p_ind + '/snot', data= [[particle.material.snot]])

                    # torque info
                    f.create_dataset(p_ind + '/torque_axis', data= [[particle.torque_axis]])
                    f.create_dataset(p_ind + '/torque_val', data= [[particle.torque_val]])

                    # extra properties
                    f.create_dataset(p_ind + '/movable', data= [[particle.movable]])
                    f.create_dataset(p_ind + '/breakable', data= [[particle.breakable]])
                    f.create_dataset(p_ind + '/stoppable', data= [[particle.stoppable]])

                    print(j, end = ' ', flush=True)
                    j = j+1

            print('\n')

            # making all scalars a rank-2 data (matrix) to make it compatible with matlab
            # f.create_dataset('total_particles_univ', data=[[j-1]])
            f.create_dataset('total_particles_univ', data=[[j]])

            # contact info
            f.create_dataset('pairwise/contact_radius', data=[[self.contact.contact_radius]])
            f.create_dataset('pairwise/normal_stiffness', data=[[self.contact.normal_stiffness]])
            f.create_dataset('pairwise/damping_ratio', data=[[self.contact.damping_ratio]])
            f.create_dataset('pairwise/friction_coefficient', data=[[self.contact.friction_coefficient]])

            # wall info
            f.create_dataset('wall/allow_wall', data = [[self.wall.allow]])
            f.create_dataset('wall/geom_wall_info', data = np.array([self.wall.get_lrtp()]).transpose() )

