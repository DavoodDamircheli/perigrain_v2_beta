import sys, os
sys.path.append(os.getcwd())
import numpy as np
import h5py
import re
# from multiprocessing import Pool
from pathos.multiprocessing import ProcessingPool as Pool
import time

from scipy.signal import savgol_filter
# yhat = savgol_filter(y, 51, 3) # window size 51, polynomial order 3
from scipy.stats import linregress

import matplotlib.pyplot as plt
import matplotlib
# matplotlib.use('Agg')
from matplotlib.collections import LineCollection

import argparse
# Instantiate the parser
parser = argparse.ArgumentParser(description='Optional app description')

parser.add_argument('--timestep_file', type=str, help='output setup file', default='data/hdf5/tc_00000.h5')
parser.add_argument('--save_to', type=str, help='output setup file', default='output/csv/arrangement.csv')
args = parser.parse_args()

filename = args.timestep_file

f = h5py.File(filename, "r")

loc = []
rad = []
for name in f:
    if re.match(r'P_[0-9]+', name):
        # pid = int(name[2:])
        currpos = np.array(f[name+'/CurrPos'])
        centroid = np.mean(currpos, axis=0)
        # print(centroid)
        loc.append(centroid)
        # r_max = np.sqrt(np.max(np.sum((currpos - centroid)**2, axis=1), axis=0))
        r = np.sqrt(np.max(np.sum((currpos - centroid)**2, axis=1)))
        # print(s)
        rad.append(r)

# loc = np.array(loc)
# rad = np.array(rad)

# print(loc)
# print(rad)
arr = np.c_[ loc, rad ]
print('saving to:', args.save_to)
np.savetxt(args.save_to, arr, delimiter=",")


