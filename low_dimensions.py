# WSPD with dimensionality reduction and evaluation of separation in orginal
# space

import numpy as np
from matplotlib import pyplot as plt
from sklearn.decomposition import PCA
from itertools import product
import os
import json



def print_status(i,n):
  end = ''
  if i == n-1: end='\n'
  print('\r%8d : %d'%(i+1,n),end=end)

# get separation given centers and radii
def get_sep(lc,lr,rc,rr):
  # loop
  d = np.linalg.norm(lc - rc)
  r = max(lr,rr)
  # possibly an error condition...
  if r == 0: 
    #print('r == 0')
    return 0
  d = max(0,d - 2*r)
  #if d == 0:
    #print('d == 0')
  return d/r

# list of dictionaries holding results of the experiments
results = []

def record_result(filename, shape, target_sep, dumbells, actual_seps):
  results.append({
    'filename' : filename,
    'shape' : shape,
    'target_sep' : target_sep,
    'dumbells' : dumbells,
    'actual_seps' : actual_seps})

# read data
#filename = '../data4Domagoj/goolam-prepare-log_count.csv'
#filename = '../data4Domagoj/biase-prepare-log_count.csv'
filenames = [
  '../data4Domagoj/biase-prepare-log_count.csv',
  '../data4Domagoj/darmanisfilledSOUP-prepare-log_count.csv',
  '../data4Domagoj/darmanis-prepare-log_count.csv',
  '../data4Domagoj/deng-reads-prepare-log_count.csv',
  '../data4Domagoj/deng-rpkms-prepare-log_count.csv',
  '../data4Domagoj/fan-prepare-log_count.csv',
  '../data4Domagoj/goolam-prepare-log_count.csv']

target_seps = [.5, 1, 2]
      
dummy = True
if dummy:
  filenames = [
    '../data4Domagoj/biase-prepare-log_count.csv',
    ]
  target_seps = [.5,1]

first = True
for filename,target_sep in product(filenames, target_seps):
  print('processing: ', filename, ', sep:', target_sep)

  # loop
  data = np.loadtxt(filename, delimiter=',')
  print('shape:', data.shape, 'max pairs:', data.shape[0]*(data.shape[0]-1)/2)

### Preprocessing

# reduce dimension
  pca = PCA(n_components=10)
  transformed_data = pca.fit_transform(data)

# output transformed data
  shape = transformed_data.shape
  header = str(shape[0]) + ' ' + str(shape[1])
  pca_out = 'pca_out.np'
  np.savetxt(pca_out, transformed_data, header=header)

### WSPD

# run wspd on transformed data
  os.system('./wsp ' + pca_out + ' ' + str(target_sep) 
                     + '1>/dev/null 2>/dev/null')

# read wspd result
  dumbells = []
  info = []

# dumbells: l1 l2 ... | r1 r2 ...
  with open('dumbells.txt') as f:
    for line in f:
      l,r = line.split('|')
      dumbells.append([[int(x) for x in l.split()],[int(y) for y in r.split()]])

# info: sep l_n r_n l_center r_center l_radius r_radius
  with open('dumbells_s_c_r.txt') as f:
    for line in f:
      info.append([float(x) for x in line.split()])

  print('wspd done. starting analysis')

### Analysis

# get bounding boxes
# TODO too slow...
# consider C++ interface / parallelization
#   really just need min/max...

#f_indices = filter((lambda x: len(x[0]) > 1 or len(x[1]) > 1), dumbells)

# TODO parallelization:
#   get_bounding_box as function
#   iterate through dumbells...

  print('getting bounding boxes')
  boxes = []
#for l,r in f_indices:
  for i in range(len(dumbells)):
    l,r = dumbells[i]
    # loop loop loop loop
    lbox = np.array([[min(x),max(x)] for x in data[l].T])
    rbox = np.array([[min(x),max(x)] for x in data[r].T])
    boxes.append([lbox,rbox])
    print_status(i,len(dumbells))

  print('getting real separations')

  balls = []
  actual_seps = []
# need to account for single points (still good...)
  for i in range(len(boxes)):
    l,r = boxes[i]
    # loop loop
    lc = (l[:,0] + l[:,1])/2
    lr = np.linalg.norm(l[:,0] - lc) 
    # loop loop
    rc = (r[:,0] + r[:,1])/2
    rr = np.linalg.norm(r[:,0] - rc)
    if lr == 0 and rr == 0:
      actual_seps.append(-1)
    else: 
      actual_seps.append(get_sep(lc,lr,rc,rr))
    # as of right now I don't think I want the balls
    # balls.append((lc,lr,rc,rr))
    print_status(i,len(boxes))

  record_result(filename, data.shape, target_sep, dumbells, actual_seps)

#   compare old and new separations:
#     histogram of new separations
#     mean of new separations
#     variance of separation change for nodes
#     -- any patterns in relation to other data

# output information:
# (( a singleton is a dumbell of the form {{a},{b}}
# data filename, n, d, target_sep
#  low_d: count of dumbells, count of singletons, non-singletons, mean size
#  (( in the high_d info, singletons are disregarded ))
#  (( let s := target_sep, s' separation in original space ))
# high_d: 
#  mean(s), (mean size, count) for dumbells 
#                   satisfying s' >= : [0, s/16, s/8, s/4, s/2, s]
# 

with open('low_dimensions.results.json','w') as f: json.dump(results,f)
