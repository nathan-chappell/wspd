# multithreaded...

# WSPD with dimensionality reduction and evaluation of separation in orginal
# space

import numpy as np
from matplotlib import pyplot as plt
from sklearn.decomposition import PCA
from itertools import product
from concurrent.futures import ThreadPoolExecutor
from concurrent.futures import ProcessPoolExecutor
import os
import json

#{{{ Globals

# here are the parameters of the experiment:
# filenames target_seps n_components
# set dummy to True to only do a test run (one file, fewer parameters)

data_dir = '../data4Domagoj/'

#filenames = [
  #'biase-prepare-log_count.csv',
  #'darmanisfilledSOUP-prepare-log_count.csv',
  #'darmanis-prepare-log_count.csv',
  #'deng-reads-prepare-log_count.csv',
  #'deng-rpkms-prepare-log_count.csv',
  #'fan-prepare-log_count.csv',
  #'goolam-prepare-log_count.csv']

filenames = [
  'yan-prepare-log_count.csv',
  'treutlein-prepare-log_count.csv'
]

# maybe too optimistic...
target_seps = [.25, .5, 1, 1.5]
pca_size = [5, 10, 20]
      
# dummy runs are to test the overall setup
dummy = False
#dummy = True
if dummy:
  filenames = [
    '../data4Domagoj/biase-prepare-log_count.csv',
    ]
  target_seps = [1, 1.5]
  pca_size = [5, 10]

#}}}

#{{{ Routines
def print_status(i,n):
  #print('\r%8d : %d'%(i+1,n),end='',flush=True)
  pass

# runs concurrently
def get_and_add_box(dumbell_indices,boxes,data, i):
  l,r = dumbell_indices[i]
  # loop loop loop loop
  lbox = np.array([[min(x),max(x)] for x in data[l].T])
  rbox = np.array([[min(x),max(x)] for x in data[r].T])
  boxes[i] = [lbox,rbox]
  #print_status(len(boxes), len(dumbell_indices))

# runs concurrently
def get_and_add_hi_d_dumbell(hi_d_dumbells_dict, boxes, i):
  l,r = boxes[i]
  # loop loop
  lc = (l[:,0] + l[:,1])/2
  lr = np.linalg.norm(l[:,0] - lc) 
  # loop loop
  rc = (r[:,0] + r[:,1])/2
  rr = np.linalg.norm(r[:,0] - rc)
  if lr == 0 and rr == 0:
    sep = -1
  else: 
    sep = get_sep(lc,lr,rc,rr)
  hi_d_dumbells_dict[i] = [sep,lr,rr]
  #print_status(len(hi_d_dumbells_dict),len(boxes))

# get separation given centers and radii
def get_sep(lc,lr,rc,rr):
  r = max(lr,rr)
  if r == 0: return 0
  # loop
  #return max(0, np.linalg.norm(lc - rc)/r - 2)
  return np.linalg.norm(lc - rc)/r - 2

def info_to_dict(info):
  return { 
     's' : info[0],
    'lr' : info[1],
    'rr' : info[2]
  }

def record_result(results, filename, shape, target_sep, n_components,
                  dumbell_indices, low_d_info, hi_d_info):
  results.append({
    'filename' : filename,
    'shape' : shape,
    'target_sep' : target_sep,
    'n_components' : n_components,
    'dumbell_indices' : dumbell_indices,
    'low_d_info' : low_d_info,
    'hi_d_info' : hi_d_info
    })

def print_index_count(i,a):
  #print(a.index(i), '/', len(a))
  pass

#}}}

def run_experiment(filename):

### Read Data
#{{{
  data = np.loadtxt(data_dir + filename, delimiter=',')
  #print('shape:', data.shape, 'max pairs:', data.shape[0]*(data.shape[0]-1)/2)
#}}}

  # list of dictionaries holding results of the experiments
  results = []
  for target_sep, n_components in product(target_seps, pca_size):
    #print('----------------------------------------\n')
    print('processing: ', filename, 'sep:', target_sep, ', n_c:', n_components)
    #print('\n----------------------------------------')

### Preprocessing
#{{{
# reduce dimension
    #pca = PCA(n_components=n_components, svd_solver='arpack')
    pca = PCA(n_components=n_components)
    transformed_data = pca.fit_transform(data)

# output transformed data
    shape = transformed_data.shape
    header = str(shape[0]) + ' ' + str(shape[1])
    pca_out = 'results/' + filename + '.pca_out.np'
    np.savetxt(pca_out, transformed_data, header=header)
#}}}

### WSPD
#{{{
# run wspd on transformed data
    os.system('./wsp ' + pca_out + ' ' + str(target_sep)
                       + '1>/dev/null 2>/dev/null')

    print(filename, 'wspd done')
# read wspd result
    dumbell_indices = []
    low_d_info = []

# dumbell_indices: s lc lr rc rr | l1 l2 ... | r1 r2 ...
    with open(pca_out +  '.wsp_out.txt') as f:
      for line in f:
        i,l,r = line.split('|')
        low_d_info.append(info_to_dict([float(x) for x in i.split()]))
        dumbell_indices.append([[int(x) for x in l.split()],[int(y) for y in r.split()]])

    print(filename, 'wspd done,', len(dumbell_indices), 'pairs')
#}}}

### Analysis
#{{{ 
    #print('getting bounding boxes')

    boxes = {}
    [get_and_add_box(dumbell_indices,boxes,data,i) for i in range(len(dumbell_indices))]
    print(filename, 'bounding boxes done')

    #print('\nboxes:', len(boxes))
    #print('\ngetting hi-d separation')

    hi_d_dumbells_dict = {}
    [get_and_add_hi_d_dumbell(hi_d_dumbells_dict, boxes, i) for i in range(len(boxes))]
    #print('...dun')
    print(filename, 'hi_d_dumbells done')
      
    # move the dict to a list, then each item of the list to a dict
    hi_d_info = []
    for i in range(len(hi_d_dumbells_dict)): 
      hi_d_info.append(info_to_dict(hi_d_dumbells_dict[i]))

    print('----------------------------------------')
    print('\nrecording result for: ', filename, target_sep, n_components)
    print('----------------------------------------')

    record_result(results, filename, data.shape, target_sep, n_components,
                  dumbell_indices, low_d_info, hi_d_info)

#}}}

### Output the results
#{{{

  print('****************************************')
  print('\noutputting results for: ', filename)
  print('****************************************')
  with open('results/' + filename + '.results.json','w') as f: json.dump(results,f)

#}}}

### Run Experiments

with ProcessPoolExecutor() as pool:
  pool.map(run_experiment, filenames)

#with ThreadPoolExecutor() as pool:
  #pool.map(run_experiment, filenames)

#for filename in filenames:
  #run_experiment(filename)

### End of Experiement

