# analysis.py

import json
import numpy as np
import matplotlib.pyplot as plt

with open('low_dimensions.results.json') as f: results = json.load(f)

# scatter the old separation vs new separation

def print_basic_info(result):
  print(result['filename'])
  print('n,d:',result['shape'],'s:',result['target_sep'],
        'n_components:', result['n_components'],
        'dumbells:',len(result['dumbell_indices']))

def get_np_for_dumbells(dumbell_info):
  return np.array([[db['s'],db['lr'],db['rr']] for db in dumbell_info])

for result in results:
  print_basic_info(result)
  low_d = get_np_for_dumbells(result['low_d_info'])
  hi_d = get_np_for_dumbells(result['hi_d_info'])
  num_dumbells = len(low_d)
  plt.scatter(np.arange(num_dumbells), low_d[:,0],c='b')
  plt.scatter(np.arange(num_dumbells), hi_d[:,0],c='r')
  plt.yscale('symlog',basey=2)
  #plt.savefig('arpack_sep_' + str(result['target_sep']) + '.png')
  plt.savefig('sep_' + 
               str(result['target_sep']) + '.' +
               str(int(result['n_components'])) + 
               '.png')
  #plt.show()


