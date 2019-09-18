# analysis.py

import json
import numpy as np
import matplotlib.pyplot as plt

results_dir = 'results/'

results_filenames = [
  'biase-prepare-log_count.csv.results.json',
  'darmanisfilledSOUP-prepare-log_count.csv.results.json',
  'fan-prepare-log_count.csv.results.json'
]

results_jsons = []

for filename in results_filenames:
  with open(results_dir + filename) as f: results_jsons.append(json.load(f))

# scatter the old separation vs new separation

def print_basic_info(result):
  print(result['filename'])
  print('n,d:',result['shape'],'s:',result['target_sep'],
        'n_components:', result['n_components'],
        'dumbells:',len(result['dumbell_indices']))

def get_np_for_dumbells(dumbell_info):
  return np.array([[db['s'],db['lr'],db['rr']] for db in dumbell_info])

target_sep_grid_map = {
  .25:0,
  .5:1,
  1:2,
  1.5:3
}

n_components_grid_map = {
  5:1,
  10:2,
  20:3
}

columns = len(target_sep_grid_map)
rows = len(n_components_grid_map)
def get_subplt_index(target_sep, n_components):
  return target_sep_grid_map[target_sep]*rows + n_components_grid_map[n_components]

for results in results_jsons:
  for result in results:
    print_basic_info(result)
    target_sep = result['target_sep']
    n_components = int(result['n_components'])
    low_d = get_np_for_dumbells(result['low_d_info'])
    hi_d = get_np_for_dumbells(result['hi_d_info'])
    num_dumbells = len(low_d)
    #plt.subplot(rows,columns,get_subplt_index(target_sep,n_components))
    plt.scatter(np.arange(num_dumbells), 
        low_d[:,0],c=[(0,0,1,.8)],label='projected')
    plt.scatter(np.arange(num_dumbells), 
        hi_d[:,0],c=[(1,0,0,.8)], label='original space')
    plt.plot(np.arange(num_dumbells),np.ones(num_dumbells)*target_sep,ls='--',c='k',label='target sep')
    plt.plot(np.arange(num_dumbells),np.ones(num_dumbells)*0,c='k',label='zero sep')
    plt.yscale('symlog',basey=2)
    plt.ylabel('actual separation')
    plt.xlabel('dumbells')
    plt.title('shape:' + str(result['shape']) + ', dumbells: ' +
                  str(len(result['dumbell_indices'])) + '\n' +
                  'target_sep: ' + str(target_sep) + ', n_components: ' + 
                  str(n_components))
    #plt.savefig('arpack_sep_' + str(result['target_sep']) + '.png')
    plt.legend()
    plt.savefig('results/' + result['filename'] + 'sep_' + 
                 str(target_sep) + '_' +
                 str(n_components) + 
                 '.png')
    plt.clf()
    #plt.show()

#plt.savefig('combined.png')


