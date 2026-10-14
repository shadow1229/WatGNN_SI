import os, glob, time,copy,random
import numpy as np
import pickle
from pathlib import Path

def read_targets(fpath):
    result = []
    f = open(fpath,'r')
    lines = f.readlines()
    for line in lines:
        if line.startswith('#'):
            continue
        else:
            #target = line.strip()[:4]
            target = line.strip()
            tsp = target.split(',')
            if len(tsp) <= 1:
                result.append(target)
            else:
                newname = '%s_%s'%(tsp[0],tsp[1])
                result.append(newname)
    f.close()
    return result


def read_and_sort_water(pdbpath,resname_only=None):
    pdb_f = open(pdbpath,'r')
    
    lines = pdb_f.readlines()
    waters = []
    for line in lines:
        lst = line.strip()
        if lst.startswith('ATOM') or lst.startswith('HETATM'):
            #if len(lst)< 66:
            #    continue
            #if len(lst)> 66:
            #    print(len(lst), pdbpath, lst)
            resname = lst[17:20]
            if resname_only != None and (resname not in resname_only):
                continue

            waters.append(line)

    #waters_sort_tmp = sorted( waters, key = lambda line: -1.0*float(line[60:66]) ) #descending order of occupancy -> will be fixed into B-factor column
    
    waters_sort = ''
    for water_i, water_line in enumerate(waters):
        waters_sort += water_line
    return waters_sort
"""
dir_in  = './gnn_result_wkgb_and_training'
dir_out = './gnn_result_wkgb_patch_sorted'
suffix_in = '_pred.pdb'
suffix_out = '.pdb'
Path(dir_out).mkdir(parents=True,exist_ok=True)

targets = read_targets('wkgb_targets_foldx.txt')
"""
values = ['002','003','004','005','006','007','008','009','010','015','020','025','030','035','040','045','050']
for v in values:
    dir_in  = './predictions_wkgb_%s'%v
    dir_out = './superwater_%s'%v
    suffix_in = '_centroid.pdb'
    suffix_out = '.pdb'
    Path(dir_out).mkdir(parents=True,exist_ok=True)
    targets = read_targets('targets_cmp')

    resname_only = ['HOH']
    for target in targets:
        pdb_in = '%s/%s/%s%s'%(dir_in,target,target,suffix_in)
        pdb_out = '%s/%s%s'%(dir_out,target,suffix_out)
        #if os.access(pdb_out,0):
        #    continue
        if not os.access(pdb_in,0):
            print('not exist - %s',pdb_in)
            continue
        print (target)
        water_pos = read_and_sort_water(pdb_in,resname_only=resname_only)
        pdb_out_f = open(pdb_out,'w')
        pdb_out_f.write(water_pos)
        pdb_out_f.close()


