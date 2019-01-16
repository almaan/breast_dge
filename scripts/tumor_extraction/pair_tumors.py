#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 14 12:36:08 2019

@author: alma
"""

from funcs.utils import *

import argparse as arp
import os.path as osp
import os

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import cdist, pdist


def get_centroids(x,y,idx,ignore_neg = True):

    centroids = []
    labels = np.unique(idx)
    if ignore_neg:
        labels = np.array([x for x in labels if x >= 0])

    for label in labels:
        x_m = np.mean(x[idx == label])
        y_m = np.mean(y[idx == label])
        
        centroids.append([x_m,y_m])
    
    centroids = np.array(centroids)    
    
    return centroids 


def make_weight_matrix(n1,n2):
    weights = np.zeros((n1.shape[0],n2.shape[0]))
    for ii in range(n1.shape[0]):
        for jj in range(n2.shape[0]):
            weights[ii,jj] = np.max([float(n1[ii])/float(n2[jj]),
                                     float(n2[jj])/float(n1[ii])])
    return weights
    
def pair_clusters_centroid(samples,
                           colname):
    
    
    debug_plot = False
    
    pcs = [samples[ii].loc[:,['xcoord','ycoord']].copy().values for ii in range(len(samples))]
    pcs = [pcs[ii] - pcs[ii].mean(axis=0) for ii in range(len(pcs))]
    
    centroids = [get_centroids(samples[ii]['xcoord'].values,samples[ii]['ycoord'].values,
                               samples[ii][colname]) for ii in range(len(samples))]

    reference = np.argmin([cent.shape[0] for cent in centroids])
    to_align = np.delete(np.arange(len(centroids)), reference)
    

    if debug_plot:
        for cs in pcs:
            plt.scatter(cs[:,0], cs[:,1])
            plt.show()  
    
    ref_labels = np.unique(samples[reference][colname])
    ref_labels = np.delete(ref_labels,np.where(ref_labels == -1))
    
    for samp in to_align:        
        dmat = cdist(centroids[reference],
                     centroids[samp])
        
        pairs = []
        for d in range(np.min(dmat.shape)):
            row,col = np.unravel_index(np.argmin(dmat),(dmat.shape))
            pairs.append([row,col])
            dmat[row,:] = np.inf
            dmat[:,col] = np.inf
        
        al_labels = np.unique(samples[samp][colname])
        al_labels = np.delete(al_labels,np.where(al_labels == -1))
        
        tmp = np.repeat(-1,samples[samp].values.shape[0])

        for p in pairs:
            tmp[samples[samp][colname].values == al_labels[p[1]]] = al_labels[p[0]]
    
        samples[samp].loc[:,colname] = tmp
        
    
    return samples
                        

def pair_clusters_members(dfs,colname):
    
    for df in dfs:
        sample_labels, counts = np.unique(df[colname][df[colname] != -1], return_counts=True)
        sample_labels = sample_labels[np.argsort(counts)]
        tmp = np.repeat(-1, df.shape[0])
        
        for (k,newlabel) in enumerate(sample_labels):
            tmp[df[colname] == newlabel] = k
            
        df[colname] = tmp
    
    return dfs

    
def main(input, odir, colname):
    
    data = []
    initial_clusters  = []
    
    for file in input:
            data.append(load_file(file))
            idx = np.unique(data[-1][colname])
            idx = np.delete(idx,np.where(idx == -1))
            nclust = idx.shape[0]
            if nclust < 1:
                data = data[0:-1]
            else:
                initial_clusters.append(idx.shape[0])
            
    initial_clusters = np.array(initial_clusters)
    
    if initial_clusters.shape[0] > 0:        
        dfs = pair_clusters_centroid(data,colname)
    else:
        dfs = data
        
    if len(data) > 0:
        suffix = '_paired'            
        mymap = ['palevioletred','cadetblue','indigo','tan','green','yellow']
            
        for num in range(len(dfs)):
            df_oname = osp.join(odir,''.join(['.'.join(osp.basename(input[num]).split('.')[0:-1]),suffix,'.tsv']))
            plt_oname = osp.join(odir,''.join(['.'.join(osp.basename(input[num]).split('.')[0:-1]),suffix,'.png']))
            
            fig, ax = plt.subplots(1,1)
            cidx = 0
            for ii in np.unique(dfs[num][colname]):
                
                ax.scatter(dfs[num].loc[:,'xcoord'].values[dfs[num][colname] == ii],
                           dfs[num].loc[:,'ycoord'].values[dfs[num][colname] == ii],
                           s = 100,
                           c = (mymap[cidx] if ii >= 0 else 'black'),
                           alpha = (1.0 if ii >= 0 else 0.3),
                           edgecolor = 'black',
                           )
                cidx = (cidx + 1 if ii > -1 else cidx)
                
            fig.savefig(plt_oname)
            
            dfs[num].to_csv(df_oname, sep = '\t', header = True, index = True)
    

#%%

if __name__ == '__main__':
    
    prs = arp.ArgumentParser()
    
    prs.add_argument('-i','--input',
                     required = True,
                     nargs = '*',
                     help = ''.join(['',
                                     ])
                     )
    
    prs.add_argument('-o','--odir',
                     required = True,
                     help = ''.join(['',
                                     ])
                     )

    
    prs.add_argument('-f','--feature',
                     required = False,
                     default = 'tumor_id',
                     help = ''.join(['',
                                     ])
                     )    
    
    args = prs.parse_args()

    arguments = dict(input = args.input,
                     odir = args.odir,
                     colname = args.feature)    
    
    main(**arguments)

