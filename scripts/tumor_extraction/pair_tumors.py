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
from scipy.spatial.distance import cdist


def get_centroids(df,colname):

    centroids = []
    labels = np.unique(df[colname].values)
    labels = np.array([x for x in labels if x >= 0])

    for label in labels:
        x_m = np.mean(df['xcoord'].values[df[colname] == label])
        y_m = np.mean(df['ycoord'].values[df[colname] == label])
        centroids.append([x_m,y_m])
    
    centroids = np.array(centroids)
    centroids = centroids - np.mean(centroids,axis=0)
    
    return centroids 

def pair_clusters_centroid(dfs,colname):

    
    def pair_clusters(sample1,
                      sample2,
                      ncentroids,
                      center = False,
                      kabsch = False):
        
        if center:
            cs1 = sample1 - np.mean(sample1,axis = 0)
            cs2 = sample2 - np.mean(sample2,axis = 0)
        else:
            cs1 = np.copy(sample1)
            cs2 = np.copy(sample2)
            
        if kabsch:
            H = np.dot(cs1.T,cs2)
            u,s,vh = np.linalg.svd(H)
            d = np.linalg.det(np.dot(vh.T,u.T))
            D = np.zeros((2,2))
            D[0,0] = 1
            D[1,1] = d
            R = np.dot(np.dot(vh.T,D),u.T)
            
        
            cs2 = np.dot(cs2,R)
        
        dmat = cdist(cs1,cs2)
        pairs = []    
        
        for d in range(ncentroids):
            row,col = np.unravel_index(np.argmin(dmat),(ncentroids,ncentroids))
            pairs.append([row,col])
            
            dmat[row,:] = np.inf
            dmat[:,col] = np.inf
        
        return pairs
    
    centroids =  [get_centroids(x,colname) for x in dfs]
        
    labels_0 = np.unique(dfs[0][colname].values)
    labels_0 = (labels_0[1::] if -1 in labels_0 else labels_0)
    
    for k in range(1,len(dfs)):
        
        pairs = pair_clusters(centroids[0],centroids[k],
                              ncentroids=centroids[0].shape[0],
                              center=False,
                              kabsch= False)
        
        labels_k = np.unique(dfs[k][colname].values)
        labels_k = (labels_k[1::] if -1 in labels_k else labels_k)
        tmp = np.repeat(-1,dfs[k].values.shape[0])
        for p in pairs:
            tmp[dfs[k][colname].values == labels_k[p[1]]] = labels_0[p[0]]
        
        dfs[k][colname] = tmp
    
    return dfs

def pair_clusters_members(dfs,colname):
    
    for df in dfs:
        sample_labels, counts = np.unique(df[colname][df[colname] != -1], return_counts=True)
        sample_labels = sample_labels[np.argsort(counts)]
        tmp = np.repeat(-1, df.shape[0])
        for (k,newlabel) in enumerate(sample_labels):
            tmp[df[colname] == newlabel] = k
            
        df[colname] = tmp
    
    return dfs

#%%
    
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
    centroids = []
    
    
    if np.sum(np.abs(np.diff(initial_clusters))) != 0 and initial_clusters.shape[0] > 0:
        
        min_clusters = np.min(initial_clusters)
        for df in data:
            centroids.append(get_centroids(df,colname))
        
        for (k,df) in enumerate(data):
            nclust = len(centroids[k])
            if nclust > min_clusters:
                Z = linkage(centroids[k], method = 'single', metric = 'euclidean')
                idx = nclust
                for ii in range(nclust - min_clusters):
                    df[colname][df[colname] == int(Z[ii,0])] = idx     
                    df[colname][df[colname] == int(Z[ii,1])] = idx
                    idx = idx + 1
                
    
    if np.min(initial_clusters) > 2 and initial_clusters.shape[0] > 0:        
        dfs = pair_clusters_centroid(data,colname)
    else:
        dfs = data
    suffix = '_paired'
    
    colormap = plt.get_cmap(plt.cm.Dark2_r)
    unique_labels = np.unique(dfs[0][colname])
    nlabels = unique_labels.shape[0]
    mymap = (colormap(np.linspace(0.1,1,nlabels+1)) if nlabels > 1 else 
                 (['red'] if unique_labels[0] >= 0 else ['green']))
        
    for num in range(len(dfs)):
        df_oname = osp.join(odir,''.join(['.'.join(osp.basename(input[num]).split('.')[0:-1]),suffix,'.tsv']))
        plt_oname = osp.join(odir,''.join(['.'.join(osp.basename(input[num]).split('.')[0:-1]),suffix,'.png']))
        
        fig, ax = plt.subplots(1,1)
        
        for (cidx,ii) in enumerate(np.unique(dfs[num][colname])):
            
            ax.scatter(dfs[num]['xcoord'].values[dfs[num][colname] == ii],
                       dfs[num]['ycoord'].values[dfs[num][colname] == ii],
                       s = 100,
                       c = (mymap[cidx] if ii >= 0 else 'black'),
                       alpha = (1.0 if ii >= 0 else 0.3),
                       edgecolor = 'black',
                       )
        
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

