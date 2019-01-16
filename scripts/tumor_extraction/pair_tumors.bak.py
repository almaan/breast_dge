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


def get_centroids(df,colname, ignore_neg = True):

    centroids = []
    labels = np.unique(df[colname].values)
    if ignore_neg:
        labels = np.array([x for x in labels if x >= 0])

    for label in labels:
        x_m = np.mean(df['xcoord'].values[df[colname] == label])
        y_m = np.mean(df['ycoord'].values[df[colname] == label])
        
        centroids.append([x_m,y_m])
    
    centroids = np.array(centroids)    
    
    return centroids 

def pair_clusters_centroid(dfs,colname):

    
    def pair_clusters(sample1,
                      sample2,
                      colname,
                      kabsch = False):
        
        
        pcs1 = sample1.loc[:,['xcoord','ycoord']].copy().values
        pcs2 = sample2.loc[:,['xcoord','ycoord']].copy().values
    
        if kabsch:
            min_members = np.min([sample1.shape[0], sample2.shape[0]])
            
            idx_1 = np.random.choice(np.arange(sample1.shape[0]),min_members, replace = False)
            idx_2 = np.random.choice(np.arange(sample2.shape[0]),min_members, replace = False)
            
            cs1 = pcs1[idx_1,:]
            cs2 = pcs2[idx_2,:]
            
            cs1 = cs1 - np.mean(cs1,axis = 0)
            cs2 = cs2 - np.mean(cs2,axis = 0)
            
            
            H = np.dot(np.transpose(cs2),cs1)
            V,S,W = np.linalg.svd(H)
            d = (np.linalg.det(V) * np.linalg.det(W)) < 0.0
            
            if d:
                V[:, -1] = -V[:, -1]

               
            R = np.dot(V,W)
        
            pcs2 = np.dot(pcs2 - np.mean(pcs2),R)
            pcs1 = pcs1 - np.mean(pcs1)
        
        
        sample1.loc[:,['xcoord','ycoord']] = pcs1 - np.mean(pcs1, axis = 0) 
        sample2.loc[:,['xcoord','ycoord']] = pcs2 - np.mean(pcs2, axis = 0 )
        
        dmat = cdist(get_centroids(sample1,colname),get_centroids(sample2,colname))
        
        pairs = []    
        
        for d in range(dmat.shape[0]):
            row,col = np.unravel_index(np.argmin(dmat),(dmat.shape[0],dmat.shape[0]))
            pairs.append([row,col])
            
            dmat[row,:] = np.inf
            dmat[:,col] = np.inf
        
        return pairs
    
                        
    labels_0 = np.unique(dfs[0][colname])
    labels_0 = (labels_0[1::] if -1 in labels_0 else labels_0)
    
    for k in range(1,len(dfs)):
        
        pairs = pair_clusters(dfs[0],dfs[k],
                              colname = colname,
                              kabsch= False)
        
        labels_k = np.unique(dfs[k][colname])
        if -1 in labels_k:
            labels_k = np.delete(labels_k,np.where(labels_k == -1))
        tmp = np.repeat(-1,dfs[k].values.shape[0])
        
        for p in pairs:
            tmp[dfs[k][colname].values == labels_k[p[1]]] = labels_0[p[0]]
        
        dfs[k].loc[:,colname] = tmp
    
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
    
    
    if (np.sum(np.abs(np.diff(initial_clusters))) != 0 if  initial_clusters.shape[0] > 0 else False):
        
        min_clusters = np.min(initial_clusters)
        for df in data:
            centroids.append(get_centroids(df,colname))
        
        for (k,df) in enumerate(data):
            nclust = centroids[k].shape[0]
            if nclust > min_clusters:
                Z = linkage(pdist(centroids[k]), method = 'average', metric = 'euclidean')
                idx = nclust
                for ii in range(nclust - min_clusters):
                    df.loc[df[colname] == int(Z[ii,0]),colname] = idx     
                    df.loc[df[colname] == int(Z[ii,1]),colname] = idx
                    idx = idx + 1
                    
            new_idx = np.repeat(-1,df.shape[0])

            for (k,ii) in enumerate(np.unique(df[colname])):
                if ii > -1:
                    new_idx[df[colname] == ii] = k

            df.loc[:,colname] = new_idx
    
    if (np.min(initial_clusters) > 2 if initial_clusters.shape[0] > 0 else False):        
        dfs = pair_clusters_centroid(data,colname)
    else:
        dfs = data
        
    if len(data) > 0:
        suffix = '_paired'            
        mymap = ['tan','green','yellow','palevioletred','cadetblue','indigo']
            
        for num in range(len(dfs)):
            df_oname = osp.join(odir,''.join(['.'.join(osp.basename(input[num]).split('.')[0:-1]),suffix,'.tsv']))
            plt_oname = osp.join(odir,''.join(['.'.join(osp.basename(input[num]).split('.')[0:-1]),suffix,'.png']))
            
            fig, ax = plt.subplots(1,1)
            
            for (cidx,ii) in enumerate(np.unique(dfs[num][colname])):
                
                ax.scatter(dfs[num].loc[:,'xcoord'].values[dfs[num][colname] == ii],
                           dfs[num].loc[:,'ycoord'].values[dfs[num][colname] == ii],
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

