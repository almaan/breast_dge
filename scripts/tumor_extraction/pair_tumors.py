#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 14 12:36:08 2019

@author: alma
"""

from funcs.utils import *

import argparse as arp
import os.path as osp

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import cdist, pdist


def get_centroids(x,y,idx,ignore_neg = True):
    """
    Returns 2d centroids for each subset of within data. Labels
    are provided by idx. x and y are respective coordinates.
    If ignore_neg is used, centroids for subsets with negative 
    (non-tumor) labels will not be computed.
    
    Centroids are given in order of increasing label id. 
    Centroids are not necessarily set members but rather
    the unweighted mean coordinates of all members.
    
    
    """
    centroids = []
    # get unqiue id's within sample. Ascending order.
    labels = np.unique(idx)
    
    # if negative labels should not be considered
    if ignore_neg:
        labels = np.array([x for x in labels if x >= 0])

    # iterate over all labels within set
    for label in labels:
        # get centroid  
        x_m = np.mean(x[idx == label])
        y_m = np.mean(y[idx == label])
        
        centroids.append([x_m,y_m])
    
    # convert to numpy array
    centroids = np.array(centroids)    
    
    return centroids 


    
def pair_clusters_centroid(samples,
                           colname):
    """
    Returns modified dataframes of provided list containing
    patient replicates within which tumors are to be paired.
    
    """
    
    #used for debug purpose. Remove later
    debug_plot = False
    
    # deep copy coordinates for all replicates
    pcs = [samples[ii].loc[:,['xcoord','ycoord']].copy().values for ii in range(len(samples))]
    
    # center all coordinates around origo
    # TODO: this could be stored and added back on before return
    pcs = [pcs[ii] - pcs[ii].mean(axis=0) for ii in range(len(pcs))]
    
    # get centroids for all replicates and all the tumors within each replicate
    # centroids will be list of 2d arrays
    centroids = [get_centroids(samples[ii]['xcoord'].values,samples[ii]['ycoord'].values,
                               samples[ii][colname]) for ii in range(len(samples))]

    # get index of sample with lowest cardinality of clsuters
    # this will be lower bound (reference) for all pairing 
    reference = np.argmin([cent.shape[0] for cent in centroids])
    
    # get indices of all replicates to _be_ aligned to reference
    # refered to as targets
    targets = np.delete(np.arange(len(centroids)), reference) 
    
    #for debugging purposes. Remove later
    if debug_plot:
        for cs in pcs:
            plt.scatter(cs[:,0], cs[:,1])
            plt.show()  
    
    # get all cluster labels of reference replicate
    ref_labels = np.unique(samples[reference][colname])
    # remove negatively labeled clusters (non-tumor)
    ref_labels = np.delete(ref_labels,np.where(ref_labels <  0))
    
    # iterate over all replicates to be aligned to reference
    for samp in targets:
        # generate distance matrix between reference and target
        # cdist uses l2 distance between centroids         
        dmat = cdist(centroids[reference],
                     centroids[samp])
        
        pairs = [] #   to store which label mapping
        
        # iterate over all clusters within reference
        # is per construction less than or equal to target 
        for d in range(centroids[reference].shape[0]):
            # pairs tumors with nearest centroids wihtin target and reference 
            row,col = np.unravel_index(np.argmin(dmat),(dmat.shape))
            pairs.append([row,col]) #  save pairing
            
            # do not "double" assign tumors to same label
            dmat[row,:] = np.inf 
            dmat[:,col] = np.inf
        
        # get labels of target
        target_labels = np.unique(samples[samp][colname])
        # remove negative labels (non-tumor)
        target_labels = np.delete(target_labels,np.where(target_labels < 0))
        
        # create temporary vector for new labels
        # TODO: currently unpaired clusters are set to -1; find better solution
        tmp = np.repeat(-1,samples[samp].values.shape[0])
        
        
        # iterate over all mapped pairs 
        for p in pairs:
            # set new labels for target samples, corresponding to 
            # label of paired tumor in reference
            tmp[samples[samp][colname].values == target_labels[p[1]]] = target_labels[p[0]]
    
        # in-place modification of sample labels
        samples[samp].loc[:,colname] = tmp
        
    
    return samples
                        


def pair_clusters_members(dfs,colname):
   """
   
   Pairing based on number of members within each cluster.
   Failed miserably
   
   """
   for df in dfs:
        sample_labels, counts = np.unique(df[colname][df[colname] != -1], return_counts=True)
        sample_labels = sample_labels[np.argsort(counts)]
        tmp = np.repeat(-1, df.shape[0])
        
        for (k,newlabel) in enumerate(sample_labels):
            tmp[df[colname] == newlabel] = k
            
        df[colname] = tmp
    
   return dfs

    
def main(input_name, odir, colname):
    
    data = [] #  will hold data frames of replicates
    initial_clusters  = [] # will hold initial 
    
    # read files within input_directory
    for file in input_name:
            data.append(load_file(file)) #  load data-frame from filename
            idx = np.unique(data[-1][colname]) #  find numbers of clusters within sample
            idx = np.delete(idx,np.where(idx == -1)) #  remove negative labels (non-tumor)
            nclust = idx.shape[0] #  number of clusters within set
            
            # if no tumor spots do not use sample
            if nclust < 1: 
                data = data[0:-1] #  removes last appended data frame
            else:
                # register number of clusters within dataset
                initial_clusters.append(idx.shape[0])
            
    # convert to numpy array
    initial_clusters = np.array(initial_clusters)
    
    # process replicates jointly
    if initial_clusters.shape[0] > 1:
        # if at least two replicates contain non-negatively labeled spots
        dfs = pair_clusters_centroid(data,colname) #  perform pairing between replicates
    else:
        # if pairing is not required. Either only negatively labeled spots or less than
        # two replicates to pair 
        dfs = data #  keeps data as it is
        
    # if any samples contain non-negative labeled spots
    if len(data) > 0:
        suffix = '_paired'
        # TODO: make palette based on number of clusters. Not hard-coded            
        # colormap
        mymap = ['palevioletred','cadetblue','indigo','tan','green','yellow']
        
        # iterate over all loaded replicates
        for num in range(len(dfs)):
            
            # filename for tsv file of results
            df_oname = osp.join(odir,''.join(['.'.join(osp.basename(input_name[num]).split('.')[0:-1]),suffix,'.tsv']))
            # filename of png file of visual result
            plt_oname = osp.join(odir,''.join(['.'.join(osp.basename(input_name[num]).split('.')[0:-1]),suffix,'.png']))
            
            # create separate image for each replicate
            fig, ax = plt.subplots(1,1)
            cidx = 0
            for ii in np.unique(dfs[num][colname]):
                
                # scatter plot for visualization
                ax.scatter(dfs[num].loc[:,'xcoord'].values[dfs[num][colname] == ii],
                           dfs[num].loc[:,'ycoord'].values[dfs[num][colname] == ii],
                           s = 100,
                           c = (mymap[cidx] if ii >= 0 else 'black'), #  non-tumors as black
                           alpha = (1.0 if ii >= 0 else 0.3), # non-tumors have lower alpha
                           edgecolor = 'black',
                           )
                
                cidx = (cidx + 1 if ii > -1 else cidx) # for color consistency
            
            # save results
            dfs[num].to_csv(df_oname, sep = '\t', header = True, index = True)
            # save visualization of results    
            fig.savefig(plt_oname)
            
if __name__ == '__main__':

# Parser ----------------------------
    
    #TODO: add help information and write this up nicer
    
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

    arguments = dict(input_name = args.input,
                     odir = args.odir,
                     colname = args.feature)    
    
    main(**arguments)

