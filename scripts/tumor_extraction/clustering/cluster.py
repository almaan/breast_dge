#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 26 17:07:28 2018

@author: alma
"""
import logging
logger = logging.getLogger('feature_extraction_log')


import matplotlib.pyplot as plt
from funcs.utils import *
import numpy as np





def featureScan(data,
              minSpt,
              maxDist,
              minFeature,
              ):
    
    """
    Clustering algorithm. Returns a list of cluster indices where -1 represents
    noise or no feature. 
    
    Arguments:
        minSpt      - minimum number of neighbours for a spot to be considered an inner point
        maxDist     - maximum distance between two spots to be considered neighbours
        minFeature  - minimum number of spots within a cluster to be considered a proper cluster.
        
    Returns:
        clusters    - numpy array of cluster indices
    
    """
    
    def propagate(idx,
                  data,
                  minSpt,
                  maxDist,
                  ori_clusters,
                  ):
        """
        Propagation of inner points as to expand a cluster.
        
        """
        
        clusters = np.copy(ori_clusters)
        in_range = distance(data[idx,:],data) < maxDist
        
        idx_list = np.where(in_range)
        idx_list = np.delete(idx_list, np.argmax(idx_list == idx))
        
        Neps = in_range.sum()
        
        if Neps >= minSpt:

            if clusters[idx] == -1:
                cidx = np.max(clusters) + 1
                clusters[idx] = cidx
            else:
                cidx = clusters[idx]
            
            clusters[idx_list] = cidx
            
            while len(idx_list) > 0:
                in_range = distance(data[idx_list[0],:],data) < maxDist
                Neps = in_range.sum()
                               
                if Neps >= minSpt:
                
                    in_range_and_new = (in_range)*(clusters != cidx)
                    idx_list = np.append(idx_list, np.where(in_range_and_new))
                    
                    clusters[np.where(in_range_and_new)] = cidx
                    
                idx_list = np.delete(idx_list, 0)
                
        return clusters
            
    clusters = np.ones(data.shape[0])*(-1)

    for ii in range(clusters.shape[0]):
        if clusters[ii] == -1:
            clusters = propagate(ii, data, minSpt, maxDist, clusters)
    
    if minFeature > 0:
        for cidx in np.unique(clusters):
            if np.sum(clusters == cidx) < minFeature:
                clusters[clusters == cidx] = -1
    
    return clusters
        
def extractFeatures(df,
                  minSpt,
                  maxDist,
                  minFeature,
                  feature,
                  select_for,
                  ):
    
    """
    Extract features from a dataframe generated from annotated data. Only passes those spots
    that are annotated with the feature of interest (to select for) to the featureScan function.
    
     Arguments:
        minSpt      - minimum number of neighbours for a spot to be considered an inner point
        maxDist     - maximum distance between two spots to be considered neighbours
        minFeature  - minimum number of spots within a cluster to be considered a proper cluster.
        feature     - name of column where annotations are found
        select_for  - feature to select for, feature of interest
        
    Returns:
        all_cluster_labels - array of clusters labels for all spots both with feature and without
        
    """
    crd = get_coordinates(df)
    
    all_cluster_labels = np.ones(crd.shape[0])*(-1)
    tumor_idx = np.where((df[feature] == select_for))
    tumor_cluster_labels = featureScan(crd[tumor_idx], minSpt=minSpt, maxDist=maxDist, minFeature=minFeature)
    
    all_cluster_labels[tumor_idx] = tumor_cluster_labels
    
    return all_cluster_labels