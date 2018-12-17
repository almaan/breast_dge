#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  5 08:03:55 2018

@author: Alma Andersson

Get indices for connected subgraphs within a given sample.

"""
import logging

import matplotlib.pyplot as plt

import networkx as nx
import numpy as np
from scipy.spatial import KDTree

logger = logging.getLogger('feature_extraction_log')

def pair_edges(arr, p, thrs):
    """
    generates a list of tuples where each item is the indices of 
    two nodes between which there should be an edge. Nodes are
    connected according if their distance is less than or
    equal to thrs measured in the p-norm.
    
    """
    assert isinstance(arr,np.ndarray)
    #generate kd-tree object from provided array 
    kd = KDTree(arr)
    ngh = []
    #make tuples of neighbouring spots
    for node in range(arr.shape[0]):
        #find all spots witin given distance thrs to a given spot
        #ignore neighbours of lower index than the node as to avoid
        #double edges
        
        ball = filter( lambda x: x > node,
                      kd.query_ball_point(arr[node,:], 
                                          r = thrs, #distance
                                          p = p, #norm
                                          ))
        
        #make list of tuples compatible with graph object edge generation
        pairs = [(node,x) for x in ball]
        ngh = ngh + pairs
        
    return ngh

def connected_graphs_extraction(crd,
                                feature_vector, 
                                select_for,
                                p_norm,
                                maxDist,
                                minFeature):
    
    #make sure that provided data is compatible
    assert crd.shape[0] == feature_vector.shape[0],' '.join(
            ["Number of spots does not",
             "match number of annotated samples"])
    
    idx = np.repeat(-1, crd.shape[0])
    #get position of spots with desired feature
    pos = np.where(feature_vector == select_for)[0]
    
    #make sure that at least one spot for given condition is present
    if sum(pos) > 0:
        #create graph object
        G = nx.Graph()
        #add some number of nodes as spots to graph, spatial information in pos
        G.add_nodes_from(np.arange(crd[pos,:].shape[0]), pos = crd)
        #generate edges between nodes according to their spatial position
        edges = pair_edges(crd[pos,:], p = p_norm, thrs=maxDist)
        G.add_edges_from(edges)
        #get index of connected component for each node
        cc = [list(x) for x in nx.connected_components(G)]
        
        #remove subgraphs that do not have sufficient number of spots
        if minFeature > 0:
            cc = list(filter(lambda x: len(x) >= minFeature,cc))
        #generate an index vector corresponding for the provided data
        for k,ii in enumerate(cc): 
            #only assing spots with non-negative cluster index
            idx[pos[ii]] = k        
    
    return(idx)