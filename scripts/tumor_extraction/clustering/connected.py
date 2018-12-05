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
    assert isinstance(arr,np.ndarray)
    kd = KDTree(arr)
    ngh = []
    for node in range(arr.shape[0]):
        hood = filter( lambda x: x > node, kd.query_ball_point(arr[node,:], r = thrs, p = p))
        pairs = [(node,x) for x in hood]
        ngh = ngh + pairs
    return ngh

def connected_graphs_extraction(crd, feature_vector, select_for, p_norm, maxDist, minFeature):
    ae1 = "Number of spots does not match number of annotated samples"
    assert crd.shape[0] == feature_vector.shape[0], ae1
    
    pos = np.where(feature_vector == select_for)[0]
    idx = np.repeat(-1, crd.shape[0])
    if sum(pos) > 0:
        edges = pair_edges(crd[pos,:], p = p_norm, thrs=maxDist)
        G = nx.Graph()
        G.add_nodes_from(np.arange(crd[pos,:].shape[0]), pos = crd)
        G.add_edges_from(edges)
        cc = [list(x) for x in nx.connected_components(G)]
        if minFeature > 0:
            cc = list(filter(lambda x: len(x) >= minFeature,cc))
        
        
        for k,ii in enumerate(cc): 
            idx[pos[ii]] = k        
    return(idx)

