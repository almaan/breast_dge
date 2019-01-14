#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 26 17:01:28 2018

@author: alma
"""

import logging
logger = logging.getLogger('feature_extraction_log')

import pandas as pd
import numpy as np
import re
from scipy.spatial.distance import cityblock as l1
import sys
import os.path as osp

def match_rownames(data,sample_id):
    """
    generates new rownames which useful for downstream processes. Rownames are on form
    sample_id_[x_coord]x[y_coord]
    
    """
    idx = np.apply_along_axis(func1d = lambda x: '_'.join([sample_id,'x'.join([str(int(round(x[0],0))),str(int(round(x[1],0)))])]),
                              axis = 1,
                              arr = data[['xcoord','ycoord']].values)
    return idx

def load_file(filename):
    """
    Function to load tumor file specified according to the predefined format
    """
    try:
        df = pd.read_csv(filename,
                         sep = '\t',
                         index_col = 0,
                         )
        return df
    except Exception as e:
        logger.error(f'Could not load file : {e}')        
        sys.exit(0)

def get_sample_tags(filename,pattern,join=False):
    """
    
    Get name of a sample based on a filename.
    If a regex pattern for sample and replicate
    name is provided, the filename 
    
    """
    res = osp.basename(filename).split('.')[0:-1]
    
    if pattern[0]:
        rep = re.search(pattern[1],filename)
        samp = re.search(pattern[0],filename)
        res = ([samp.group(0),rep.group(0)] if (samp and res) else res)
    
    if len(res) > 1 and join:
        res = '_'.join(res)
    
    elif len(res) < 2:
        res = res[0]
        
    return res

    
def get_coordinates(df):
    """
    Extract coordinates from dataframe and return as a
    Nx2 numpy array, x-coordinates in first column and
    y-coordinates in second.
    
    Arguments:
        df - pandas datafram
    Returns:
        crd - numpy array of size Nx2
    
    """
    assert 'xcoord' in df.columns and 'ycoord' in df.columns, "Feature file not properly formatted"
    
    x_coord = df['xcoord'].values.reshape(-1,1)
    y_coord = df['ycoord'].values.reshape(-1,1)
    
    crd = np.hstack((x_coord,y_coord))
    
    return crd
    

def label_to_feature_id(idx_list,
                   prefix = '',
                   non_label = '-1',
                   ):
    
    """
    provided an list of feature indices where the indices are in the range [0,infty]
    and non-feature items are indicated with -1; a label with specified notation
    of feature number and non-feature labeling according to
    
    feature : [prefix][idx]
    non-feature : [non_label]
    
    Arguments:
        idx_list : list or array of indices to be converted to labels
        prefix : prefix to add before index of feature. Default is NULL aka. same number as index.
        non_label : label to mark "non-features". Default is "non"
    
    Returns:
        feature_labels : list of modified tumor labels
    
    """
    
    feature_labels = []
    for label in idx_list:
        if label >= 0:
            feature_labels.append(''.join([prefix,str(int(label))]))
        else:
            feature_labels.append(non_label)
    
    return feature_labels

def distance(x,y_vec):
    """
    computes the manhattan distance between a point x and an array/list
    of points y_vec.
    
    Arguments:
        x     - point to compute distances from
        y_vec - list/array of points to compute distances from x to
    
    Return:
        dist  - array of distances between x and all points in y
        
    """
    dist = np.array([l1(x,y) for y in y_vec])
    
    return dist