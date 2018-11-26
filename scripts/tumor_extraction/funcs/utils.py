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


def load_file(filename):
    """
    Function to load tumor file specified according to the predefined format
    """
    df = pd.read_csv(filename,
                     sep = '\t',
                     index_col = 0,
                     )
    df.columns = ['xcoord'] + df.columns.tolist()[0:-1]
    return df


def get_sample_name(filename):
    """
    
    Get name of a sample based on a filename.
    Uses regex-based patterns. Sample must be
    names in format "XY#####" where represents
    two letters in upper or lowercase and # represents
    a number between 0-9.
    
    """
    pattern = '[A-Z]{2}\d{5}[_][A-Z]\d'
    return re.search(pattern, filename.upper())[0]

    
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
    x_coord = df['xcoord'].values.reshape(-1,1)
    y_coord = df['ycoord'].values.reshape(-1,1)
    
    crd = np.hstack((x_coord,y_coord))
    
    return crd
    

def label_to_feature_id(idx_list,
                   prefix = '',
                   non_label = 'non',
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