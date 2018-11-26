#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 26 17:10:03 2018

@author: alma
"""
import logging
logger = logging.getLogger('feature_extraction_log')

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from funcs.utils import *

plt.style.use('ggplot')


def visualizeTumorSelection(data,
                            labels,
                            sample_name,
                            feature,
                            annotation,
                            figsize = (10,5),
                            marker_size = 80,
                            ):
    """
    Function to visualize the clustering result using matplotlib. Generates
    a figure with two subplots; the first one illustrating the clusters and the
    second illustrating the original annotation.
    
    Arguments:
        data - pandas dataframe containing data to be visualized
        labels - feature labels (cluster indicies) corresponding to data object
        feature - name of column where annotated features are listed
        annotation - feature of interst (to be contrasted)
        figsize - size of plot, given in hundreds of pixels. Default is (10,5)
        marker_size - size of marker to be used in scatter-plot. Default is 80
    
    Returns:
        fig - matplotlib figure containing the axes element of the two subplots
        ax - list of 2 axes objects representing the two subplots
    
    """
    
    unique_labels = np.unique(labels)
    colormap = plt.cm.Dark2(np.arange(unique_labels.shape[0]))
    crd = get_coordinates(data)

    fig, ax = plt.subplots(1,2, figsize = figsize)
    
    for (color, label_idx) in zip(colormap,np.unique(labels)):
            
            pos = (labels == label_idx)
            if label_idx >= 0:
                ax[0].scatter(x = crd[pos,0],
                              y = crd[pos,1],
                              s = marker_size,
                              edgecolors = 'black',
                              c = color)
            else:
                ax[0].scatter(x = crd[pos,0],
                              y = crd[pos,1],
                              c = 'k',
                              alpha = 0.3,
                              s = marker_size)
    
    ax[0].set_aspect("equal")
    ax[0].grid(False)
    ax[0].set_title(f"{sample_name:s} {annotation:s} sepration result")
    ax[0].set_ylim([0,35])
    ax[0].set_xlim([0,33])
    
    for annotated_label in  np.unique(data[feature].values):
        pos = (data[feature] == annotated_label)
        spot_color = ('red' if annotated_label == annotation else 'green')
        
        ax[1].scatter(x= crd[pos,0],
                      y = crd[pos,1],
                      c = spot_color,
                      s = marker_size)    
    
    ax[1].set_aspect("equal")
    ax[1].grid(False)
    ax[1].set_title(f"{annotation:s} vs. non-{annotation:s}")
    ax[1].set_ylim([0,35])
    ax[1].set_xlim([0,33])
    
    fig.tight_layout()
    
    return fig, ax