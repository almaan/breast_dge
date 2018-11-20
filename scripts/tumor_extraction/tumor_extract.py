#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 20 09:28:23 2018

@author: alma
"""

import pandas as pd
import numpy as np
from scipy.spatial.distance import cityblock as l1


import matplotlib.pyplot as plt

import os.path as osp
import os
import re

import argparse as arp
#import logger

plt.style.use('ggplot')

#TODO : Add docstrings to functions

def _load_file(filename):
    df = pd.read_csv(filename,
                     sep = '\t',
                     index_col = 0,
                     )
    df.columns = ['xcoord'] + df.columns.tolist()[0:-1]
    return df

def _get_sample_name(filename):
    pattern = '[A-Z]{2}\d{5}[_][A-Z]\d'
    return re.search(pattern, filename)[0]

    
def get_coordinates(df):
    x_coord = df['xcoord'].values.reshape(-1,1)
    y_coord = df['ycoord'].values.reshape(-1,1)
    
    return np.hstack((x_coord,y_coord))
    

def label_to_tumor_id(idx_list,
                   prefix = 'T',
                   non_label = 'non',
                   ):
    
    tumor_labels = []
    for label in idx_list:
        if label >= 0:
            tumor_labels.append(''.join([prefix,str(int(label))]))
        else:
            tumor_labels.append(non_label)
    
    return tumor_labels

def plot_section(dataframe):
    x_coord = dataframe['xcoord'].values.reshape(-1,1)
    y_coord = dataframe['ycoord'].values.reshape(-1,1)
    cmap = dict(tumor = 'red', non = 'green')
    fig, ax = plt.subplots(1,1)
    
    for idx in ['tumor','non']:
        pos = (dataframe['tumor'] == idx)
        ax.scatter(x_coord[pos],
                   y_coord[pos],
                   c = cmap[idx],
                   marker = 'o',
                   )
    
    fig.show()

def distance(x,y_vec):
    return np.array([l1(x,y) for y in y_vec])

def TumorScan(data,
              minSpt,
              maxDist,
              minTumor = 5,
              ):
    
    
    def Propagate(idx,
                  data,
                  minSpt,
                  maxDist,
                  ori_clusters,
                  ):
        
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
            
    rnd = np.arange(data.shape[0])
    np.random.shuffle(rnd)
    data = data[rnd, :]
    
    clusters = np.ones(data.shape[0])*(-1)

    for ii in range(clusters.shape[0]):
        if clusters[ii] == -1:
            clusters = Propagate(ii, data, minSpt, maxDist, clusters)
    
    clusters = clusters[np.argsort(rnd)]
    
    if minTumor > 0:
        for cidx in np.unique(clusters):
            if np.sum(clusters == cidx) < minTumor:
                clusters[clusters == cidx] = -1
    
    return clusters
        
def ExtractTumors(df,
                  minSpt,
                  maxDist,
                  minTumor,
                  feature = 'tumor',
                  select_for = 'tumor',
                  ):
    
    crd = get_coordinates(df)
    
    all_cluster_labels = np.ones(crd.shape[0])*(-1)
    tumor_idx = np.where((df[feature] == select_for))
    tumor_cluster_labels = TumorScan(crd[tumor_idx], minSpt=minSpt, maxDist=maxDist, minTumor=minTumor)
    
    all_cluster_labels[tumor_idx] = tumor_cluster_labels
    
    return all_cluster_labels

    
def VisualizeTumorSelection(data,
                            labels,
                            sample_name,
                            figsize = (10,5),
                            marker_size = 80,
                            annotations = ['tumor', 'non']):
    
    unique_labels = np.unique(labels)
    annotated_cmap = {annotations[0] : 'red', annotations[1] : 'green'}
    colormap = plt.cm.Dark2(np.arange(unique_labels.shape[0]))
    
    fig, ax = plt.subplots(1,2, figsize = figsize)
    
    for (color, label_idx) in zip(colormap,np.unique(labels)):
            
            pos = (labels == label_idx)
            if label_idx >= 0:
                ax[0].scatter(x = data[pos,0],
                              y = data[pos,1],
                              s = marker_size,
                              edgecolors = 'black',
                              c = color)
            else:
                ax[0].scatter(x = data[pos,0],
                              y = data[pos,1],
                              c = 'k',
                              alpha = 0.3,
                              s = marker_size)
    
    ax[0].set_aspect("equal")
    ax[0].grid(False)
    ax[0].set_title(' '.join([sample_name, 'Tumor Separation Result']))
    ax[0].set_ylim([0,35])
    ax[0].set_xlim([0,33])
    
    for annotated_label in annotations:
        pos = (df['tumor'] == annotated_label)
        ax[1].scatter(crd[pos,0],crd[pos,1], c = annotated_cmap[annotated_label], s = marker_size)    
    
    ax[1].set_aspect("equal")
    ax[1].grid(False)
    ax[1].set_title('Tumor vs. non-tumor')
    ax[1].set_ylim([0,35])
    ax[1].set_xlim([0,33])
    
    fig.tight_layout()
    
    return fig, ax




if __name__ == '__main__':

# TODO : Add help information to all arguments
    
    prs = arp.ArgumentParser()

    prs.add_argument('-i','--input',
                     required = True,
                     type = str,
                     help = '',
                     )
    
    prs.add_argument('-o', '--output',
                     required = True,
                     type = str,
                     help = '',
                     )
    
    prs.add_argument('-s', '--min_spot',
                     required = False,
                     type = int,
                     default = 3,
                     help = '',
                     )
    
    prs.add_argument('-d','--max_dist',
                     required = False,
                     type = float,
                     default = 2,
                     help = '',
                     )
    
    prs.add_argument('-t', '--min_tumor_spots',
                     required = False,
                     type = int,
                     default = -1,
                     help = '',
                     )
    
    prs.add_argument('-p', '--plot',
                     required = False,
                     default = False,
                     action = 'store_true',
                     help = '',
                     )
    
    prs.add_argument('-ps', '--save_plot',
                     required = False,
                     default = False,
                     action = 'store_true',
                     help = '',
                     )
    
    args = prs.parse_args()
    
    if osp.isdir(args.input):
        all_files = os.listdir(args.input)
        
        for single_file in all_files:
            
            df = _load_file(osp.join(args.input, single_file))
            labels = ExtractTumors(df, 
                                   minSpt= args.min_spot,
                                   maxDist = args.max_dist,
                                   minTumor = args.min_tumor_spots)
            
            tumor_labels = label_to_tumor_id(labels)
            
            df['tumor_id'] = tumor_labels
            
            if not osp.isdir(args.output):
                os.mkdir(args.output)
                
            output = osp.join(args.output,'_'.join([_get_sample_name(single_file), 'tumor_separation.tsv']))
            df.to_csv(output, sep = '\t', header = True, index = True)
            
            if args.plot:
                crd = get_coordinates(df)
                fig, ax = VisualizeTumorSelection(crd, labels, sample_name=_get_sample_name(single_file))
                plt.show()
            
            if args.save_plot:
                crd = get_coordinates(df)
                fig, ax = VisualizeTumorSelection(crd, labels, sample_name=_get_sample_name(single_file))
                
                img_output = osp.join(args.output,'_'.join([_get_sample_name(single_file), 'tumor_separation.png']))
                fig.savefig(img_output)
                
            plt.close('all')
            
    elif osp.isfile(args.input):
        
        df = _load_file(args.input)
        labels = ExtractTumors(df, 
                                   minSpt= args.min_spot,
                                   maxDist = args.max_dist,
                                   minTumor = args.min_tumor_spots)
            
        tumor_labels = label_to_tumor_id(labels)
        df['tumor_id'] = tumor_labels
        
        output = osp.join('_'.join([_get_sample_name(args.input), 'tumor_separation.tsv']))
        df.to_csv(output, sep = '\t', header = True, index = True)
        
        if args.plot:
                crd = get_coordinates(df)
                fig, ax = VisualizeTumorSelection(crd, labels, sample_name = _get_sample_name(osp.basename(args.input)))
                plt.show()
            
        if args.save_plot:
            crd = get_coordinates(df)
            fig, ax = VisualizeTumorSelection(crd, labels, sample_name = _get_sample_name(osp.basename(args.input)))
            
            img_output = osp.join('_'.join([_get_sample_name(args.input), 'tumor_separation.png']))
            fig.savefig(img_output)
        
        plt.close('all')
        

