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

import subprocess as sp
import sys

import os.path as osp
import os
import re

import argparse as arp
import logging

import datetime

plt.style.use('ggplot')

#TODO : Add docstrings to functions
#TODO : Add logger througout procedure


def _load_file(filename):
    df = pd.read_csv(filename,
                     sep = '\t',
                     index_col = 0,
                     )
    df.columns = ['xcoord'] + df.columns.tolist()[0:-1]
    return df


def _version_control():
    file_dir = osp.dirname(osp.abspath(__file__))
    current_dir = os.getcwd()
    
    commands = [f"cd {file_dir:s}",
                f"git describe --tags",
                f"cd {current_dir:s}"]
    
    
    return_values = []
    for cmd in commands:
        proc = sp.Popen(cmd, shell = True, stdout=sp.PIPE)
        return_values.append(proc.communicate()[0])
    
    version = return_values[1].decode('utf-8').replace('\n','')
    
    return version


def _get_sample_name(filename):
    pattern = '[A-Z]{2}\d{5}[_][A-Z]\d'
    return re.search(pattern, filename)[0]

    
def _get_coordinates(df):
    x_coord = df['xcoord'].values.reshape(-1,1)
    y_coord = df['ycoord'].values.reshape(-1,1)
    
    return np.hstack((x_coord,y_coord))
    

def log_header():
    
    time = str(datetime.datetime.now())
    logger.info(f"time of execution : {time:s}")
    
    try:
        cl = ' '.join(sys.argv)
        logger.info(f"shell command : {cl:s} ")
    except Exception as e:
        logger.error(f"could not capture shell command : {e:s}")
    
    try:
        version = _version_control()
        logger.info(f"git version : {version:s}")
    except Exception as e:
        logger.error(f"could not stat git version : {e:s}")
    
    logger.info(f"software file-path : {osp.abspath(__file__):s}")
    logger.info(f"current working directory : {os.getcwd():s}")

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

def configure_logger(logger):
        
        timestamp = '-'.join(str(datetime.datetime.now()).split(' ')).split('.')[0]
    
        c_handler = logging.StreamHandler()
        f_handler = logging.FileHandler(osp.join(os.getcwd(),''.join(['tumor_extraction-', timestamp, '.log'])))
        
        c_handler.setLevel(logging.INFO)
        f_handler.setLevel(logging.INFO)
        
        c_format = logging.Formatter('%(name)s - %(levelname)s - %(message)s')
        f_format = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        
        c_handler.setFormatter(c_format)
        f_handler.setFormatter(f_format)
        
        logger.addHandler(c_handler)
        logger.addHandler(f_handler)
        

def distance(x,y_vec):
    return np.array([l1(x,y) for y in y_vec])

def tumorScan(data,
              minSpt,
              maxDist,
              minTumor = 5,
              ):
    
    
    def propagate(idx,
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
            
    clusters = np.ones(data.shape[0])*(-1)

    for ii in range(clusters.shape[0]):
        if clusters[ii] == -1:
            clusters = propagate(ii, data, minSpt, maxDist, clusters)
    
    if minTumor > 0:
        for cidx in np.unique(clusters):
            if np.sum(clusters == cidx) < minTumor:
                clusters[clusters == cidx] = -1
    
    return clusters
        
def extractTumors(df,
                  minSpt,
                  maxDist,
                  minTumor,
                  feature = 'tumor',
                  select_for = 'tumor',
                  ):
    
    crd = _get_coordinates(df)
    
    all_cluster_labels = np.ones(crd.shape[0])*(-1)
    tumor_idx = np.where((df[feature] == select_for))
    tumor_cluster_labels = tumorScan(crd[tumor_idx], minSpt=minSpt, maxDist=maxDist, minTumor=minTumor)
    
    all_cluster_labels[tumor_idx] = tumor_cluster_labels
    
    return all_cluster_labels

    
def visualizeTumorSelection(data,
                            labels,
                            sample_name,
                            figsize = (10,5),
                            marker_size = 80,
                            annotations = ['tumor', 'non']):
    
    unique_labels = np.unique(labels)
    annotated_cmap = {annotations[0] : 'red', annotations[1] : 'green'}
    colormap = plt.cm.Dark2(np.arange(unique_labels.shape[0]))
    crd = _get_coordinates(data)

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
    ax[0].set_title(' '.join([sample_name, 'Tumor Separation Result']))
    ax[0].set_ylim([0,35])
    ax[0].set_xlim([0,33])
    
    for annotated_label in annotations:
        pos = (data['tumor'] == annotated_label)
        
        ax[1].scatter(x= crd[pos,0],
                      y = crd[pos,1],
                      c = annotated_cmap[annotated_label],
                      s = marker_size)    
    
    ax[1].set_aspect("equal")
    ax[1].grid(False)
    ax[1].set_title('Tumor vs. non-tumor')
    ax[1].set_ylim([0,35])
    ax[1].set_xlim([0,33])
    
    fig.tight_layout()
    
    return fig, ax


def main(input_name,
         output_name,
         plot,
         save_plot,
         min_spot,
         max_dist,
         min_tumor_spots,
         ):

    if osp.isdir(input_name):
        
        logger.info("reading multiple files from {input_name:s} : MULTIPLE FILE MODE")
        
        all_files = os.listdir(input_name)
        n_files = len(all_files)
        
        logger.info(f"{n_files:d} were found. ")
        
        for (num,single_file) in enumerate(all_files):
            try:
                df = _load_file(osp.join(input_name, single_file))
                
                labels = extractTumors(df, 
                                       minSpt= min_spot,
                                       maxDist = max_dist,
                                       minTumor = min_tumor_spots)
                
                tumor_labels = label_to_tumor_id(labels)
                
                df['tumor_id'] = tumor_labels
                
                if not osp.isdir(output_name):
                    os.mkdir(output_name)
                    
                file_output = osp.join(output_name,
                                  '_'.join([_get_sample_name(single_file), 'tumor_separation.tsv']))
                
                df.to_csv(file_output,
                          sep = '\t', 
                          header = True,
                          index = True)
                
                logger.info(f"successfully processed file {single_file:s}")
            except Exception as e:
                logger.error(f"could not process file {single_file:s} : {e:s}")
            
            if plot:
                try:
                    fig, ax = visualizeTumorSelection(df,
                                                      labels,
                                                      sample_name=_get_sample_name(single_file),
                                                      )
                    plt.show()
                except Exception as e:
                    logger.error(f"could not plot cluster result of {single_file:s}")
            
            if save_plot:
                try:
                    fig, ax = visualizeTumorSelection(df,
                                                      labels,
                                                      sample_name=_get_sample_name(single_file),
                                                      )
                    
                    img_output = osp.join(output_name,
                                          '_'.join([_get_sample_name(single_file), 'tumor_separation.png']),
                                          )
                    fig.savefig(img_output)
                    logger.info(f"successfully saved image of {single_file:s} clustering at : {img_output:s})
                except Exception as e:
                    logger.error(f"could not save image of {single_file:s} : {e:s}")
                    
                    
                
            plt.close('all')
            
    elif osp.isfile(input_name):
        logger.info(f"reading file from {input_name:s} : SINGLE FILE MODE")
        df = _load_file(input_name)
        
        labels = extractTumors(df, 
                               minSpt= min_spot,
                               maxDist = max_dist,
                               minTumor = min_tumor_spots)
            
        tumor_labels = label_to_tumor_id(labels)
        
        df['tumor_id'] = tumor_labels
        
        if osp.isdir(output_name):
            file_output = osp.join(output_name,'_'.join([_get_sample_name(input_name), 'tumor_separation.tsv']))
        
        else:
            file_output = osp.join('_'.join([_get_sample_name(input_name), 'tumor_separation.tsv']))

        df.to_csv(file_output, sep = '\t', header = True, index = True)
        
        if plot:
            try:
                fig, ax = visualizeTumorSelection(df, 
                                                  labels,
                                                  sample_name = _get_sample_name(osp.basename(input_name)),
                                                  )
                plt.show()
            except Exception as e:
                logger.info()
            
        if save_plot:
            try:
                fig, ax = visualizeTumorSelection(df,
                                                  labels,
                                                  sample_name = _get_sample_name(osp.basename(input_name)),
                                                  )
                
                if osp.isdir(output_name):
                    img_output = osp.join(output_name,
                                          '_'.join([_get_sample_name(input_name), 'tumor_separation.png']),
                                          )
                else:
                    img_output = osp.join('_'.join([_get_sample_name(input_name), 'tumor_separation.png']))
                    
                fig.savefig(img_output)
            
            logger.info(f"successfully saved image of {input_name:s} at : {img_output:s}")
        
        except Exception as e:
            logger.error(f"could not save image of {input_name:s} cluster result : {e:s}")
       
        plt.close('all')

#%%
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
    
    
    logger = logging.basicConfig(level = logging.INFO, filename = os.devnull)
    logger = logging.getLogger('te_log')
    
    configure_logger(logger)
    log_header()
    
    arguments = dict(input_name = args.input,
                     output_name = args.output,
                     plot = args.plot,
                     save_plot = args.save_plot,
                     min_spot = args.min_spot,
                     max_dist = args.max_dist,
                     min_tumor_spots = args.min_tumor_spots,
                     )
    
    main(**arguments)

