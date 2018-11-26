#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 20 09:28:23 2018

@author: Alma Anderssson

Program to identify isolated regions (such as tumors) from annotated ST-data. Clusters of spatially coherent
and isolated spots will be marked accordingly to cluster index. The number of clusters
within the sample does not have to be prespecified (and cannot be either).

Each cluster consists of two type of spots, inner spots and boundary spots. A spot is considered an inner spot
if it has at least min_spot neighbours. An inner spot can "propagate" the cluster and will assign all of it's 
neighbors to the same cluster as itself. Boundary spots belong to a tumor but has less than min_spot neighbors
meaning that they cannot propagate a cluster. 

An option to mask clusters below a given size (number of spots) is given, which can be specified with 
the argument     .

Input can be either a directory of multiple files to be processed simultaneously or a single file. For
multiple-file mode enter a directory as output, if directory does not exist it will be created in current
working directory. For single file mode enter a file name as output.

Each file to be processed must have the sample name specified somewhere within the file according to the pattern
"XY#####" where X and Y are arbitrary alphabetical characters within the range [A-Z] (upper or lowercase) and "#"
represents a digit in the range[0-9]. The x-coordinates and y-coordinates should be given in columns named "xcoord"
and "ycoords" respectively. The column containing the annotated feature can be named arbitrarily (default is "tumor"),
and should be passed as an argument (feature) if not "tumor". The annotation of interest (only support) for one as of now
should be passed as an argument if other than default ("tumor".)


Three parameters are used for cluster generation:
    
    max_distance    - the maximum manhattan distance for two spots to be considered neihbours 
    min_spots       - the number of neihbours a spot must have to be considered an inner point of the cluster
    min_total_spots - the minimum number of spots that a cohort of spots must have to be considered a tumor


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


def _load_file(filename):
    """
    Function to load tumor file specified according to the predefined format
    """
    df = pd.read_csv(filename,
                     sep = '\t',
                     index_col = 0,
                     )
    df.columns = ['xcoord'] + df.columns.tolist()[0:-1]
    return df


def _version_control():
    """
    Returns the git-version of the program
    """
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
    """
    
    Get name of a sample based on a filename.
    Uses regex-based patterns. Sample must be
    names in format "XY#####" where represents
    two letters in upper or lowercase and # represents
    a number between 0-9.
    
    """
    pattern = '[A-Z]{2}\d{5}[_][A-Z]\d'
    return re.search(pattern, filename.upper())[0]

    
def _get_coordinates(df):
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
    

def _log_header():
    """
    Adds basic information to top of log-file for reproducibility and
    easy backtracking of action. The elements added are:
        
        - time and date of execution
        - arguments passed (copy of sys.argv)
        - git-version
        - location of software executed
        - current working directory
        
    """
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


def configure_logger(logger, log_name):
        """
        Configures logger. Messages will be passed to STDOUT and a log-file. If
        no log-file name is defined by the user the default name in the form 
        "feature_extraction-[timestamp].log" will be used and saved to the
        current working directory. 
        
        If the file-name provided by the user does not have a ".log" extension
        this will be appended to it.
        
        Arguments:
            logger - logging object to be configured
            log_name - log-file name. If set to False default is used.
        
        """
    
        prefix = 'feature_extraction-'    
    
        if isinstance(log_name, bool):
            timestamp = str(datetime.datetime.now()).replace(':','_').replace(' ','_')
            log_out = ''.join([prefix, timestamp, '.log'])
            
        else:
            if not (log_name.splot('.')[-1] == 'log'):
                log_out = '.'.join(log_name,'log')
            
        
        
        c_handler = logging.StreamHandler()
        f_handler = logging.FileHandler(osp.join(os.getcwd(),log_out))
        
        c_handler.setLevel(logging.INFO)
        f_handler.setLevel(logging.INFO)
        
        c_format = logging.Formatter('%(name)s - %(levelname)s - %(message)s')
        f_format = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        
        c_handler.setFormatter(c_format)
        f_handler.setFormatter(f_format)
        
        logger.addHandler(c_handler)
        logger.addHandler(f_handler)
        
        if isinstance(log_name, bool):
            logger.info(f'using default log-name. Log is saved in : {log_out:s}')
        else:
            logger.info(f'using user-defined log-name. Log is saved in : {log_out:s}')

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
    crd = _get_coordinates(df)
    
    all_cluster_labels = np.ones(crd.shape[0])*(-1)
    tumor_idx = np.where((df[feature] == select_for))
    tumor_cluster_labels = featureScan(crd[tumor_idx], minSpt=minSpt, maxDist=maxDist, minFeature=minFeature)
    
    all_cluster_labels[tumor_idx] = tumor_cluster_labels
    
    return all_cluster_labels

    
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


def main(input_name,
         output_name,
         plot,
         save_plot,
         min_spot,
         max_dist,
         min_total_spots,
         feature,
         select_for,
         ):

    if osp.isdir(input_name):
        
        suffix = 'feature_separation'
        logger.info("reading multiple files from {input_name:s} : MULTIPLE FILE MODE")
        
        all_files = os.listdir(input_name)
        all_files = list(filter(lambda x: x.split('.')[-1] == 'tsv',all_files))
        n_files = len(all_files)
        
        logger.info(f"{n_files:d} files were found. ")
        
        for (num,single_file) in enumerate(all_files):
            try:
                df = _load_file(osp.join(input_name, single_file))
                
                labels = extractFeatures(df, 
                                         minSpt= min_spot,
                                         maxDist = max_dist,
                                         minFeature = min_total_spots,
                                         feature = feature,
                                         select_for = select_for,
                                         )
                
                feature_labels = label_to_feature_id(labels)
                
                df['feature_id'] = feature_labels
                
                if not osp.isdir(output_name):
                    os.mkdir(osp.join(os.getcwd(),output_name))
                    
                file_output = osp.join(output_name,
                                  '_'.join([_get_sample_name(single_file), '.'.join([suffix,'tsv'])]))
                
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
                                                      feature=feature,
                                                      annotation=select_for,
                                                      )
                    plt.show()
                except Exception as e:
                    logger.error(f"could not plot cluster result of {single_file:s} : {e:s}")
            
            if save_plot:
                try:
                    fig, ax = visualizeTumorSelection(df,
                                                      labels,
                                                      sample_name=_get_sample_name(single_file),
                                                      feature=feature,
                                                      annotation=select_for,
                                                      )
                    
                    img_output = osp.join(output_name,
                                          '_'.join([_get_sample_name(single_file),'.'.join([suffix,'png'])]))

                    fig.savefig(img_output)
                    logger.info(f"successfully saved image of {single_file:s} clustering at : {img_output:s}")
                except Exception as e:
                    logger.error(f"could not save image of {single_file:s} : {e:s}")
                    
                    
                
            plt.close('all')
            
    elif osp.isfile(input_name):
        logger.info(f"reading file from {input_name:s} : SINGLE FILE MODE")
        df = _load_file(input_name)
        
        labels = extractFeatures(df, 
                                 minSpt= min_spot,
                                 maxDist = max_dist,
                                 minFeature = min_total_spots,
                                 feature = feature,
                                 select_for = select_for,
                                 )
            
        feature_labels = label_to_feature_id(labels)
        
        df['feature_id'] = feature_labels
        
        if osp.isdir(output_name):
            file_output = osp.join(output_name,'_'.join([_get_sample_name(input_name), '.'.join([suffix,'tsv'])]))
        
        else:
            file_output = osp.join('_'.join([_get_sample_name(input_name), '.'.join([suffix,'tsv'])]))

        df.to_csv(file_output, sep = '\t', header = True, index = True)
        
        if plot:
            try:
                fig, ax = visualizeTumorSelection(df, 
                                                  labels,
                                                  sample_name = _get_sample_name(osp.basename(input_name)),
                                                  feature = feature,
                                                  annotation=select_for,
                                                  )
                plt.show()
            except Exception as e:
                logger.info(f"could not generate plot : {e:s}")
            
        if save_plot:
            try:
                fig, ax = visualizeTumorSelection(df,
                                                  labels,
                                                  sample_name = _get_sample_name(osp.basename(input_name)),
                                                  feature=feature,
                                                  annotation=select_for,
                                                  )
                
                if osp.isdir(output_name):
                    img_output = osp.join(output_name,
                                          '_'.join([_get_sample_name(input_name), '.'.join([suffix,'png'])]))
                else:
                    img_output = osp.join('_'.join([_get_sample_name(input_name), '.'.join([suffix,'png'])]))
                    
                fig.savefig(img_output)
            
                logger.info(f"successfully saved image of {input_name:s} at : {img_output:s}")
        
            except Exception as e:
                logger.error(f"could not save image of {input_name:s} cluster result : {e:s}")

        if plot or save_plot:
            plt.close('all')

if __name__ == '__main__':

    
    prs = arp.ArgumentParser()

    prs.add_argument('-i','--input',
                     required = True,
                     type = str,
                     help = ' '.join(['input file or directory. If single filename is',
                                     'provided then all ".tsv" files within that directory',
                                     'will be processed. Names of files must contain sample id',
                                     'in format "XY#####" where X and Y are arbitrary alphabetic',
                                     'characters and "#" are digits in the range [0-9]',
                                     ]),
                     )
    
    prs.add_argument('-o', '--output',
                     required = True,
                     type = str,
                     help = ' '.join(['name of output. If input is a single file then output will',
                                      'be interpreted as a single-file name. If input is a directory',
                                      'output will be interpreted as a output directory where the result',
                                      'will be saved to. Each filename will be given a standard name based',
                                      'on sample id.',
                                      ])
                     )
    
    prs.add_argument('-s', '--min_spot',
                     required = False,
                     type = int,
                     default = 3,
                     help = ' '.join(['minimum number of neighbouring spots required for a point',
                                      'to be considered an inner point.',
                                      ]),
                     )
    
    prs.add_argument('-d','--max_dist',
                     required = False,
                     type = float,
                     default = 2,
                     help = ' '.join(['maximum distance between two spots',
                                      'for them to be considered as neighbours',
                                      ]),
                     )
    
    prs.add_argument('-t', '--min_total_spots',
                     required = False,
                     type = int,
                     default = -1,
                     help = ' '.join(['minimum number of spots within a cluster',
                                      'for the cluster to be considered a proper cluster',
                                      'i.e. not only noise.'
                                      ]),
                     )
    
    prs.add_argument('-p', '--plot',
                     required = False,
                     default = False,
                     action = 'store_true',
                     help = 'plot result and display result',
                     )
    
    prs.add_argument('-ps', '--save_plot',
                     required = False,
                     default = False,
                     action = 'store_true',
                     help = 'save plots of result. Does not display result',
                     )
    
    prs.add_argument('-l','--logname',
                     required = False,
                     default = False,
                     help = ' '.join(['name of log-file.',
                                      'if non is provided a default log-file name',
                                      'will be used.',
                                      ]),
                     )

    prs.add_argument('-ft','--feature',
                 required = False,
                 default = 'tumor',
                 help = ' '.join(['name of column with annotated features.',
                                  'if not specified "tumor" is used default',
                                  ]),
                 )
    
    prs.add_argument('-se','--select_for',
                 required = False,
                 default = 'tumor',
                 help = ' '.join(['feature of interest to be clustered',
                                  'if non is specified "tumor" is used',
                                  ]),
                 )
    
    args = prs.parse_args()
    
    
    logger = logging.basicConfig(level = logging.INFO, filename = os.devnull)
    logger = logging.getLogger('feature_extraction_log')
    
    configure_logger(logger, log_name = args.logname)
    _log_header()
    
    arguments = dict(input_name = args.input,
                     output_name = args.output,
                     plot = args.plot,
                     save_plot = args.save_plot,
                     min_spot = args.min_spot,
                     max_dist = args.max_dist,
                     min_total_spots = args.min_total_spots,
                     feature = args.feature,
                     select_for = args.select_for,
                     )
    
    main(**arguments)

