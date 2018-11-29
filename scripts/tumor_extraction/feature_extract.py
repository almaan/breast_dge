#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Tue Nov 20 09:28:23 2018

@author: Alma Anderssson

Script to identify isolated regions (such as tumors) from annotated ST-data. Clusters of spatially coherent
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
and "ycoord" respectively. The column containing the annotated feature can be named arbitrarily (default is "tumor"),
and should be passed as an argument (feature) if not "tumor". The annotation of interest (only support) for one as of now
should be passed as an argument if other than default ("tumor".)


Three parameters are used for cluster generation:
    
    max_distance    - the maximum manhattan distance for two spots to be considered neihbours 
    min_spots       - the number of neihbours a spot must have to be considered an inner point of the cluster
    min_total_spots - the minimum number of spots that a cohort of spots must have to be considered a tumor


"""

import os
import logging

logger = logging.basicConfig(level = logging.INFO, filename = os.devnull)
logger = logging.getLogger('feature_extraction_log')



import matplotlib.pyplot as plt

import sys

import os.path as osp

import argparse as arp

from funcs.log_configs import log_header, configure_logger
from funcs.utils import *
from clustering.cluster import *
from clustering.visual import *

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
        logger.info(f"reading multiple files from {input_name:s} : MULTIPLE FILE MODE")
        
        all_files = os.listdir(input_name)
        all_files = list(filter(lambda x: x.split('.')[-1] == 'tsv',all_files))
        n_files = len(all_files)
        
        logger.info(f"{n_files:d} files were found. ")
        
        for (num,single_file) in enumerate(all_files):
            try:
                print(get_sample_name(single_file))
                df = load_file(osp.join(input_name, single_file))
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
                                  '_'.join([get_sample_name(single_file), '.'.join([suffix,'tsv'])]))
                
                df.to_csv(file_output,
                          sep = '\t', 
                          header = True,
                          index = True)
                
                logger.info(f"successfully processed file {single_file:s}")
            except Exception as e:
                logger.error(f"could not process file {single_file:s} : {e}")
            
            if plot:
                try:
                    fig, ax = visualizeTumorSelection(df,
                                                      labels,
                                                      sample_name=get_sample_name(single_file),
                                                      feature=feature,
                                                      annotation=select_for,
                                                      )
                    plt.show()
                except Exception as e:
                    logger.error(f"could not plot cluster result of {single_file:s} : {e}")
            
            if save_plot:
                try:
                    fig, ax = visualizeTumorSelection(df,
                                                      labels,
                                                      sample_name=get_sample_name(single_file),
                                                      feature=feature,
                                                      annotation=select_for,
                                                      )
                    
                    img_output = osp.join(output_name,
                                          '_'.join([get_sample_name(single_file),'.'.join([suffix,'png'])]))

                    fig.savefig(img_output)
                    logger.info(f"successfully saved image of {single_file:s} clustering at : {img_output:s}")
                except Exception as e:
                    logger.error(f"could not save image of {single_file:s} : {e}")
                    
                    
                
            plt.close('all')
            
    elif osp.isfile(input_name):
        logger.info(f"reading file from {input_name:s} : SINGLE FILE MODE")
        df = load_file(input_name)
        
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
            file_output = osp.join(output_name,'_'.join([get_sample_name(input_name), '.'.join([suffix,'tsv'])]))
        
        else:
            file_output = osp.join('_'.join([get_sample_name(input_name), '.'.join([suffix,'tsv'])]))

        df.to_csv(file_output, sep = '\t', header = True, index = True)
        
        if plot:
            try:
                fig, ax = visualizeTumorSelection(df, 
                                                  labels,
                                                  sample_name = get_sample_name(osp.basename(input_name)),
                                                  feature = feature,
                                                  annotation=select_for,
                                                  )
                plt.show()
            except Exception as e:
                logger.info(f"could not generate plot : {e}")
            
        if save_plot:
            try:
                fig, ax = visualizeTumorSelection(df,
                                                  labels,
                                                  sample_name = get_sample_name(osp.basename(input_name)),
                                                  feature=feature,
                                                  annotation=select_for,
                                                  )
                
                if osp.isdir(output_name):
                    img_output = osp.join(output_name,
                                          '_'.join([get_sample_name(input_name), '.'.join([suffix,'png'])]))
                else:
                    img_output = osp.join('_'.join([get_sample_name(input_name), '.'.join([suffix,'png'])]))
                    
                fig.savefig(img_output)
            
                logger.info(f"successfully saved image of {input_name:s} at : {img_output:s}")
        
            except Exception as e:
                logger.error(f"could not save image of {input_name:s} cluster result : {e}")

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
    
    configure_logger(logger, log_name = args.logname)
    log_header(logger)
    
    try:
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
            
    except KeyboardInterrupt:
        from time import sleep
        logger.info(f"Keyboard Interrupt. Will wait for potential writing of files to terminate")
        sleep(3)
        sys.exit(0)
