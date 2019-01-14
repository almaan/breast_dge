#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Tue Nov 20 09:28:23 2018

@author: Alma Anderssson

Script to identify isolated regions (such as tumors) from annotated ST-data. Clusters of spatially coherent
and isolated spots will be marked accordingly to cluster index. The number of clusters
within the sample does not have to be prespecified (and cannot be either).

Three parameters are used for cluster generation:
    
    max_distance    - the maximum manhattan distance for two spots to be considered neihbours 
    min_total_spots - the minimum number of spots that a cohort of spots must have to be considered a tumor
    norm            - minkowski norm to be used in distance esimation 

"""

import os
import logging

logger = logging.basicConfig(level = logging.INFO, filename = os.devnull)
logger = logging.getLogger('feature_extraction_log')



import matplotlib.pyplot as plt

import sys
import time

import os.path as osp
import numpy as np
import argparse as arp

from funcs.log_configs import log_header, configure_logger
from funcs.utils import *
from clustering.connected import *
from clustering.visual import *

def main(input_names,
         output_name,
         plot,
         save_plot,
         max_dist,
         min_total_spots,
         p_norm,
         feature,
         select_for,
         tagged,
         ):

    #if p_norm is negative set to infinity
    p_norm = (np.inf if p_norm == -1 else p_norm)    
    
    for input_name in input_names:    
        logger.info(f"reading file from {input_name:s}")
        #load and prepare data
        df = load_file(input_name)
        #get coordinates from feature file
        crd = get_coordinates(df)
        
        #get indices of isolated feature regions
        try:
            labels = connected_graphs_extraction(crd, 
                                                 maxDist = max_dist,
                                                 minFeature = min_total_spots,
                                                 feature_vector = df[feature].values,
                                                 select_for = select_for,
                                                 p_norm = p_norm,
                                                 )
            
        except Exception as e:
            logger.error(f"Could not cluster sample. {e}")
            sys.exit(0)
        
        feature_labels = label_to_feature_id(labels)
        
        df['_'.join([feature,'id'])] = feature_labels
        
        #if regex patterns for id and replicate are provided, add as columns to output
        if tagged[0]:
            sid, replicate = get_sample_tags(input_name, pattern = tagged, join = False)  
            df['id'] = sid
            df['replicate'] = replicate
        
        
        #suffic of output files
        suffix = 'feature_separation'
        #basename for output files to be generated
        ofile_base = '_'.join([get_sample_tags(input_name, pattern = tagged, join = True), suffix])
        #save results
        if not osp.isdir(output_name):
            os.mkdir(output_name)
                
        file_output = osp.join(output_name,'.'.join([ofile_base,'tsv']))
        
        #save result as modified feature-file using newly defined rownames and header
        try:
            df.to_csv(file_output, sep = '\t', header = True, index = True)
        
        except: 
            logger.error(f'could not save file : {file_output:s}')
            
        #if flag for interactinve plotting is included
        if plot:
            try:
                fig, ax = visualizeTumorSelection(df, 
                                                  labels,
                                                  sample_name = get_sample_tags(osp.basename(input_name),pattern=tagged,join=True),
                                                  feature = feature,
                                                  annotation=select_for,
                                                  )
                plt.show()
            except Exception as e:
                
                logger.info(f"could not generate plot : {e}")
        #if flag to save plot-result is included. Note how this does _not_ display the graph but only saves it    
        if save_plot:
            try:
                fig, ax = visualizeTumorSelection(df,
                                                  labels,
                                                  sample_name = get_sample_tags(osp.basename(input_name),pattern = tagged, join=True),
                                                  feature=feature,
                                                  annotation=select_for,
                                                  )
                
                #if output is given as directory
                img_output = osp.join(output_name,'.'.join([ofile_base,'png']))
                    
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
                     nargs = '*',
                     help = ' '.join(['input file name. If multiple file names',
                                      'are given each of them will be read separately',
                                      'and placed in a common output folder.'
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
    
    prs.add_argument('-tg','--tagged',
                     required = False,
                     default = ['',''],
                     nargs = '*',
                     help = ' '.join(['use in order to extract sample id',
                                      'and replicate. Takes two arguments',
                                      '(1) regex for sample id and (2) regex',
                                      'for replicate. If only flag is included',
                                      'default patterns will be used, being : \n',
                                      'sample id: [A-Z]{2}[0-9]{5} \n',
                                      'replicate: (?<=_)([A-Z]{1}[0-9]{1})'
                                    ]))
    
    prs.add_argument('-n', '--norm',
                     required = False,
                     type = int,
                     default = -1,
                     help = ' '.join(['which form of Minkowski norm that',
                                      'should be used when computing distance',
                                      'to neighbours. Default is -1, here being',
                                      'representing the inf-norm.'
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
    
    #configure logger
    configure_logger(logger, 
                     log_name = '_'.join([(args.logname if args.logname else 'feature_extraction'),time.strftime('%d%m%y_%H:%M')]), 
                     output_dir = args.output,
                     )
    
    log_header(logger)
    
    if (args.tagged != [] and len(args.tagged) != 2):
        logger.error('Provide two patterns to the "tagged" flag')
        sys.exit(0)
    
    #if tagged flag is used without any arguments. Mimics "const" argument when nargs = '?'    
    args.tagged = (['[A-Z]{2}[0-9]{5}','(?<=_)([A-Z]{1}[0-9]{1})'] if args.tagged == [] else args.tagged)
    
    arguments = dict(input_names = args.input,
                     output_name = args.output,
                     plot = args.plot,
                     save_plot = args.save_plot,
                     max_dist = args.max_dist,
                     min_total_spots = args.min_total_spots,
                     p_norm = args.norm,
                     feature = args.feature,
                     select_for = args.select_for,
                     tagged = args.tagged
                     )
    try:
        main(**arguments)
            
    except KeyboardInterrupt:
        from time import sleep
        logger.info(f"Keyboard Interrupt. Will wait for potential writing of files to terminate")
        sleep(2)
        sys.exit(0)
