#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  6 13:08:48 2019

@author: alma

"""

import argparse as arp

import os.path as osp
import os 

import pandas as pd
import numpy as np

import matplotlib.pyplot as plt

from funcs.pairing import pair_multiple_sections as pms


# Parser --------------------------


prs = arp.ArgumentParser()

prs.add_argument('-f','--feature_file',
                 required = True,
                 nargs = '*',
                 help = ''.join(['',
                                 ])
                 )

prs.add_argument('-o','--output',
                 required = True,
                 default = '',
                 help = ''.join(['',
                                 ])
                 ) 


prs.add_argument('-c','--colname',
                 required = False,
                 default = 'tumor_id',
                 help = ''.join(['',
                                 ])
                 )    

prs.add_argument('-is','--imgsave',
                 required = False,
                 action = 'store_true',
                 default = False,
                 help = ''.join(['',
                                 ])
                 )    


prs.add_argument('-e','--exclude',
                 action = 'store_const',
                 default = np.inf,
                 const = -1)

prs.add_argument('-x','--xlab',
                 default = 'xcoord')

prs.add_argument('-y','--ylab',
                 default = 'ycoord')



args = prs.parse_args()

# Main ------------------------

read_feature_file = lambda x: pd.read_csv(x, header = 0, index_col = 0, sep = "\t")

feature_files = [read_feature_file(f) for f in args.feature_file]

n_members = [np.unique(f[args.colname][f[args.colname] != args.exclude]).shape[0] for f in feature_files]

if any([x > 1 for x in n_members]):
    new_feature_files = pms(*feature_files, colname = args.colname,xlab =args.xlab,ylab = args.ylab)

else:
    new_feature_files = feature_files

if args.output:
    onames = [osp.join(args.output,
                       ''.join(['mapped_id.',osp.basename(f)])) for f in args.feature_file]
else:
    onames = [osp.join(osp.dirname(f),
                       ''.join(['mapped_id.',osp.basename(f)])) for f in args.feature_file]

for filenum in range(len(new_feature_files)):
    feature_files[filenum].to_csv(onames[filenum], sep = '\t',header = True, index = True)
    
    if args.imgsave:
        plt.clf()
        plt.figure(1)
        for lab in np.unique(new_feature_files[filenum][args.colname]):
            pos = new_feature_files[filenum][args.colname] == lab
            plt.scatter(new_feature_files[filenum][args.xlab][pos],
                        new_feature_files[filenum][args.ylab][pos])
        
        plt.savefig(onames[filenum].replace('.tsv','.png'))