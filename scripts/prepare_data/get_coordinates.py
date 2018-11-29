#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 29 16:05:30 2018

@author: alma
"""

import os
import os.path as osp
import pandas as pd
import re


import argparse as arp


prs = arp.ArgumentParser()

prs.add_argument('-a','--all_spots',required = True, type = str)

prs.add_argument('-t','--tumor_spots', required = True, type = str)

prs.add_argument('-o', '--out_dir', required = True, type = str)

prs.add_argument('-p','--prefix',
                 required = False,
                 type = str,
                 default = 'std_tumor_annotation-BC0',
                 help = ''.join(['prefix to append to experiment id'
                                 'in the output',
                                 'file names. if non provided',
                                 'default will be used.',
                                 ]))

args = prs.parse_args()

pth_all = args.all_spots
pth_tumor = args.tumor_spots
pth_out = args.out_dir

prefix = args.prefix

all_files = os.listdir(pth_all)
tumor_files = os.listdir(pth_tumor)
single_file = all_files[0]

for single_file in all_files:

    all_spots = list((pd.read_csv(osp.join(pth_all,single_file), compression = 'gzip', sep = '\t', nrows = 0, index_col = 0)).columns)
    pid = re.search('\d{4}', single_file)[0]
    tumor_file = list(filter(lambda x: pid in x, tumor_files))[0]
    tumor_spots = list((pd.read_csv(osp.join(pth_tumor,tumor_file), sep = '\t', compression = 'gzip', nrows = 0, index_col = 0).columns))
    labels = [('tumor' if x in tumor_spots else 'non') for x in all_spots]
    
    xcoord = [float(x.split('x')[0]) for x in all_spots]
    ycoord = [float(x.split('x')[1]) for x in all_spots]
    data = dict(xcoord = xcoord, ycoord = ycoord, tumor = labels)
    
    
    df = pd.DataFrame(data = data, index = all_spots)
    df.to_csv(osp.join(pth_out,''.join([prefix,pid,'_X0.tsv'])), sep = '\t', index = True, header = True)
