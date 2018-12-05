#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  5 07:43:07 2018

@author: alma
"""

import sys
import os.path as osp

filename = sys.argv[1]
home = osp.dirname(osp.abspath(filename))
with open(filename, "r+") as fopen:
    result_files = fopen.readlines()

result_files = [x.replace("\n","") for x in result_files]
n_runs = len(result_files)
topgenes = []
for result_single in result_files:
    with open(osp.join(home,result_single)) as fopen:
        tmp = [x.replace("\n","") for x in fopen.readlines()]
        topgenes.append(set(tmp))

n_genes = len(tmp)
inter = topgenes[0]

for ii in range(1,len(topgenes)):
    inter = inter.intersection(topgenes[ii])

inter = list(inter)
n_inter = len(inter)

print(f"Runs : {n_runs:d}  | Number of genes : {n_genes:d}  | Intersection : {n_inter} ")