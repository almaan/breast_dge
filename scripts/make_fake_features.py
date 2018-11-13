#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 12 17:18:31 2018

@author: alma
"""

import pandas as pd
import numpy as np
import os.path as osp

pth = '/home/alma/ST-2018/CNNp/data/pre-tile-extraction'
out_pth = '/home/alma/ST-2018/CNNp/DGE/data/test_data'
sample_name = 'count_data-23209_C1.tsv'
sample_id = (sample_name.split('-')[-1]).split('.')[0]
#%%
df = pd.read_csv(osp.join(pth,sample_name), sep = '\t', index_col = 0)
#%%
coords = np.array([ x.split('x') for x in df.index ])
f1 = np.random.choice(['tumor','non-tumor'], size = coords.shape[0], replace = True, p = [0.7,0.3])
f2 = np.array(['luma'] * coords.shape[0])
f2[f1 == 'non-tumor'] = 'non'

feature_data = dict(x_coord = coords[:,0], y_coord= coords[:,1], f1 = f1, f2 = f2)
fdf = pd.DataFrame(feature_data)
fdf.to_csv(osp.join(out_pth, '-'.join(['feature_data','.'.join([sample_id,'tsv'])])), index = False, sep = '\t')

