#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec  1 14:33:35 2018

@author: alma
"""

import pandas as pd
import os
import os.path as osp
import re
import numpy as np

pth = "/home/alma/ST-2018/CNNp/data/YAO_feature_files"
pth_out = "/home/alma/ST-2018/CNNp/DGE/data/YAO_feature_files"

files = os.listdir(pth)

for filename in files:
    sample_id = re.search(pattern = '[0-9]{5}_\w{2}', string = files[0])[0]
    df = pd.read_csv(osp.join(pth,filename),
                     sep = '\t',
                     index_col = 0,
                     )
    
    df.columns = ['xcoord'] + df.columns.tolist()[0:-1]
    data = dict(xcoord = df['xcoord'].values,ycoord = df['ycoord'].values,tumor = df['tumor'].values)
    arr = np.round(df[['xcoord','ycoord']].values,0).astype(int)

    idx = [str(arr[k,0]) + 'x' + str(arr[k,1]) for k in range(arr.shape[0])]
    
    
    new_df = pd.DataFrame(data,index = idx)
    new_df = new_df.drop_duplicates()
    new_df.to_csv(osp.join(pth_out,filename), sep = '\t', index = True, header = True)