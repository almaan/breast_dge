#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 23 17:58:19 2018

@author: alma
"""

import PIL
import PIL.Image as IMG 

import os.path as osp
import os

import numpy as np

pth = "/home/alma/ST-2018/CNNp/DGE/scripts/tumor_extraction/output_results2"

# TODO : Make this clean

def get_format(N,w,h, desired_format):
    r1 = desired_format[0]/desired_format[1]
    
    n_row = np.sqrt(N / r1)
    n_col = r1 * n_row
    
    
    if np.floor(n_row)*np.ceil(n_col) >= N:
        n_row = np.floor(n_row)
        n_col = np.ceil(n_col)
    else:
        n_col = np.ceil(n_col)
        n_row = np.ceil(n_row)
        
    return (n_row.astype(int), n_col.astype(int))        
    
    

file_extension = ".png"
all_images = list(filter(lambda x: file_extension in x, os.listdir(pth)))
all_images = sorted(all_images)
n_images = len(all_images)


desired_format = (1024, 768)

img = IMG.open(osp.join(pth, all_images[0]))
width, height = img.size
ratio = (height / width )
new_width = 800
new_height = np.ceil(new_width * ratio).astype(int)



grid_size = get_format(n_images, width, height, desired_format)

k = 0
background = IMG.new('RGBA',(new_width*(grid_size[0]+1), new_height*(grid_size[1]+1)), (255, 255, 255, 255))

for jj in range(0,int(grid_size[1]+1)):
    for ii in range(0,int(grid_size[0]+1)):
        if k < n_images:
            img = IMG.open(osp.join(pth, all_images[k]))
            img = img.resize((new_width,new_height),IMG.ANTIALIAS)
            background.paste(img,(ii*new_width,jj*new_height))
            k = k +1
        else:
            break
        
        
background.save('clustering_result.png')