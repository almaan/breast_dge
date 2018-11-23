#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 23 17:00:28 2018

@author: alma
"""
import os
import sys
import os.path as osp
import subprocess as sp

def _version_control():
    file_dir = osp.dirname(os.path.abspath(__file__))
    current_dir = os.getcwd()
    
    commands = [f"cd {file_dir:s}",
                f"git describe --tags",
                f"cd {current_dir:s}"]
    
    
    return_values = []
    for cmd in commands:
        proc = sp.Popen(cmd, shell = True, stdout=sp.PIPE)
        return_values.append(proc.communicate()[0])
    
    version = return_values[1].decode('utf-8').replace('\n','')
    
    return(f"git version : {version:s}")


print(_version_control())

command = ' '.join(sys.argv)

print(command)
