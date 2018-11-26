#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 26 16:41:34 2018

@author: alma
"""

import logging
import datetime

import os.path as osp
import os
import subprocess as sp
import sys

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



def log_header(logger):
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