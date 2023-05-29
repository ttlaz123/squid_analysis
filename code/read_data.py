import glob 
import os
import numpy as np
import pandas as pd
import configparser

import mce_data

def get_rsservo_data(dir_path, bias_suffix='_rsservo.bias', run_suffix='_rsservo.run'):
    '''
    input: path/to/mce_folder
            grabs the first file that is found with the given suffixes
    output: retrieves rsservo data
    '''
    bias_suffix = '_rsservo_sa.bias'
    run_suffix = '_rsservo_sa.run'
    return get_bias_run_data(dir_path, bias_suffix, run_suffix)  
    

def get_sq1_tune_data(dir_path):
    '''
    input: path/to/mce_folder
            grabs the first file that is found with the given suffixes
    output: retrieves sq1 tuning data
    '''
    bias_suffix = '_sq1servo_sa.bias'
    run_suffix = '_sq1servo_sa.run'
    return get_bias_run_data(dir_path, bias_suffix, run_suffix)  

def get_bias_run_data(dir_path, bias_suffix='_sq1servo_sa.bias', run_suffix = '_sq1servo_sa.run'):
    '''
    input: path/to/mce_folder
            grabs the first file that is found with the given suffixes
    output: the requested .run and .bias data
    '''
    bias_file = glob.glob(f'{dir_path}/*{bias_suffix}')[0]
    bias_path = os.path.join(bias_file)
    print('Reading: ' + bias_path)
    try:
        bias_df = pd.read_csv(bias_path, delim_whitespace=True,
                    on_bad_lines='warn', index_col=False)
    except TypeError:
        print('Using old version of pandas:')
        bias_df = pd.read_csv(bias_path, delim_whitespace=True,
                    error_bad_lines=False, index_col=False)
        
    run_file = glob.glob(f'{dir_path}/*{run_suffix}')[0]

    print('Reading: ' + run_file)
    mce_runfile = mce_data.MCERunfile(run_file)
    return bias_df, mce_runfile  


def get_config_file(file_path='tune_cfg/slac_cd19.cfg'):
    '''
    returns config file containing SSA, SQ1, and other parameters
    '''
    cfg = configparser.ConfigParser()
    print('Reading: ' + file_path)
    cfg.read(file_path)
    return cfg 

def get_ssa_tune_data(dir_path, suffix='_ssa', run_ext = '.run'):
    '''
    Retrieves the relevant parameters from the SSA tuning data
    '''
    sa_tune = glob.glob(f'{dir_path}/*{suffix}')[0] 
    sa_data = os.path.join(sa_tune)
    print('Reading: ' + str(sa_data))
    mcefile = mce_data.MCEFile(sa_data)
    sa_data = mcefile.Read(field='error', row_col=True)
    print('Reading: ' + sa_tune+run_ext)
    sa_runfile = mce_data.MCERunfile(sa_tune+run_ext)

    return sa_data, sa_runfile 


