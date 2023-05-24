import glob 
import os
import numpy as np
import pandas as pd
import configparser

import mce_data

def read_ivsrc(srcn):
    '''
    copied over from showiv.py
    srcn: path/to/bias_script.scr file
    output: list of all biases
    '''
    bias_marker = 'wb'
    ball  = []
    with open(srcn, 'r') as f:
        array = [[x for x in line.split()] for line in f]
        
        for i in range(len(array)):
            l = array[i]
            if l[0]==bias_marker:
                bcol = np.empty(len(l)-3, dtype=int)
                for j in range(3,len(l)):
                    bcol[j-3] = int(l[j])
                ball.append(bcol)
            else:
                continue
        ball = np.array(ball)
    return ball

def get_rsservo_data(dir_path):
    '''
    input: path/to/mce_file
    output: retrieves rsservo data
    '''
    rsservo = glob.glob(f'{dir_path}/*_rsservo.bias')[0]
    rsservo_path = os.path.join(rsservo)
    print('Reading: ' + rsservo_path)
    try:
        rsservo_df = pd.read_csv(rsservo_path, delim_whitespace=True,
                    on_bad_lines='warn', index_col=False)
    except TypeError:
        rsservo_df = pd.read_csv(rsservo_path, delim_whitespace=True,
                    error_bad_lines=False, index_col=False)
    rsservo_file = rsservo_path.replace('.bias', '.run')
    print('Reading: ' + rsservo_file)
    rsservo_runfile = mce_data.MCERunfile(rsservo_file)
    return rsservo_df, rsservo_runfile  

def get_sq1_tune_data(dir_path):
    '''
    input: path/to/mce_file
    output: retrieves sq1 tuning data
    '''
    sq1_tune = glob.glob(f'{dir_path}/*_sq1servo_sa.bias')[0]
    sq1_data = os.path.join(sq1_tune)
    print('Reading: ' + sq1_data)
    try:
        sq1df = pd.read_csv(sq1_data, delim_whitespace=True,
                    on_bad_lines='warn', index_col=False)
    except TypeError:
        sq1df = pd.read_csv(sq1_data, delim_whitespace=True,
                    error_bad_lines=False, index_col=False)
    sq1_file = sq1_tune.replace('.bias', '.run')
    print('Reading: ' + sq1_file)
    sq1_runfile = mce_data.MCERunfile(sq1_file)
    return sq1df, sq1_runfile  


def get_config_file(file_path='tune_cfg/slac_cd19.cfg'):
    '''
    returns config file containing SSA, SQ1, and other parameters
    '''
    cfg = configparser.ConfigParser()
    print('Reading: ' + file_path)
    cfg.read(file_path)
    return cfg 

def get_ssa_tune_data(dir_path):
    '''
    Retrieves the relevant parameters from the SSA tuning data
    '''
    sa_tune = glob.glob(f'{dir_path}/*_ssa')[0] 
    sa_data = os.path.join(sa_tune)
    print('Reading: ' + str(sa_data))
    mcefile = mce_data.MCEFile(sa_data)
    sa_data = mcefile.Read(field='error', row_col=True)
    print('Reading: ' + sa_tune+'.run')
    sa_runfile = mce_data.MCERunfile(sa_tune+'.run')

    return sa_data, sa_runfile 


