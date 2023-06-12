'''
Written by Tom Liu, Documentation last updated 2023 June 9

Contains the Python functions that handle all the IO for 
Squid tuning. 
All the sample files are provided in the Github Repo
'''
import glob
import os
import pandas as pd
import csv
import configparser

import mce_data


def get_rsservo_data(dir_path, bias_suffix=None, run_suffix=None):
    '''
    Retrieves the RSSERVO data from the specified directory.

    Inputs:
        dir_path (str): Path to the directory containing the data files.
        bias_suffix (str, optional): Suffix for the .bias file. Defaults to None.
        run_suffix (str, optional): Suffix for the .run file. Defaults to None.

    Returns:
        bias_df (DataFrame): DataFrame containing the .bias data.
        mce_runfile (MCERunfile object): MCERunfile object for the .run file.
    '''
    # The default parameters are written on separate lines in the function
    # rather than in the parameters inputs
    # for easier readability
    if(bias_suffix is None):
        bias_suffix = '_rsservo_sa.bias'
    if(run_suffix is None):
        run_suffix = '_rsservo_sa.run'
    bias_df, mce_runfile = get_bias_run_data(dir_path, bias_suffix, run_suffix)
    return bias_df, mce_runfile


def get_sq1_tune_data(dir_path, bias_suffix=None, run_suffix=None, fast_csv_reading=False):
    '''
    Retrieves the SQ1 tuning data from the specified directory.

    Inputs:
        dir_path (str): Path to the directory containing the data files.
        bias_suffix (str, optional): Suffix for the .bias file. Defaults to None.
        run_suffix (str, optional): Suffix for the .run file. Defaults to None.
        fast_csv_reading (bool, optional): Flag to enable fast CSV reading. Defaults to False.

    Returns:
        bias_df (DataFrame): DataFrame containing the .bias data.
        mce_runfile (MCERunfile object): MCERunfile object for the .run file.
    '''

    # The default parameters are written on separate lines in the function
    # rather than in the parameters inputs
    # for easier readability
    if(bias_suffix is None):
        bias_suffix = '_sq1servo_sa.bias'
    if(run_suffix is None):
        run_suffix = '_sq1servo_sa.run'
    bias_df, mce_runfile = get_bias_run_data(
        dir_path, bias_suffix, run_suffix, fast_csv_reading=fast_csv_reading)
    return bias_df, mce_runfile


def get_bias_run_data(dir_path, bias_suffix='_sq1servo_sa.bias', run_suffix='_sq1servo_sa.run',
                      fast_csv_reading=False):
    '''
    Retrieves the requested .run and .bias data from the specified directory.

    Inputs:
        dir_path (str): Path to the directory containing the data files.
        bias_suffix (str, optional): Suffix for the .bias file. Defaults to '_sq1servo_sa.bias'.
        run_suffix (str, optional): Suffix for the .run file. Defaults to '_sq1servo_sa.run'.
        fast_csv_reading (bool, optional): Flag to enable fast CSV reading. Defaults to False.

    Returns:
        bias_df (DataFrame): DataFrame containing the .bias data.
        mce_runfile (MCERunfile object): MCERunfile object for the .run file.
    '''
    
    search_str = f'{dir_path}/*{bias_suffix}'
    # spaces or commas
    separator = r'\s*,\s*|\s+'
    try:
        bias_file = glob.glob(search_str)[0]
    except IndexError as e:
        print(e)
        raise FileNotFoundError("No files found with pattern: " + search_str)

    bias_path = os.path.join(bias_file)
    print('Reading: ' + bias_path)
    min_column_num = 3
    
    if(fast_csv_reading):
        fast_sep = ','
        with open(bias_path, newline='') as f:
            row1 = ''
            while(row1 == ''):
                row1 = f.readline().strip()
            num_commas = row1.count(',')
        if(num_commas <= min_column_num):
            fast_sep = '\s+'
        bias_df = pd.read_csv(bias_path,  sep=fast_sep,
                              index_col=False, engine='c')

    else:
        try:
            bias_df = pd.read_csv(bias_path,  sep=separator,
                                  on_bad_lines='warn', index_col=False, skipinitialspace=True)
        except TypeError:
            print('Using old version of pandas:')
            bias_df = pd.read_csv(bias_path,  sep=separator,
                                  error_bad_lines=False, index_col=False, skipinitialspace=True)

    print('Columns in .bias file: ' + str(bias_df.columns))
    assert len(
        bias_df.columns) > min_column_num, "Not enough columns in .bias file, data is improperly formatted"
    search_str = f'{dir_path}/*{run_suffix}'
    try:
        run_file = glob.glob(search_str)[0]
    except IndexError as e:
        print(e)
        raise FileNotFoundError("No files found with pattern: " + search_str)

    print('Reading: ' + run_file)
    mce_runfile = mce_data.MCERunfile(run_file)

    return bias_df, mce_runfile


def get_config_file(file_path='tune_cfg/slac_cd19.cfg'):
    '''
    Returns the configuration file containing SSA, SQ1, and other parameters.

    Inputs:
        file_path (str, optional): Path to the configuration file. Defaults to 'tune_cfg/slac_cd19.cfg'.

    Returns:
        cfg (ConfigParser object): Configuration file object.
    '''
    cfg = configparser.ConfigParser()
    print('Reading: ' + file_path)
    cfg.read(file_path)
    return cfg


def get_ssa_tune_data(dir_path, suffix='_ssa', run_ext='.run'):
    '''
    Retrieves the relevant parameters from the SSA tuning data.

    Inputs:
        dir_path (str): Path to the directory containing the SSA tuning data files.
        suffix (str, optional): Suffix used to identify the relevant files. Defaults to '_ssa'.
        run_ext (str, optional): File extension for the run files. Defaults to '.run'.

    Returns:
        sa_data (object): SSA tuning data object.
        sa_runfile (object): SSA tuning run file object.
    '''
    search_str = f'{dir_path}/*{suffix}'
    try:
        sa_tune = glob.glob(search_str)[0]
    except IndexError as e:
        print(e)
        raise FileNotFoundError("No files found with pattern: " + search_str)
    sa_data = os.path.join(sa_tune)
    print('Reading: ' + str(sa_data))
    mcefile = mce_data.MCEFile(sa_data)
    sa_data = mcefile.Read(field='error', row_col=True)
    print('Reading: ' + sa_tune+run_ext)
    sa_runfile = mce_data.MCERunfile(sa_tune+run_ext)

    return sa_data, sa_runfile


def write_optimal_bias_data(cols, biases, savename, savedir):
    '''
    Writes a list of optimal biases for each column into a CSV file.

    Inputs:
        cols (list): List of column values.
        biases (list): List of optimal bias values.
        savename (str): Name of the file to be saved (without file extension).
        savedir (str): Directory where the file should be saved.

    Output:
        None
    '''
    assert(len(cols) == len(biases))
    col_name = 'column'
    bias_name = 'optimal_bias'

    savepath = os.path.join(savedir, savename + '_col_biases.csv')
    d = {col_name: cols,
         bias_name: biases}
    df = pd.DataFrame(d)
    df.to_csv(savepath, index=False)


def read_optimal_bias_data(data_path):
    '''
    Reads a CSV file for the column and optimal biases.
    Returns a dictionary from column to bias.

    Input:
        data_path (str): Path to the CSV file.

    Output:
        col_bias_dict (dict): Dictionary mapping column (int) to optimal bias (float).
    '''
    col_name = 'column'
    bias_name = 'optimal_bias'
    df = pd.read_csv(data_path)

    biases = df[bias_name].to_list()
    cols = df[col_name].to_list()
    col_bias_dict = {int(cols[i]): float(biases[i]) for i in range(len(cols))}
    return col_bias_dict
