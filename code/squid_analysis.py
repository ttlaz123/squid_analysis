'''
Written by Tom Liu, adapting code from David Goldfinger
2023 July 26 last documentation update
'''
import os
import sys

if sys.version_info.major != 3:
    sys.exit("This software is only supported on Python 3")
if sys.version_info.minor < 6:
    sys.exit("This software is expecting Python 3.6+.")
print(f"Running Python version {sys.version_info.major}.{sys.version_info.minor}")

import argparse
import time
import numpy as np

import read_data as rd
import plot_data as pd
import calc_parameters as cp
import calculate_slopes as cs

# we should figure out where these warnings are coming from some day
import warnings
warnings.filterwarnings("ignore")


class NoColumnException(Exception):
    pass


def fill_grid_data(value, row, col, grid=None, max_rows=41, max_cols=32):
    '''
    Fills a grid with a value at a specific location.

    Inputs:
        value: The value to be filled in the grid.
        row: The row index of the location.
        col: The column index of the location.
        grid: Optional. The grid to be filled. If not provided, a grid of size [max_rows, max_cols] will be initialized.
        max_rows: The maximum number of rows in the grid. Defaults to 41.
        max_cols: The maximum number of columns in the grid. Defaults to 32.

    Returns:
        grid: The updated grid with the value filled at the specified location.
    '''
    if(grid is None):
        grid = np.zeros((max_rows, max_cols))
    grid[row, col] = value
    return grid


def get_icmaxcolmod(ic_params_rson, ic_params_rsoff, manual_bias_idx=None):
    '''
    Given ic_params, obtain standard squid parameters
    Args:
        ic_params_rson (dict): Dictionary containing the parameters for the "on" state.
        ic_params_rsoff (dict): Dictionary containing the parameters for the "off" state.
        manual_bias (int, optional): Index of the manual bias. Defaults to None.

    Returns:
        tuple: A tuple containing the following values:
            ic_col (float): The value of ic_col.
            ic_min (float): The value of ic_min.
            ic_max (float): The value of ic_max.
            mod (float): The value of mod.
            optimal_bias (float): The value of optimal_bias.
            crosstalk_bias (float): The value of crosstalk_bias.
            manual_mod_idx (float): The index of chosen bias.
    '''
    sq1_safb_biases = ic_params_rson['bias']
    sq1_safb_min = ic_params_rson['fb_min']
    sq1_safb_max = ic_params_rson['fb_max']
    max_sq1imod_idx = ic_params_rson['bias_max_idx']
    start_sq1imod_idx = ic_params_rson['bias_min_idx']

    ic_min = sq1_safb_min[start_sq1imod_idx]
    ic_max = sq1_safb_max[max_sq1imod_idx]
    mod = sq1_safb_max[max_sq1imod_idx] - sq1_safb_min[max_sq1imod_idx]
    optimal_bias = sq1_safb_biases[max_sq1imod_idx]

    manual_mod = -1
    if(manual_bias_idx is not None):
        manual_mod = (sq1_safb_max[manual_bias_idx] -
                      sq1_safb_min[manual_bias_idx])

    ic_col = -1
    crosstalk_bias = -1
    if(ic_params_rsoff is not None):
        sq1_safb_biases = ic_params_rsoff['bias']
        sq1_safb_min = ic_params_rsoff['fb_min']
        start_sq1imod_idx = ic_params_rsoff['bias_min_idx']

        ic_col = sq1_safb_min[start_sq1imod_idx]
        crosstalk_bias = sq1_safb_biases[start_sq1imod_idx]

    return ic_col, ic_min, ic_max, mod, optimal_bias, crosstalk_bias, manual_mod


def initialize_grid_color_range():
    # Will need adjustment based on the calibration of the particular system
    vmin_vmax_dict = {
        'DAC': {
            'vmin': {
                'Ic,col': 0,
                'Ic,max': 0,
                'Optimal Modulation': 0,
                'Manually chosen Modulation': 0,
                'Optimal Bias': 0,
                'Crosstalk Bias Limit': 0,
                'Optimal Bias - Crosstalk Bias Limit': -5000,
                'Ic,max - Ic,col': -1000,
            },
            'vmax': {
                'Ic,col': 8000,
                'Ic,max': 8000,
                'Optimal Modulation': 2000,
                'Manually chosen Modulation': 2000,
                'Optimal Bias': 20000,
                'Crosstalk Bias Limit': 20000,
                'Optimal Bias - Crosstalk Bias Limit': 5000,
                'Ic,max - Ic,col': 1000,
            }
        },
        'uA': {
            'vmin': {
                'Ic,col': 0,
                'Ic,max': 0,
                'Optimal Modulation': 0,
                'Manually chosen Modulation': 0,
                'Optimal Bias': 0,
                'Crosstalk Bias Limit': 0,
                'Optimal Bias - Crosstalk Bias Limit': -25,
                'Ic,max - Ic,col': -5,
            },
            'vmax': {
                'Ic,col': 40,
                'Ic,max': 40,
                'Optimal Modulation': 10,
                'Manually chosen Modulation': 10,
                'Optimal Bias': 100,
                'Crosstalk Bias Limit': 100,
                'Optimal Bias - Crosstalk Bias Limit': 25,
                'Ic,max - Ic,col': 5,
            }
        },
    }
    return vmin_vmax_dict


def make_grids(all_grids, rows, cols, ctime, show_plot, savedir, convert_units):
    """
    Plot grids based on the provided data.

    Args:
        all_grids (dict): Dictionary containing all grid data.
        rows (list): List of row values.
        cols (list): List of column values.
        ctime: Current time.
        show_plot (bool): Flag to indicate whether to show the plot.
        savedir: Directory to save the plots.
        convert_units (bool): Flag to indicate whether to convert units.

    Returns:
        None
    """
    rows = range(max(rows))
    cols = range(max(cols))
    ic_max_grid = all_grids['ic_max']
    ic_col_grid = all_grids['ic_col']
    ic_maxcol_diff_grid = all_grids['ic_maxcol_diff']
    mod_grid = all_grids['mod']
    optimal_bias_grid = all_grids['opt_bias']
    crosstalk_bias_grid = all_grids['cross_bias']
    bias_crosstalk_diff_grid = all_grids['opt_cross_diff']
    chosen_mod_grid = all_grids['chosen']

    vmin_vmax_dict = initialize_grid_color_range()
    if convert_units:
        uname = 'uA'
        vmin = vmin_vmax_dict['uA']['vmin']
        vmax = vmin_vmax_dict['uA']['vmax']
    else:
        uname = 'DAC'
        vmin = vmin_vmax_dict['DAC']['vmin']
        vmax = vmin_vmax_dict['DAC']['vmax']

    # TODO: there's probably a way to do this in a loop but too lazy to implement
    print('plotting grids...')
    pd.tile_plot(len(rows), len(cols), ic_col_grid,
                 'Ic,col ('+uname+')', str(ctime)+'_Ic_col'+'_units'+uname,
                 show_plot=show_plot, savedir=savedir,
                 vmin=vmin['Ic,col'], vmax=vmax['Ic,col'])

    pd.tile_plot(len(rows), len(cols), ic_max_grid,
                 'Ic,max ('+uname+')', str(ctime)+'_Ic_max'+'_units'+uname,
                 show_plot=show_plot, savedir=savedir,
                 vmin=vmin['Ic,max'], vmax=vmax['Ic,max'])

    pd.tile_plot(len(rows), len(cols), mod_grid,
                 'Optimal Modulation ('+uname+')', str(ctime) +
                 '_optmod'+'_units'+uname,
                 show_plot=show_plot, savedir=savedir,
                 vmin=vmin['Optimal Modulation'], vmax=vmax['Optimal Modulation'])
    pd.tile_plot(len(rows), len(cols), chosen_mod_grid,
                 'Chosen Modulation ('+uname+')', str(ctime) +
                 '_chosenmod'+'_units'+uname,
                 show_plot=show_plot, savedir=savedir,
                 vmin=vmin['Manually chosen Modulation'], vmax=vmax['Manually chosen Modulation'])

    pd.tile_plot(len(rows), len(cols), optimal_bias_grid,
                 'Optimal Bias ('+uname+')', str(ctime) +
                 '_optbias'+'_units'+uname,
                 show_plot=show_plot, savedir=savedir,
                 vmin=vmin['Optimal Bias'], vmax=vmax['Optimal Bias'])
    pd.tile_plot(len(rows), len(cols), crosstalk_bias_grid,
                 'Crosstalk Bias Limit ('+uname+')', str(ctime) +
                 '_crosstalk'+'_units'+uname,
                 show_plot=show_plot, savedir=savedir,
                 vmin=vmin['Crosstalk Bias Limit'], vmax=vmax['Crosstalk Bias Limit'])

    pd.tile_plot(len(rows), len(cols), bias_crosstalk_diff_grid,
                 'Optimal Bias - Crosstalk Bias Limit ('+uname+')', str(
                     ctime)+'_optbias_crosstalk_diff'+'_units'+uname,
                 show_plot=show_plot, savedir=savedir,
                 vmin=vmin['Optimal Bias - Crosstalk Bias Limit'], vmax=vmax['Optimal Bias - Crosstalk Bias Limit'])

    pd.tile_plot(len(rows), len(cols), ic_maxcol_diff_grid,
                 'Ic,max - Ic,col ('+uname+')', str(ctime) +
                 '_Ic_maxcol_diff'+'_units'+uname,
                 show_plot=show_plot, savedir=savedir,
                 vmin=vmin['Ic,max - Ic,col'], vmax=vmax['Ic,max - Ic,col'])
    return


def find_bias_idx(sq1_runfile, chosen_bias):
    '''
    Finds the index corresponding to the chosen bias value in the bias DACs.

    Parameters:
        - sq1_runfile (MCEFile): SQ1 runfile.
        - chosen_bias (float): The chosen optimal bias value.

    Returns:
        - idx (int): The index corresponding to the chosen bias value.

    Note:
        - The function assumes that the bias DACs and feedback DACs ranges are obtained
          using the `cp.get_bias_fb_range` function from the same runfile.
    '''
    bias_dacs, fb_dacs = cp.get_bias_fb_range(sq1_runfile)
    if(chosen_bias > max(bias_dacs) or chosen_bias < min(bias_dacs)):
        print("Warning: chosen optimal bias out of range: " + str(chosen_bias))
        print("Bias Dacs: " + str(bias_dacs))
    idx = (chosen_bias <= bias_dacs).argmax()
    return idx


def fill_all_ic_grids(all_grids, col, ic_params_rson_allrows,
                      ic_params_rsoff_allrows, chosen_bias_idx,
                      max_rows=41, max_cols=32):
    """
    Fills the IC grids with data based on the provided IC parameters and chosen bias.

    Args:
        all_grids (dict): A dictionary containing the IC grids.
        col (int): The column index.
        ic_params_rson_allrows (dict): A dictionary containing IC parameters for each row in the 'rson' dataset.
        ic_params_rsoff_allrows (dict): A dictionary containing IC parameters for each row in the 'rsoff' dataset.
        chosen_bias_idx (int): The chosen bias index.

    Returns:
        dict: The updated dictionary of IC grids.
    """
    ic_max_grid = all_grids['ic_max']
    ic_col_grid = all_grids['ic_col']
    ic_maxcol_diff_grid = all_grids['ic_maxcol_diff']
    mod_grid = all_grids['mod']
    optimal_bias_grid = all_grids['opt_bias']
    crosstalk_bias_grid = all_grids['cross_bias']
    bias_crosstalk_diff_grid = all_grids['opt_cross_diff']
    chosen_mod_grid = all_grids['chosen']

    for row in ic_params_rson_allrows:
        ic_params_rson = ic_params_rson_allrows[row]
        if(ic_params_rsoff_allrows is not None):
            ic_params_rsoff = ic_params_rsoff_allrows[row]
        else:
            ic_params_rsoff = None
        (ic_col, ic_min, ic_max, mod,
         optimal_bias, crosstalk_bias, manual_mod) = get_icmaxcolmod(
            ic_params_rson, ic_params_rsoff, manual_bias_idx=chosen_bias_idx)

        ic_col_grid = fill_grid_data(
            ic_col, row, col, grid=ic_col_grid,
            max_rows=max_rows, max_cols=max_cols)
        ic_max_grid = fill_grid_data(
            ic_max, row, col, grid=ic_max_grid,
            max_rows=max_rows, max_cols=max_cols)
        ic_maxcol_diff_grid = fill_grid_data(
            ic_max-ic_col, row, col, grid=ic_maxcol_diff_grid,
            max_rows=max_rows, max_cols=max_cols)

        mod_grid = fill_grid_data(
            mod, row, col, grid=mod_grid,
            max_rows=max_rows, max_cols=max_cols)
        optimal_bias_grid = fill_grid_data(
            optimal_bias, row, col, grid=optimal_bias_grid,
            max_rows=max_rows, max_cols=max_cols)
        crosstalk_bias_grid = fill_grid_data(
            crosstalk_bias, row, col, grid=crosstalk_bias_grid,
            max_rows=max_rows, max_cols=max_cols)
        bias_crosstalk_diff_grid = fill_grid_data(
            optimal_bias-crosstalk_bias, row, col, grid=bias_crosstalk_diff_grid,
            max_rows=max_rows, max_cols=max_cols)

        chosen_mod_grid = fill_grid_data(
            manual_mod, row, col, grid=chosen_mod_grid,
            max_rows=max_rows, max_cols=max_cols)

    all_grids['ic_max'] = ic_max_grid
    all_grids['ic_col'] = ic_col_grid
    all_grids['mod'] = ic_maxcol_diff_grid
    all_grids['ic_maxcol_diff'] = mod_grid
    all_grids['opt_bias'] = optimal_bias_grid
    all_grids['cross_bias'] = crosstalk_bias_grid
    all_grids['opt_cross_diff'] = bias_crosstalk_diff_grid
    all_grids['chosen'] = chosen_mod_grid

    return all_grids


def get_icparams_squid_column(col, sa_data, sa_runfile, cfg,
                              sq1df_rson, sq1_runfile_rson,
                              sq1df_rsoff=None, sq1_runfile_rsoff=None,
                              convert_units=False, flip_signs=False, mod_thresh=20,
                              verbose=False):
    '''
    Retrieves the IC parameters for a specific column in an SQ1 array.

    Parameters:
        col (int): Column number to analyze.
        sa_data (MCEFile): SA data file.
        sa_runfile (MCEFile): SA run file.
        cfg (CFGFile): Cconfiguration file.
        sq1df_rson (pandas DataFrame): Dataframe containing SQ1 data with row select on.
        sq1_runfile_rson (MCEFile): SQ1 run file with row select on.
        sq1df_rsoff (pandas DataFrame, optional): Dataframe containing SQ1 data with row select off. Defaults to None.
        sq1_runfile_rsoff (MCEFile, optional): SQ1 run file with row select off. Defaults to None.
        convert_units (bool, optional): Flag to convert units of IC parameters. Defaults to False.
        flip_signs (bool, optional): Flag to flip signs of IC parameters. Defaults to False.
        verbose (bool, optional): Flag to enable verbose output. Defaults to False.

    Returns:
        ic_params_rson_allrows (dict): Dictionary containing IC parameters for each row with row select on.
        ic_params_rsoff_allrows (dict): Dictionary containing IC parameters for each row with row select off.

    Raises:
        NoColumnException: If the specified column is missing or cannot be processed.
    '''
    bname = '<bias>'
    fluxname = '<flux>'
    rowname = '<row>'
    colname = '<safb' + str(str(col).zfill(2)) + '>'

    print("Analyzing Column: " + str(col))
    try:
        ssa_params = cp.calculate_ssa_parameters(
            sa_data, sa_runfile, cfg, col)
    except TypeError as e:
        print('Skipping Column: ' + str(col))
        raise NoColumnException(e)
    if(ssa_params is None):
        print('Skipping Column: ' + str(col))
        raise NoColumnException("Cannot process SA for column " + str(col))

    sq1df_col = sq1df_rson.filter([bname, fluxname, rowname, colname], axis=1)
    if(sq1df_col.shape[1] < 4):
        print('Skipping Column: ' + str(col))
        raise NoColumnException("No column data for column " + str(col))

    if(sq1df_rsoff is not None):
        sq1df_off_col = sq1df_rsoff.filter(
            [bname, fluxname, rowname, colname], axis=1)

    rows = np.unique(sq1df_rson[rowname].tolist())
    ic_params_rson_allrows = {}
    ic_params_rsoff_allrows = {}
    if(sq1df_rsoff is None):
        ic_params_rsoff_allrows = None
    for row in rows:
        sq1df_row = sq1df_col[sq1df_col[rowname] == row]
        if(sq1df_row.shape[0] < 5):
            if(verbose):
                print("Not enough data for row " +
                      str(row), '. Skipping row.')
            continue
        ic_params = cp.calculate_ic_params(sq1df_row, sq1_runfile_rson,
                                           col, mod_thresh=mod_thresh,
                                           convert_units=convert_units,
                                           cfg=cfg, ssa_params=ssa_params,
                                           flip_signs=flip_signs, verbose=verbose)
        ic_params_rson_allrows[row] = ic_params
        if(sq1df_rsoff is not None):
            sq1df_off_row = sq1df_off_col[sq1df_off_col[rowname] == row]
            if(sq1df_row.shape[0] < 5):
                if(verbose):
                    print("Not enough data for row " +
                          str(row), '. Skipping row.')
                continue
            ic_params_rsoff = cp.calculate_ic_params(sq1df_off_row, sq1_runfile_rsoff,
                                                     col, mod_thresh=mod_thresh,
                                                     convert_units=convert_units,
                                                     cfg=cfg, ssa_params=ssa_params,
                                                     flip_signs=flip_signs, verbose=verbose)
            ic_params_rsoff_allrows[row] = ic_params_rsoff

    return ic_params_rson_allrows, ic_params_rsoff_allrows


def choose_bias(ic_params_rson_allrows, ic_params_rsoff_allrows=None,
                method='naive'):
    '''
    Determines the best bias current based on the given IC parameters.

    Inputs:
        ic_params_rson_allrows: A dictionary containing IC parameters with row select on.
        ic_params_rsoff_allrows: Optional. A dictionary containing IC parameters with row select off.
                                 Required for methods 'device_current' and 'bias_current'.
        method: The method to choose the bias current. Options: 'naive', 'device_current', 'bias_current'.
                Defaults to 'naive'.

    Returns:
        chosen_bias_idx: The index of the chosen bias current.
        chosen_bias: The chosen bias current.

    Raises:
        ValueError: If the chosen method is not one of the available methods.

    Note:
        The IC parameters dictionaries should have the following keys:
        - 'bias_max_idx': The index of the maximum modulation bias.
        - 'bias_min_idx': The index of the minimum modulation bias.
        - 'fb_max': A list of feedback currents for different biases.
        - 'fb_min': A list of feedback currents for different biases.
        - 'bias': A list of bias currents.
    '''
    methods = ['naive', 'device_current', 'bias_current']
    if(method == methods[0]):
        best_bias_per_row = []
        for row in ic_params_rson_allrows:
            max_mod_bias_idx = ic_params_rson_allrows[row]['bias_max_idx']
            best_bias_per_row.append(max_mod_bias_idx)
        chosen_bias_idx = int(np.mean(best_bias_per_row))
        chosen_bias = ic_params_rson_allrows[row]['bias'][chosen_bias_idx]
        return chosen_bias_idx, chosen_bias

    elif(method == methods[1]):
        assert ic_params_rsoff_allrows is not None, \
            ("Row Select Off data is required to choose bias with method:" + str(method))
        best_bias_per_row = []
        for row in ic_params_rson_allrows:
            max_mod_bias_idx = ic_params_rson_allrows[row]['bias_max_idx']

            ic_max = ic_params_rson_allrows[row]['fb_max'][max_mod_bias_idx]

            crosstalk_bias_limit_idx = ic_params_rsoff_allrows[row]['bias_min_idx']
            ic_col = ic_params_rsoff_allrows[row]['fb_min'][crosstalk_bias_limit_idx]

            # if ic_max is too large, iteratively decrease bias current
            while(ic_max > ic_col):
                max_mod_bias_idx -= 1
                if(max_mod_bias_idx == 0):
                    break
                ic_max = ic_params_rson_allrows[row]['fb_max'][max_mod_bias_idx]

            best_bias_per_row.append(max_mod_bias_idx)

        chosen_bias_idx = min(best_bias_per_row)
        chosen_bias = ic_params_rson_allrows[row]['bias'][chosen_bias_idx]
        return chosen_bias_idx, chosen_bias

    elif(method == methods[2]):
        assert ic_params_rsoff_allrows is not None, \
            ("Row Select Off data is required to choose bias with method:" + str(method))
        best_bias_per_row = []
        all_bias_limits = []
        for row in ic_params_rson_allrows:
            max_mod_bias_idx = ic_params_rson_allrows[row]['bias_max_idx']
            best_bias_per_row.append(max_mod_bias_idx)
            crosstalk_bias_limit_idx = ic_params_rsoff_allrows[row]['bias_min_idx']
            all_bias_limits.append(crosstalk_bias_limit_idx)

        naive_chosen_bias = np.mean(best_bias_per_row)
        bias_limit = min(all_bias_limits)
        if(naive_chosen_bias < bias_limit):
            chosen_bias_idx = int(naive_chosen_bias)
        else:
            chosen_bias_idx = int(bias_limit - 1)
        chosen_bias = ic_params_rson_allrows[row]['bias'][chosen_bias_idx]
        return chosen_bias_idx, chosen_bias
    else:
        raise ValueError('Chosen method ' + str(method) + ' is not ' +
                         'one of the availabel methods:' + str(methods))


def setup_grids():
    """
    Initializes and returns a dictionary containing various grids for IC parameters.

    Returns:
        dict: A dictionary containing the initialized IC parameter grids.
    """
    all_grids = {}
    ic_max_grid = None
    ic_col_grid = None
    ic_maxcol_diff_grid = None
    mod_grid = None
    optimal_bias_grid = None
    crosstalk_bias_grid = None
    bias_crosstalk_diff_grid = None
    chosen_mod_grid = None
    all_grids['ic_max'] = ic_max_grid
    all_grids['ic_col'] = ic_col_grid
    all_grids['mod'] = ic_maxcol_diff_grid
    all_grids['ic_maxcol_diff'] = mod_grid
    all_grids['opt_bias'] = optimal_bias_grid
    all_grids['cross_bias'] = crosstalk_bias_grid
    all_grids['opt_cross_diff'] = bias_crosstalk_diff_grid
    all_grids['chosen'] = chosen_mod_grid
    return all_grids


def setup_directories(savedir):
    """
    Creates and sets up output directories for saving different types of files.

    Parameters:
        savedir (str): The base directory where subdirectories will be created.

    Returns:
        tuple: A tuple containing three directory paths corresponding to the created subdirectories:
            - savedir_cols (str): Directory path for 'col_summary' subdirectory.
            - savedir_rows (str): Directory path for 'all_rows' subdirectory.
            - savedir_grids (str): Directory path for 'gridplots' subdirectory.

    Note:
        This function creates three subdirectories inside the specified 'savedir'
        if they do not already exist. If the subdirectories already exist, it will not modify them.
    """
    col_summary_name = 'col_summary'
    all_rows_name = 'all_rows'
    grid_name = 'gridplots'
    savedir_cols = os.path.join(savedir, col_summary_name)
    while not os.path.exists(savedir_cols):
        os.makedirs(savedir_cols)
    savedir_rows = os.path.join(savedir, all_rows_name)
    while not os.path.exists(savedir_rows):
        os.makedirs(savedir_rows)
    savedir_grids = os.path.join(savedir, grid_name)
    while not os.path.exists(savedir_grids):
        os.makedirs(savedir_grids)
    return savedir_cols, savedir_rows, savedir_grids

def ic_driver(sq1df_rson, sq1_runfile_rson, ctime=None,
              sq1df_off=None,  sq1_runfile_off=None,
              cols=range(0, 16), flip_signs=False, manually_chosen_biases=None,
              plot_all_rows=False, savedir='test_output',
              # For unit conversion
              convert_units=False, cfg=None, sa_data=None, sa_runfile=None,
              # Debug options
              verbose=False, show_plot=False):
    """
    Perform analysis on the provided data and generate plots.

    Parameters:
        sq1df_rson (DataFrame): Dataframe containing resonant data.
        sq1_runfile_rson (str): Path to the runfile used for resonant data.
        ctime (str or None): Custom timestamp for plot filenames.
        sq1df_off (DataFrame or None): Dataframe containing off-resonance data.
        sq1_runfile_off (str or None): Path to the runfile used for off-resonance data.
        cols (range or list): Columns to analyze (default is range from 0 to 15).
        flip_signs (bool): Whether to flip the signs of certain parameters (default is False).
        manually_chosen_biases (dict or None): Manually chosen biases for each column (optional).
        plot_all_rows (bool): Whether to plot individual rows (default is False).
        savedir (str): Base directory to save output files (default is 'test_output').
        convert_units (bool): Whether to convert units (default is False).
        cfg (object or None): Configuration object (optional).
        sa_data (object or None): SA data object (optional).
        sa_runfile (object or None): SA runfile object (optional).
        verbose (bool): Whether to print detailed information (default is False).
        show_plot (bool): Whether to display plots interactively (default is False).

    Returns:
        None
    """

    # Some Options that can be changed

    # bias_current, device_current, or naive
    bias_choose_method = 'naive'
    mod_thresh = 30
    max_rows = None

    
    savedir_cols, savedir_rows, savedir_grids = setup_directories(savedir)
    # Start Analysis
    all_grids = setup_grids()
    fig = None
    ax = None
    optimal_col_biases = np.zeros(len(cols))
    count = 0
    rows = None
    max_cols = max(cols)+1
    for col in cols:
        count += 1
        try:
            (ic_params_rson_allrows,
             ic_params_rsoff_allrows) = get_icparams_squid_column(
                col, sa_data, sa_runfile, cfg,
                sq1df_rson, sq1_runfile_rson,
                sq1df_rsoff=sq1df_off, sq1_runfile_rsoff=sq1_runfile_off,
                convert_units=convert_units, flip_signs=flip_signs, mod_thresh=mod_thresh,
                verbose=verbose)
        except NoColumnException as e:
            print("Skipping column " + str(col))
            optimal_col_biases[count-1] = -1
            continue

        rows = ic_params_rson_allrows.keys()
        if(max_rows is None):
            max_rows = max(rows)+1
        if(manually_chosen_biases is None):
            chosen_bias_idx, chosen_bias = choose_bias(
                ic_params_rson_allrows, ic_params_rsoff_allrows,
                method=bias_choose_method)

        else:
            chosen_bias = manually_chosen_biases[col]
            chosen_bias_idx = find_bias_idx(sq1_runfile_rson, chosen_bias)
        optimal_col_biases[count-1] = chosen_bias

        all_grids = fill_all_ic_grids(all_grids, col, ic_params_rson_allrows,
                                      ic_params_rsoff_allrows, chosen_bias_idx,
                                      max_rows=max_rows, max_cols=max_cols)
        if(plot_all_rows):
            for row in ic_params_rson_allrows:
                time1 = time.time()
                ic_params_rson = ic_params_rson_allrows[row]
                if(ic_params_rsoff_allrows is not None):
                    ic_params_rsoff = ic_params_rsoff_allrows[row]
                else:
                    ic_params_rsoff = None
                print("reading: " + str(time1-time.time()))
                fig, ax = pd.plot_icminmax(col, row, ic_params_rson,
                                           ic_params_rsoff=ic_params_rsoff,
                                           ctime=ctime, convert_units=convert_units,
                                           chosen_bias_idx=chosen_bias_idx,
                                           savedir=savedir_rows,
                                           fig=fig, ax=ax,
                                           show_plot=show_plot, verbose=verbose)
                print("total: " + str(time1-time.time()))
        print('Plotting Summary for col ' + str(col))
        fig, ax = pd.plot_icminmax_column(col, ic_params_rson_allrows,
                                          ic_params_rsoff_allrows=ic_params_rsoff_allrows,
                                          ctime=ctime, convert_units=convert_units,
                                          chosen_bias_idx=chosen_bias_idx,
                                          savedir=savedir_cols,
                                          fig=fig, ax=ax,
                                          show_plot=show_plot, verbose=verbose)

    rd.write_optimal_bias_data(cols, optimal_col_biases,
                               ctime, savedir)

    make_grids(all_grids, rows, cols, ctime,
               show_plot, savedir_grids, convert_units)


def rs_driver(cfg, sa_data, sa_runfile, rsdf, rs_runfile, ctime=None,
              rsdf_off=None,  rs_runfile_off=None, filter_sq1=True, flip_signs=False,
              cols=range(0, 16), rows=range(0, 40)):
    """
    TODO: Needs maintenance
    """
    chip_starts = [0, 11, 21, 31, 41]
    sq1_sgfilter_window_length = 1
    sq1_sgfilter_poly_deg = 2
    calc_slopes = False
    show_plot = False
    bname = '<bias>'
    fluxname = '<flux>'
    rowname = '<row>'
    
    for col in cols:
        colname = '<safb' + str(str(col).zfill(2)) + '>'
        sq1df_col = rsdf.filter([bname, fluxname, rowname, colname], axis=1)
        
        try:
            print("Analyzing Column: " + str(col))
            ssa_params = cp.calculate_ssa_parameters(
                sa_data, sa_runfile, cfg, col)
        except TypeError:
            print('Skipping Column: ' + str(col))
            continue
        if(ssa_params is None):
            print('No SSA params for Column: ' + str(col))
            
        fig = None
        ax = None

        for row in rows:
            chip_num = -1
            if(row >= chip_starts[0] and row < chip_starts[1]):
                chip_num = 0
            elif(row >= chip_starts[1] and row < chip_starts[2]):
                chip_num = 1
            elif(row >= chip_starts[2] and row < chip_starts[3]):
                chip_num = 2
            elif(row >= chip_starts[3] and row < chip_starts[4]):
                chip_num = 3
            else:
                raise ValueError(
                    "Row does not correspond to chip: " + str(row))
            last_fig = (row == rows[-1])
            sq1df_row = sq1df_col[sq1df_col[rowname] == row]
            ic_params = cp.calculate_ic_params(sq1df_row, rs_runfile,col, flip_signs=flip_signs,
                                               filter_sq1=filter_sq1, sq1_sgfilter_window_length=sq1_sgfilter_window_length)

            sq1_params2 = None
            if(rsdf_off is not None):
                sq1_params2 = cp.calculate_ic_params(rsdf_off, rs_runfile_off, cfg, col, row, 
                                                          ssa_params,  flip_signs=flip_signs,filter_sq1=filter_sq1, sq1_sgfilter_window_length=sq1_sgfilter_window_length)

           
            fig, ax = pd.plot_rsservo_col(last_fig, col, chip_num, ic_params,
                                                                sq1_params2=sq1_params2, ctime=ctime,
                                                                ax=ax,
                                                                fig=fig)


def save_subset(df, savename, rows=[28, 29, 30, 31, 32, 33], cols=[8, 9, 10, 11, 12]):
    '''
    saves a small subset of the dataframe
    DEPRECATED, DOES WEIRD THINGS DO NOT USE
    '''
    bname = '<bias>'
    fluxname = '<flux>'
    rowname = '<row>'

    all_cols = [bname, fluxname, rowname]
    for col in cols:
        colname = '<safb' + str(str(col).zfill(2)) + '>'
        all_cols.append(colname)
    cols_df = df.filter(all_cols, axis=1)
    rows_df = cols_df[(cols_df['<row>'].isin(rows))]

    rows_df.to_csv(savename)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--ctime', help='/path/to/ctime/directory')
    parser.add_argument(
        '-u', '--ctime_off', help='/path/to/ctime/directory when the row selects are turned off')
    parser.add_argument('-c', '--config', help='/path/to/config/file')
    parser.add_argument('-i', '--dev_cur', action='store_true',
                        help='whether to perform device current analysis')
    parser.add_argument('-r', '--rsservo', action='store_true',
                        help='whether to perform rsservo analysis')
    parser.add_argument('-s', '--flip_signs', action='store_true',
                        help='whether to flip signs for safb')
    parser.add_argument('-f', '--fast_csv_reading', action='store_true',
                        help='Read csvs faster')
    parser.add_argument('-p', '--chosen_biases', default=None,
                        help='path/to/csv/file with chosen biases')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Whether to print out debug statements')
    parser.add_argument('-l', '--columns', type=int, default=16,
                        help='Number of columns to plot data for')
    args = parser.parse_args()
    flip_signs = args.flip_signs
    numcols = int(args.columns)
    cols = range(0, numcols)
    fast_csv_reading = args.fast_csv_reading
    print('Reading in files:' + str(args.ctime))
    ctime = os.path.basename(os.path.dirname(args.ctime))
    print('ctime: ' + str(ctime))
    cfg = rd.get_config_file()

    if(args.dev_cur):
        time0 = time.time()
        sa_data, sa_runfile = rd.get_ssa_tune_data(args.ctime)
        sq1df, sq1_runfile = rd.get_sq1_tune_data(
            args.ctime, fast_csv_reading=fast_csv_reading)
        all_rows = sq1df['<row>'].astype(int)

        rows = np.unique(all_rows)
        sq1df_off = None
        sq1_runfile_off = None
        plot_all_rows = False
        if(args.ctime_off is not None):
            sq1df_off, sq1_runfile_off = rd.get_sq1_tune_data(
                args.ctime_off, fast_csv_reading=fast_csv_reading)
            # save_subset(sq1df_off, 'rowsel_off_small_sq1servo_sa.bias')
        time1 = time.time()

        output_dir = './output_data/'
        print('Done reading files, time elapsed (s):' + str(time1-time0))
        for convert_units in [False, True, ]:
            if(args.chosen_biases is not None):
                manual_optbias_filepath = args.chosen_biases
                col_bias_dict = rd.read_optimal_bias_data(
                    manual_optbias_filepath)
                manually_chosen_biases = col_bias_dict
            else:
                manually_chosen_biases = None

            if(convert_units):
                branch = 'units_ua'
            else:
                branch = 'units_dac'
            savedir = os.path.join(output_dir, ctime, branch)
            while not os.path.exists(savedir):
                os.makedirs(savedir)

            ic_driver(sq1df, sq1_runfile, ctime=ctime,
                      sq1df_off=sq1df_off,  sq1_runfile_off=sq1_runfile_off,
                      manually_chosen_biases=manually_chosen_biases, cols=cols,
                      savedir=savedir,  plot_all_rows=plot_all_rows, flip_signs=flip_signs,
                      convert_units=convert_units, cfg=cfg, sa_data=sa_data, sa_runfile=sa_runfile,
                      verbose=args.verbose)

    if(args.rsservo):
        time0 = time.time()
        sa_data, sa_runfile = rd.get_ssa_tune_data(args.ctime)
        rsservo_df, rsservo_runfile = rd.get_rsservo_data(args.ctime)
        numrows = max(rsservo_df['<row>'].astype(int))+1
        rows = range(0, numrows)
        rsservo_df_off = None
        rsservo_runfile_off = None
        if(args.ctime_off is not None):
            rsservo_df_off, rsservo_runfile_off = rd.get_rsservo_data(
                args.ctime_off)

        time1 = time.time()
        print('Done reading files, time elapsed (s):' + str(time1-time0))
        ctime = ctime+'_rsservo'
        rs_driver(cfg, sa_data, sa_runfile, rsservo_df, rsservo_runfile,
                  rsdf_off=rsservo_df_off,  rs_runfile_off=rsservo_runfile_off,
                  filter_sq1=False, ctime=ctime, flip_signs=flip_signs,
                  cols=cols, rows=rows)


if __name__ == '__main__':
    time0 = time.time()
    main()
    time1 = time.time()
    print("Analysis complete, time elapsed (s): " + str(time1-time0))
