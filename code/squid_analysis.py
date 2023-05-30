'''
Written by Tom Liu, adapting code from David Goldfinger
2023 May 
'''
import os
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



def fill_grid_data(value, row, col, grid=None, max_rows=41, max_cols=32):
    '''
    Fills in grid with a value
    '''
    if(grid is None):
        grid = np.zeros((max_rows, max_cols))
    grid[row, col] = value
    return grid


def get_icmaxcolmod(ic_params, ic_params2,manual_bias=None):
    '''
    Gets ic col, ic max, and modulation given the ic_params
    '''
    sq1_safb_servo_biases_uA = ic_params['bias']
    sq1_safb_servo_mins_sa_in_uA = ic_params['fb_min']
    sq1_safb_servo_maxs_sa_in_uA = ic_params['fb_max'] 
    max_sq1imod_idx = ic_params['bias_max_idx']  
    start_sq1imod_idx = ic_params['bias_min_idx']  

    ic_min = sq1_safb_servo_mins_sa_in_uA[start_sq1imod_idx]
    ic_max = sq1_safb_servo_maxs_sa_in_uA[max_sq1imod_idx]
    mod = -(sq1_safb_servo_mins_sa_in_uA[max_sq1imod_idx] -
            sq1_safb_servo_maxs_sa_in_uA[max_sq1imod_idx])
    manual_mod = None
    if(manual_bias is not None):
        manual_mod = -(sq1_safb_servo_mins_sa_in_uA[manual_bias] -
            sq1_safb_servo_maxs_sa_in_uA[manual_bias])
    optimal_bias = sq1_safb_servo_biases_uA[max_sq1imod_idx]

    ic_col = -1
    if(ic_params2 is not None):
        sq1_safb_servo_biases_uA = ic_params2['bias']
        sq1_safb_servo_mins_sa_in_uA = ic_params2['fb_min']
        sq1_safb_servo_maxs_sa_in_uA = ic_params2['fb_max'] 
        max_sq1imod_idx = ic_params2['bias_max_idx']  
        start_sq1imod_idx = ic_params2['bias_min_idx']  

        ic_col = sq1_safb_servo_maxs_sa_in_uA[start_sq1imod_idx]
        crosstalk_bias = sq1_safb_servo_biases_uA[start_sq1imod_idx]

    return ic_col, ic_min, ic_max, mod, optimal_bias, crosstalk_bias, manual_mod


def make_grids(rows, cols, ctime, show_plot, savedir, convert_units,
               ic_col_grid, ic_max_grid, mod_grid, optimal_bias_grid, crosstalk_bias_grid,
               bias_crosstalk_diff_grid, ic_maxcoldiff_grid, manual_mod):
    rows = range(max(rows))
    cols = range(max(cols))
    if(convert_units):
        uname = 'uA'
        vmin = 5
        vmax = 15
    else:
        uname = 'DAC'
        vmin = 2000
        vmax = 6000


    print('plotting grids...')
    pd.tile_plot(len(rows), len(cols), ic_col_grid,
                 'Ic,col ('+uname+')', str(ctime)+'_Ic_col'+'_units'+uname, 
                 show_plot=show_plot, savedir=savedir, vmin =vmin, vmax=vmax)
    
    pd.tile_plot(len(rows), len(cols), ic_max_grid,
                 'Ic,max ('+uname+')', str(ctime)+'_Ic_max'+'_units'+uname, 
                 show_plot=show_plot, savedir=savedir, vmin =vmin, vmax=vmax)
    
    if(convert_units):
        vmin = 0
        vmax = 5
    else:
        vmin = 0
        vmax = 2000
    pd.tile_plot(len(rows), len(cols), mod_grid,
                 'Optimal Modulation ('+uname+')', str(ctime)+'_optmod'+'_units'+uname,
                   show_plot=show_plot, savedir=savedir, vmin =vmin, vmax=vmax)
    pd.tile_plot(len(rows), len(cols), manual_mod,
                 'Manullay Picked Modulation ('+uname+')', str(ctime)+'_manualmod'+'_units'+uname,
                   show_plot=show_plot, savedir=savedir, vmin =vmin, vmax=vmax)
    
    if(convert_units):
        vmin = 1000
        vmax = 3000
    else:
        vmin = 5000
        vmax = 15000
    pd.tile_plot(len(rows), len(cols), optimal_bias_grid,
                 'Optimal Bias ('+uname+')', str(ctime)+'_optbias'+'_units'+uname, 
                 show_plot=show_plot, savedir=savedir, vmin =vmin, vmax=vmax)
    pd.tile_plot(len(rows), len(cols), crosstalk_bias_grid,
                 'Crosstalk Bias Limit ('+uname+')', str(ctime)+'_crosstalk'+'_units'+uname, 
                 show_plot=show_plot, savedir=savedir, vmin =vmin, vmax=vmax)
    
    if(convert_units):
        vmin = -1000
        vmax = 1000
    else:
        vmin = -5000
        vmax = 5000

    pd.tile_plot(len(rows), len(cols), bias_crosstalk_diff_grid,
                 'Optimal Bias - Crosstalk Bias Limit ('+uname+')', str(ctime)+'_optbias_crosstalk_diff'+'_units'+uname, 
                 show_plot=show_plot, savedir=savedir, vmin =vmin, vmax=vmax)
    
    if(convert_units):
        vmin = -5
        vmax = 5
    else:
        vmin = -2000
        vmax = 2000

    pd.tile_plot(len(rows), len(cols), ic_maxcoldiff_grid,
                 'Ic,max - Ic,col ('+uname+')', str(ctime)+'_Ic_maxcol_diff'+'_units'+uname, 
                 show_plot=show_plot, savedir=savedir, vmin =vmin, vmax=vmax)
    return 

def ic_driver(sq1df, sq1_runfile, ctime=None,
              sq1df_off=None,  sq1_runfile_off=None,
              cols=range(0, 16), rows=range(0, 40),
              plot_all_rows=False, savedir = 'output_data',
              convert_units=False, cfg=None, sa_data=None, sa_runfile=None,
              verbose=False):
    ## TODO: make it automatically pick if there's no provided manually picked file
    ## TODO: make script auto generate pager
    manually_picked_biases = np.array([
        6000,7500,6000, 8000, 
        5000,7000, 7500, 6500,
        9000, 7500, 8000, 8000,
        8000, 8000, 8000, 8000,
        9000, 9000, 9000, 9000,
        9000, 9000, 9000, 9000,
        9000, 9000, 9000, 9000,
        9000, 9000, 9000, 9000,
    ])
 
    show_plot = False
    
    bname = '<bias>'
    fluxname = '<flux>'
    rowname = '<row>'
    
    
    ic_max_grid = None
    ic_col_grid = None
    ic_maxcoldiff_grid = None
    mod_grid = None
    optimal_bias_grid = None
    crosstalk_bias_grid = None
    bias_crosstalk_diff_grid = None
    manual_mod_grid = None
    fig = None
    ax = None
    s1b_minmax_fig = None
    s1b_minmax_ax = None
    savedir_cols = os.path.join(savedir, 'col_summary')
    while not os.path.exists(savedir_cols):
                os.makedirs(savedir_cols)
    savedir_rows = os.path.join(savedir, 'all_rows')
    while not os.path.exists(savedir_rows):
                os.makedirs(savedir_rows)          
    for col in cols:
        colname = '<safb' + str(str(col).zfill(2)) + '>'
        
        print("Analyzing Column: " + str(col))
        try:
            ssa_params = cp.calculate_ssa_parameters(
                sa_data, sa_runfile, cfg, col)
        except TypeError as e:
            print('Skipping Column: ' + str(col))
            print("Error: " + str(e))
            continue
        if(ssa_params is None):
            print('Skipping Column: ' + str(col))
            continue
       
        sq1_b0, d_sq1_b, n_sq1_b = tuple([
        int(i) for i in sq1_runfile.Item('par_ramp', 'par_step loop1 par1')])
        sq1_b = sq1_b0 + d_sq1_b*np.arange(n_sq1_b)
        manual_bias = manually_picked_biases[col]
        manual_bias_idx = (manual_bias <= sq1_b).argmax()
        
        sq1df_col = sq1df.filter([bname,fluxname,rowname, colname], axis=1)
        if(sq1df_col.shape[1] < 4):
            
            print("No column data, skpping column " + str(col))
            continue
        if(sq1df_off is not None):
            sq1df_off_col = sq1df_off.filter([bname,fluxname,rowname, colname], axis=1)
        for row in rows:
            sq1df_row = sq1df_col[sq1df_col[rowname] == row]
            if(sq1df_row.shape[0] < 5):
                if(verbose):
                    print("Not enough data for row " + str(row), '. Skipping row.')
                continue
            
            last_fig = (row == rows[-1])
            ic_params = cp.calculate_ic_params(sq1df_row, sq1_runfile, col, mod_thresh = 20,
                        convert_units = False, cfg=None, ssa_params=None)
            if(sq1df_off is not None):
                sq1df_off_row = sq1df_off_col[sq1df_off_col[rowname] ==row]
           
                ic_params2 = cp.calculate_ic_params(sq1df_off_row, sq1_runfile_off, col, mod_thresh = 20,
                        convert_units = False, cfg=None, ssa_params=None)
            else:
                ic_params2 = None
            
            if(plot_all_rows):
                fig, ax= pd.plot_icminmax(col, row, ic_params, ic_params2=ic_params2, 
                             ctime=ctime, convert_units=convert_units, s1b_minmax_ax=ax, 
                             s1b_minmax_fig=fig,
                            savedir = savedir_rows, 
                            show_plot=show_plot)
            else: 
                print('plotting summary: ' + str(last_fig))
                s1b_minmax_fig, s1b_minmax_ax = pd.plot_icminmax_col(last_fig, col, ic_params,
                                                                 ic_params2=ic_params2, ctime=ctime,
                                                                 s1b_minmax_ax=s1b_minmax_ax,
                                                                 s1b_minmax_fig=s1b_minmax_fig, 
                                                                 convert_units=convert_units,
                                                                 show_plot=show_plot, savedir=savedir_cols,
                                                                 manual_bias_idx = manual_bias_idx)
            (ic_col, ic_min, ic_max, mod,
             optimal_bias, crosstalk_bias, manual_mod) = get_icmaxcolmod(
                ic_params, ic_params2, manual_bias=manual_bias_idx)
            ic_col_grid = fill_grid_data(ic_col, row, col, grid=ic_col_grid)
            ic_max_grid = fill_grid_data(ic_max, row, col, grid=ic_max_grid)
            ic_maxcoldiff_grid = fill_grid_data(ic_max-ic_col, row, col, grid=ic_maxcoldiff_grid)
            
            mod_grid = fill_grid_data(mod, row, col, grid=mod_grid)
            optimal_bias_grid = fill_grid_data(
                optimal_bias, row, col, grid=optimal_bias_grid)
            crosstalk_bias_grid = fill_grid_data(
                crosstalk_bias, row, col, grid=crosstalk_bias_grid)
            bias_crosstalk_diff_grid = fill_grid_data(
                optimal_bias-crosstalk_bias, row, col, grid=bias_crosstalk_diff_grid)

            manual_mod_grid = fill_grid_data(
                manual_mod, row, col, grid=manual_mod_grid)
    savedir_grids = os.path.join(savedir, 'gridplots')
    while not os.path.exists(savedir_grids):
                os.makedirs(savedir_grids)
    make_grids(rows, cols, ctime, show_plot, savedir_grids, convert_units,
               ic_col_grid, ic_max_grid, mod_grid, optimal_bias_grid, crosstalk_bias_grid,
               bias_crosstalk_diff_grid, ic_maxcoldiff_grid, manual_mod_grid)
   

def rs_driver(cfg, sa_data, sa_runfile, rsdf, rs_runfile, ctime=None,
              rsdf_off=None,  rs_runfile_off=None, filter_sq1=True,
              cols=range(0, 16), rows=range(0, 40)):
    chip_starts = [0, 10, 20, 30, 41]
    sq1_sgfilter_window_length = 5
    sq1_sgfilter_poly_deg = 2
    calc_slopes = False
    show_plot = False

    for col in cols:

        print("Analyzing Column: " + str(col))
        try:
            ssa_params = cp.calculate_ssa_parameters(
                sa_data, sa_runfile, cfg, col, show_plot=show_plot)
        except TypeError:
            print('Skipping Column: ' + str(col))
            continue
        if(ssa_params is None):
            print('Skipping Column: ' + str(col))
            continue
        s1b_minmax_fig = None
        s1b_minmax_ax = None

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

            sq1_params = cp.calculate_sq1_parameters(rsdf, rs_runfile, cfg, col, row,
                                                  ssa_params, filter_sq1=filter_sq1, calc_slopes=calc_slopes,
                                                  sq1_sgfilter_window_length=sq1_sgfilter_window_length,
                                                  sq1_sgfilter_poly_deg=sq1_sgfilter_poly_deg)

            # ic_params=calculate_icminmax(cfg, filter_sq1, row, col,  sq1_params, ssa_params,
            #                             sq1_sgfilter_window_length, sq1_sgfilter_poly_deg)
            sq1_params2 = None
            if(rsdf_off is not None):
                sq1_params2 = cp.calculate_sq1_parameters(rsdf_off, rs_runfile_off, cfg, col, row,
                                                       ssa_params, filter_sq1=filter_sq1)

                # ic_params2=calculate_icminmax(cfg, filter_sq1, row, col,  sq1_params2, ssa_params,
                #                              sq1_sgfilter_window_length, sq1_sgfilter_poly_deg)
            # else:
            #    ic_params2 = None
            # pd.plot_icminmax(sq1_safb_servo_biases_uA, sq1_safb_servo_mins_sa_in_uA, sq1_safb_servo_maxs_sa_in_uA,
             #     max_sq1imod_idx, max_sq1imod_uA)#, tune_ctime, col, row)
            s1b_minmax_fig, s1b_minmax_ax = pd.plot_rsservo_col(last_fig, col, chip_num, sq1_params,
                                                                sq1_params2=sq1_params2, ctime=ctime,
                                                                s1b_minmax_ax=s1b_minmax_ax,
                                                                s1b_minmax_fig=s1b_minmax_fig)

def save_subset(df, savename, rows = [28,29,30,31,32,33], cols = [8,9,10,11,12]):
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
    parser.add_argument('-v', '--verbose', action='store_true', 
                        help='Whether to print out debug statements')
    args = parser.parse_args()

    numcols = 32
    cols = range(0, numcols)

    print('Reading in files:' + str(args.ctime))
    ctime = os.path.basename(os.path.dirname(args.ctime))
    print('ctime: ' + str(ctime))
    cfg = rd.get_config_file()

    if(args.dev_cur):
        time0 = time.time()
        sa_data, sa_runfile = rd.get_ssa_tune_data(args.ctime)
        sq1df, sq1_runfile =  rd.get_sq1_tune_data(args.ctime)
        all_rows = sq1df['<row>'].astype(int)
        
        rows = np.unique(all_rows)
        sq1df_off = None
        sq1_runfile_off = None
        if(args.ctime_off is not None):
            sq1df_off, sq1_runfile_off = rd.get_sq1_tune_data(args.ctime_off)
            #save_subset(sq1df_off, 'rowsel_off_small_sq1servo_sa.bias')
        time1 = time.time()

        output_dir = './output_data/'
        print('Done reading files, time elapsed (s):' + str(time1-time0))
        for convert_units in [True, False]:
            if(convert_units):
                branch = 'units_ua'
            else:
                branch = 'units_dac'
            savedir = os.path.join(output_dir,ctime, branch)
            while not os.path.exists(savedir):
                os.makedirs(savedir)
            ic_driver(sq1df, sq1_runfile, ctime=ctime,
              sq1df_off=sq1df_off,  sq1_runfile_off=sq1_runfile_off,
              savedir = savedir,  plot_all_rows=True,
              convert_units=convert_units, cfg=cfg, sa_data=sa_data, sa_runfile=sa_runfile,
              verbose=args.verbose)
            ic_driver(sq1df, sq1_runfile, ctime=ctime,
              sq1df_off=sq1df_off,  sq1_runfile_off=sq1_runfile_off,
              savedir = savedir,  plot_all_rows=False,
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
                  filter_sq1=True, ctime=ctime,
                  cols=cols, rows=rows)


if __name__ == '__main__':
    time0 = time.time()
    main()
    time1 = time.time()
    print("Analysis complete, time elapsed (s): " + str(time1-time0))
