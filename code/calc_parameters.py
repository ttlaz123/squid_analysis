import numpy as np 
import scipy

#import calculate_slopes as cs 
from estimate_phi0 import estimate_phi0
from scipy.signal import savgol_filter

def get_fb_dac_to_in_ua(ssa_params, cfg):
    '''
    Converts from SAFB to SAIN to get the device current
    '''
    sa_fb_dac_to_uA, M_ssa_fb_pH = ssa_params
    turns_ratio = M_ssa_fb_pH/float(cfg['SSA']['SSA_M_IN_PICOHENRY'])
    sa_fb_dac_to_sa_in_uA = sa_fb_dac_to_uA * turns_ratio
    return sa_fb_dac_to_sa_in_uA

def get_fb_dac_to_uA(cfg):
    '''
    Get conversion from DAC units to microamps for Feedback current
    '''
     # the 1.e6 converts from A->uA
    sa_fb_dac_to_uA = 1.e6*(float(cfg['SSA']['SSA_FB_DAC_MAX_VOLTAGE_VOLTS'])/( 
        (np.power(2, int(cfg['SSA']['SSA_FB_DAC_NBITS'])))*(
            float(cfg['CRYOCABLE']['CRYOCABLE_ROUNDTRIP_RESISTANCE_OHMS']) +
            float(cfg['SSA']['SSA_FB_BACKPLANE_RESISTANCE_OHMS']) +
            float(cfg['SSA']['SSA_FB_BC_RESISTANCE_OHM']))))
    
    return sa_fb_dac_to_uA

def get_bias_dac_to_uA(cfg):
    '''
    Get conversion from DAC units to microamps for Bias current
    '''
     # the 1.e6 converts from A->uA
    sa_bias_dac_to_uA = 1.e6*float(cfg['SSA']['SSA_BIAS_DAC_MAX_VOLTAGE_VOLTS'])/(
        (np.power(2,float(cfg['SSA']['SSA_BIAS_DAC_NBITS'])))*(
        float(cfg['SSA']['SSA_BIAS_RC_RESISTANCE_OHMS'])+
        float(cfg['CRYOCABLE']['CRYOCABLE_ROUNDTRIP_RESISTANCE_OHMS'])))

    
    return sa_bias_dac_to_uA


def get_bias_fb_range(runfile):
    '''
    Input: 
            runfile -- .run file with bias and feedback run parameters
    Output:
            bias_dacs -- range of biases sweeped in DAC units
            fb_dacs -- range of fbs sweeped in DAC units
    '''
    ramp_name = 'par_ramp'
    bias_loop = 'par_step loop1 par1'
    fb_loop =  'par_step loop2 par1'

    # What range of SQ1 biases did we sweep?
    bias_start, bias_step, num_bias_steps = tuple([
        int(i) for i in runfile.Item(ramp_name, bias_loop)])
    bias_dacs = bias_start + bias_step*np.arange(num_bias_steps)

    # What range of SQ1 feedbacks did we sweep?
    fb_start, fb_step, num_fb_steps = tuple([
        int(i) for i in runfile.Item(ramp_name, fb_loop)])
    fb_dacs = fb_start + fb_step*np.arange(num_fb_steps)

    return bias_dacs, fb_dacs

def get_rms_noise(df, bias_colname, fb_colname, zero_bias = 0):
    '''
    gets the rms noise at zero bias
    '''
    sq1df_rowbias =  df[df[bias_colname]==zero_bias]
    sq1_safb = sq1df_rowbias[fb_colname].values
    sq1_safb_mean = np.mean(sq1_safb)
    rms = np.sqrt(np.mean(np.square(sq1_safb-sq1_safb_mean)))
    return rms


def calculate_ic_params(sq1_rowcol_df, sq1_runfile, col, mod_thresh = 20,
                        convert_units = False, cfg=None, ssa_params=None, 
                        filter_sq1=True, sq1_sgfilter_window_length=5, sq1_sgfilter_poly_deg=2):
    '''
    Calculates the important parameters for each row/col squid
    Input: 
            sq1_rowcol_df -- dataframe containing bias and fb data for the row and columns
            sq1_runfile -- file containing info about the run, used to get bias and fb DAC ranges
            col -- column number
            mod_thresh -- threshold for determining whether modulation has occurred,
                            as multiple of rms noise at zero bias
            convert_units -- whether to convert from DAC to microamps
            cfg -- config file for calibrating DAC to microamps
            ssa_params -- parameters for calibrating DAC to microamps
            filter parameters -- for applying the savgol filters
    Output: An IC_params dictionary with 
            'fb_min' -- list of lower bound for fb curve 
            'fb_max' -- list of upper bound for fb curve
            'bias_max_idx' -- index of bias point of maximum modulation
            'bias_min_idx' -- index of bias point of first modulation
            
    '''
    bias_colname = '<bias>'
    fb_colname = f'<safb{col:02}>'

    bias_dacs, fb_dacs = get_bias_fb_range(sq1_runfile)
    icmin_thresh = mod_thresh * get_rms_noise(sq1_rowcol_df, bias_colname, fb_colname)

    
    # Getting bias currents at min and max modulations
    # Min is min bias above the mod_thresh
    sq1_safb_curves = []
    sq1_biases = []

    max_sq1_safb_mod = -1
    max_sq1_safb_mod_bias_idx = -1
    min_sq1_safb_mod_bias_idx = len(bias_dacs)-1

    for b_index, grp in sq1_rowcol_df.groupby(bias_colname):
        bias = bias_dacs[b_index]
        sq1_safb = grp[fb_colname].values
        sq1_safb_curves.append(sq1_safb)
        sq1_biases.append(bias)
        
        sq1_safb_span = np.max(sq1_safb)-np.min(sq1_safb)

        if sq1_safb_span > max_sq1_safb_mod:
            max_sq1_safb_mod_bias_idx = b_index
            max_sq1_safb_mod = sq1_safb_span
        
        if(sq1_safb_span > icmin_thresh and b_index < min_sq1_safb_mod_bias_idx):
            min_sq1_safb_mod_bias_idx = b_index

    # np.max/np.min reversed because of SQUID coil polarity; must flip to get physical current
    sq1_safb_mins_dac = np.array([np.max(sq1_safb_curve_dac)
                                        for sq1_safb_curve_dac in sq1_safb_curves])
    sq1_safb_maxs_dac = np.array([np.min(sq1_safb_curve_dac)
                                        for sq1_safb_curve_dac in sq1_safb_curves])
    
     # To convert to current, need to flip and zero starting value, then
    # scale to SSA in current units

    safb_zero_offset = np.mean(sq1_safb_curves[0])
    sq1_safb_servo_mins = (safb_zero_offset-sq1_safb_mins_dac)
    sq1_safb_servo_maxs = (safb_zero_offset-sq1_safb_maxs_dac)

    if(convert_units):
        sq1_bias_dac_to_uA =  get_bias_dac_to_uA(cfg)
        sa_fb_dac_to_sa_in_uA = get_fb_dac_to_in_ua(cfg, ssa_params)
        sq1_biases = np.array(sq1_biases)*sq1_bias_dac_to_uA
        sq1_safb_servo_mins = np.array(sq1_safb_servo_mins)*sa_fb_dac_to_sa_in_uA
        sq1_safb_servo_maxs = np.array(sq1_safb_servo_maxs)*sa_fb_dac_to_sa_in_uA

        
    if filter_sq1:
        sq1_safb_servo_mins = savgol_filter(
            sq1_safb_servo_mins, sq1_sgfilter_window_length, sq1_sgfilter_poly_deg)
        sq1_safb_servo_maxs = savgol_filter(
            sq1_safb_servo_maxs, sq1_sgfilter_window_length, sq1_sgfilter_poly_deg)



 

    ic_params = {'bias': sq1_biases,
                'fb_min':sq1_safb_servo_mins, 
                 'fb_max':sq1_safb_servo_maxs, 
                 'bias_max_idx':max_sq1_safb_mod_bias_idx, 
                 'bias_min_idx': min_sq1_safb_mod_bias_idx
                 }

    return ic_params


def calculate_ssa_parameters(sa_data, sa_runfile, cfg, col):
    '''
    Takes in the ssa data to plot all the parameters for each ssa column 
    Adapted from David Goldfinger's script
    '''
    sa_bias = sa_runfile.Item('HEADER', f'RB sa bias', type='int')
    # Which SSAs were biased with nonzero bias?
    sa_biased_columns = np.nonzero(sa_bias)[0]
    # What range of SSA feedback did we sweep?
    sa_fb0, d_sa_fb, n_sa_fb = tuple([
        int(i) for i in sa_runfile.Item('par_ramp', 'par_step loop1 par1')])
    sa_fb = sa_fb0 + d_sa_fb*np.arange(n_sa_fb)
    # Convert to physical current using configuration file info
    
    sa_fb_dac_to_uA = get_fb_dac_to_uA(cfg)
    sa_fb = np.multiply(sa_fb, sa_fb_dac_to_uA)
    nrow, ncol, sa_npt = sa_data.data.shape
    sa_coadd_adu = np.average([sa_data.data[row][col]
                               for row in range(nrow)], axis=0)

    rc = int(col/8)+1
    sample_num = sa_runfile.Item(
        'HEADER', f'RB rc{rc} sample_num', type='int')[0]
    sa_adu = sa_coadd_adu/sample_num
    # Convert sa_fb and sa_adu to physical units.
    sa_bias_dac = sa_bias[col]
    sa_bias_dac_to_uA = get_bias_dac_to_uA(cfg)
    sa_bias_ua = sa_bias_dac*sa_bias_dac_to_uA
    # Estimate the SSA phi0 from autocorrelation
    try:
        sa_phi0_est = estimate_phi0(sa_fb, sa_adu, debug=False)
    except ValueError:
        #pd.plot_ssa_debug(sa_fb, sa_adu)
        return
    M_ssa_fb_pH = 1.e12 * \
        scipy.constants.value(u'mag. flux quantum')/(sa_phi0_est*1.e-6)
    #if(show_plot):
    #    sa_ax = cs.calc_ssa_slopes(sa_adu, sa_fb, sa_phi0_est, cfg, d_sa_fb, sa_fb_dac_to_uA,
    #                col, sa_bias_ua, M_ssa_fb_pH, show_plot)
    # DONE PLOTTING SSA
    ssa_params = sa_fb_dac_to_uA, M_ssa_fb_pH
    return ssa_params
#################
# Everything below this is 
# DEPRECATED
####################


def calculate_sq1_parameters(sq1df, sq1_runfile, cfg, col, row, ssa_params,
                             filter_sq1=True, sq1_sgfilter_window_length=5, calc_slopes=False,
                             sq1_sgfilter_poly_deg=2, s1_slope_npts_to_fit=6):
    '''
    Takes in the sq1 data to plot all the parameters for each sq1 row 
    Adapted from David Goldfinger's script
    '''
    sa_fb_dac_to_uA, M_ssa_fb_pH = ssa_params
    sq1_b0, d_sq1_b, n_sq1_b = tuple([
        int(i) for i in sq1_runfile.Item('par_ramp', 'par_step loop1 par1')])
    sq1_b = sq1_b0 + d_sq1_b*np.arange(n_sq1_b)
    # What range of SQ1 feedbacks did we sweep?
    sq1_fb0, d_sq1_fb, n_sq1_fb = tuple([
        int(i) for i in sq1_runfile.Item('par_ramp', 'par_step loop2 par1')])
    sq1_fb = sq1_fb0 + d_sq1_fb*np.arange(n_sq1_fb)

    # Find bias corresponding to max SQ1 response
    max_sq1_safb_servo_span = 0
    max_sq1_safb_servo_span_bias = None
    sq1_safb_servo_curves_dac = []
    sq1_safb_servo_biases_dac = []
    sq1df_row = sq1df[(sq1df['<row>'] == row)]
    
    if(sq1df_row.shape[0] < 5):
        raise ValueError("Not enough data for row " + str(row), '. Skipping row.')
    for b_index, grp in sq1df_row.groupby('<bias>'):
        b = sq1_b[b_index]
        sq1_safb_servo = grp[f'<safb{col:02}>'].values
        sq1_safb_servo_curves_dac.append(sq1_safb_servo)
        sq1_safb_servo_biases_dac.append(b)
        sq1_safb_servo_span = np.max(sq1_safb_servo)-np.min(sq1_safb_servo)
        if sq1_safb_servo_span > max_sq1_safb_servo_span:
            max_sq1_safb_servo_span = sq1_safb_servo_span
            max_sq1_safb_servo_span_bias = b
    # print(
    #    f'max_sq1_safb_servo_span_bias={max_sq1_safb_servo_span_bias}\tmax_sq1_safb_servo_span={max_sq1_safb_servo_span}')

    #if(calc_slopes):
    #    print("Calculating Slopes")
    #    cs.calculate_sq1_slopes(cfg, sq1df, sq1_fb, sq1_b, max_sq1_safb_servo_span_bias, row, col,
    #                         sa_fb_dac_to_uA, s1_slope_npts_to_fit, M_ssa_fb_pH,
    #                         filter_sq1, sq1_sgfilter_window_length=sq1_sgfilter_window_length, sq1_sgfilter_poly_deg=sq1_sgfilter_poly_deg)

    sq1_params = (sq1_safb_servo_curves_dac, sq1_safb_servo_biases_dac,
                  max_sq1_safb_servo_span_bias, max_sq1_safb_servo_span, sq1_fb
                  )

    return sq1_params
