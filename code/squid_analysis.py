'''
Written by Tom Liu, adapting code from David Goldfinger
2023 May 
'''
import os
import argparse
import time

import numpy as np


import scipy
from scipy import stats
from scipy.signal import savgol_filter

import read_data as rd
import plot_data as pd


# we should figure out where these warnings are coming from some day
import warnings
warnings.filterwarnings("ignore")


def serial_corr(wave, lag=1):
    n = len(wave)
    y1 = wave[lag:]
    y2 = wave[:n-lag]
    corr = np.corrcoef(y1, y2)[0, 1]
    return corr


def autocorr(wave):
    lags = range(len(wave))
    corrs = np.array([serial_corr(wave, lag) for lag in lags])
    return lags, corrs


def estimate_phi0(phi, squid_curve, min_acorr_dist_from_zero_frac=0.1, debug=False):
    """Estimate SQUID curve phi0
    phi : numpy.ndarray
       Array of flux bias.
    squid_curve : numpy.ndarray
       Array of SQUID voltages, one for each flux bias.
    min_acorr_dist_from_zero_frac : float, optional, default 0.1
       ???
    debug : bool, optional, default False
       Whether or not to show debug plots/print-out.

    Returns
    -------
    phi0 : float
       Estimated phi0.  Returns None if unable to estimate
       the phi0 using lag and autocorrelation.
    """
    min_acorr_dist_from_zero = len(phi)*min_acorr_dist_from_zero_frac

    if debug:
        print(
            f'-> min_acorr_dist_from_zero_frac={min_acorr_dist_from_zero_frac}')
        print(f'-> min_acorr_dist_from_zero={min_acorr_dist_from_zero}')

    # find period from autocorrelation
    lags, corrs = autocorr(squid_curve)

    # find peaks in autocorrelation vs lag
    from scipy.signal import find_peaks
    peaks, _ = find_peaks(corrs, height=0)

    sorted_peaks = sorted([pk for _, pk in zip(corrs[peaks], peaks)])

    if debug:
        print(f'-> sorted_peaks[:4]={sorted_peaks[:4]}')

    try:
        phi0_idx = next(pk for pk in sorted_peaks if pk >
                        min_acorr_dist_from_zero)
    except:
        raise ValueError('No peaks found')

    phi0 = np.abs(phi[phi0_idx]-phi[0])

    if debug:
        print(f'-> phi0={phi0}')

    return phi0


def get_midpoints(x):
    # X should be a numpy.ndarray
    return (x[1:] + x[:-1])/2.


def calculate_ssa_parameters(sa_data, sa_runfile, cfg, col, show_plot=False):
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
    sa_fb_dac_to_uA = 1.e6*(float(cfg['SSA']['SSA_FB_DAC_MAX_VOLTAGE_VOLTS'])/(  # the 1.e6 converts from A->uA
        (np.power(2, int(cfg['SSA']['SSA_FB_DAC_NBITS'])))*(
            float(cfg['CRYOCABLE']['CRYOCABLE_ROUNDTRIP_RESISTANCE_OHMS']) +
            float(cfg['SSA']['SSA_FB_BACKPLANE_RESISTANCE_OHMS']) +
            float(cfg['SSA']['SSA_FB_BC_RESISTANCE_OHM']))))
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
    sa_bias_ua = 1.e6*sa_bias_dac*float(cfg['SSA']['SSA_BIAS_DAC_MAX_VOLTAGE_VOLTS'])/(
        (2**(float(cfg['SSA']['SSA_BIAS_DAC_NBITS'])))*(float(cfg['SSA']['SSA_BIAS_RC_RESISTANCE_OHMS'])+float(cfg['CRYOCABLE']['CRYOCABLE_ROUNDTRIP_RESISTANCE_OHMS'])))

    # Estimate the SSA phi0 from autocorrelation
    try:
        sa_phi0_est = estimate_phi0(sa_fb, sa_adu, debug=False)
    except ValueError:
        pd.plot_ssa_debug(sa_fb, sa_adu)
        return
    # Find the maximum in the first phi0
    sa_max_idx = np.argmax(
        sa_adu[np.where(sa_fb < (np.min(sa_fb)+sa_phi0_est))])

    # Find max and min slope in first full phi0
    sa_fb_midpoints = get_midpoints(sa_fb)
    sa_adu_midpoints = get_midpoints(sa_adu)
    sa_adu_diff = np.diff(sa_adu)
    # Max upslope
    sa_adu_max_upslope_idx = np.argmax(
        sa_adu_diff[np.where(sa_fb < sa_phi0_est)])
    # Max downslope
    sa_adu_max_downslope_idx = np.argmin(
        sa_adu_diff[np.where(sa_fb < sa_phi0_est)])

    dVdADU = float(cfg['PREAMPADC']['ADU_TO_VOLTS_AT_PREAMP_INPUT'])

    Vmod_adu = np.max(sa_adu)-np.min(sa_adu)
    Vmod_mV = 1.e3*Vmod_adu*dVdADU
    M_ssa_fb_pH = 1.e12 * \
        scipy.constants.value(u'mag. flux quantum')/(sa_phi0_est*1.e-6)
    m_sa_adu_upslope = sa_adu_diff[sa_adu_max_upslope_idx] / \
        (d_sa_fb*sa_fb_dac_to_uA)
    b_sa_adu_upslope = (sa_adu_midpoints[sa_adu_max_upslope_idx] -
                        m_sa_adu_upslope*sa_fb_midpoints[sa_adu_max_upslope_idx])

    # Downslopes
    m_sa_adu_downslope = sa_adu_diff[sa_adu_max_downslope_idx] / \
        (d_sa_fb*sa_fb_dac_to_uA)
    b_sa_adu_downslope = (sa_adu_midpoints[sa_adu_max_downslope_idx] -
                          m_sa_adu_downslope*sa_fb_midpoints[sa_adu_max_downslope_idx])
    dV_ADU_dSAFB_DAC_downslope = m_sa_adu_downslope*sa_fb_dac_to_uA
    dV_nV_dSAFB_uphi0_downslope = m_sa_adu_downslope * (
        # converts from ADU/uA -> V/uA
        float(cfg['PREAMPADC']['ADU_TO_VOLTS_AT_PREAMP_INPUT']) *
        # converts from V/uA to V/phi0
        sa_phi0_est *
        # converts from V/phi0 to nV/uphi0
        1.e6/1000.
    )
    dV_nV_dSAIN_pA_downslope = dV_nV_dSAFB_uphi0_downslope / (
        1.e6*scipy.constants.value(u'mag. flux quantum') /
        (float(cfg['SSA']['SSA_M_IN_PICOHENRY'])*1.e-12)
    )

    # Upslopes
    dV_ADU_dSAFB_DAC_upslope = m_sa_adu_upslope*sa_fb_dac_to_uA
    dV_nV_dSAFB_uphi0_upslope = m_sa_adu_upslope * (
        # converts from ADU/uA -> V/uA
        float(cfg['PREAMPADC']['ADU_TO_VOLTS_AT_PREAMP_INPUT']) *
        # converts from V/uA to V/phi0
        sa_phi0_est *
        # converts from V/phi0 to nV/uphi0
        1.e6/1000.
    )
    dV_nV_dSAIN_pA_upslope = dV_nV_dSAFB_uphi0_upslope / (
        1.e6*scipy.constants.value(u'mag. flux quantum') /
        (float(cfg['SSA']['SSA_M_IN_PICOHENRY'])*1.e-12)
    )

    sa_ax = pd.plot_ssa(cfg, sa_fb, sa_adu, sa_phi0_est, sa_max_idx,
                        sa_fb_midpoints, sa_adu_midpoints,
                        sa_adu_max_upslope_idx, sa_adu_max_downslope_idx,
                        sa_fb_dac_to_uA, b_sa_adu_upslope, m_sa_adu_upslope,
                        b_sa_adu_downslope, m_sa_adu_downslope,
                        dV_ADU_dSAFB_DAC_downslope, dV_nV_dSAFB_uphi0_downslope, dV_nV_dSAIN_pA_downslope,
                        dV_ADU_dSAFB_DAC_upslope, dV_nV_dSAFB_uphi0_upslope, dV_nV_dSAIN_pA_upslope,
                        col, sa_bias_ua, Vmod_mV, M_ssa_fb_pH, show_plot=show_plot)

    # DONE PLOTTING SSA
    ssa_params = sa_fb_dac_to_uA, M_ssa_fb_pH, sa_ax
    return ssa_params


def calculate_sq1_slopes(cfg, sq1df, sq1_fb, sq1_b, max_sq1_safb_servo_span_bias, row, col,
                         sa_fb_dac_to_uA, s1_slope_npts_to_fit, M_ssa_fb_pH, sa_ax,
                         filter_sq1, sq1_sgfilter_window_length=5, sq1_sgfilter_poly_deg=2):
    '''
    Calculates the up and downslope

    '''
    max_sq1_safb_servo = sq1df[(sq1df['<row>'] == row) & (
        sq1df['<bias>'] == np.where(sq1_b == max_sq1_safb_servo_span_bias)[0][0])][f'<safb{col:02}>'].values

    # Convert SSA fb and SQ1 fb to physical units
    sq1fb_dac_to_nA = 1.e9*float(cfg['SQ1']['SQ1_FB_DAC_MAX_CURRENT_A'])/(
        2**(int(cfg['SQ1']['SQ1_FB_DAC_NBITS'])))*(
            float(cfg['SQ1']['SQ1_FB_RC_R33_OHM'])/(
                float(cfg['SQ1']['SQ1_FB_RC_R33_OHM']) +
                float(cfg['SQ1']['SQ1_FB_RC_R28_OHM']) +
                float(cfg['CRYOCABLE']['CRYOCABLE_ROUNDTRIP_RESISTANCE_OHMS']) +
                float(cfg['SQ1']['SQ1_FB_BACKPLANE_RESISTANCE_OHMS'])
            )
    )

    sq1_fb_uA = sq1fb_dac_to_nA*sq1_fb/1000.
    max_sq1_safb_servo_uA_unfilt = max_sq1_safb_servo*sa_fb_dac_to_uA

    max_sq1_safb_servo_uA = max_sq1_safb_servo_uA_unfilt

    if filter_sq1:
        max_sq1_safb_servo_uA = savgol_filter(
            max_sq1_safb_servo_uA, sq1_sgfilter_window_length, sq1_sgfilter_poly_deg)

    # Estimate the SSA phi0 from autocorrelation
    s1_phi0_est = estimate_phi0(sq1_fb_uA, max_sq1_safb_servo_uA, debug=False)
    M_s1_fb = 1.e12 * \
        scipy.constants.value(u'mag. flux quantum')/(s1_phi0_est*1.e-6)

    # Find the maximum in the first phi0
    s1_max_idx = np.argmax(
        max_sq1_safb_servo_uA[np.where(sq1_fb_uA < (np.min(sq1_fb_uA)+s1_phi0_est))])

    # Find max and min slope in first full phi0
    sq1_fb_uA_midpoints = get_midpoints(sq1_fb_uA)
    max_sq1_safb_servo_uA_midpoints = get_midpoints(max_sq1_safb_servo_uA)
    max_sq1_safb_servo_uA_diff = np.diff(max_sq1_safb_servo_uA)
    # Max upslope
    max_sq1_safb_servo_uA_max_upslope_idx = np.argmax(
        max_sq1_safb_servo_uA_diff[np.where(sq1_fb_uA < (np.min(sq1_fb_uA)+s1_phi0_est))])
    # Max downslope
    max_sq1_safb_servo_uA_max_downslope_idx = np.argmin(
        max_sq1_safb_servo_uA_diff[np.where(sq1_fb_uA < (np.min(sq1_fb_uA)+s1_phi0_est))])

    # SQ1 curves tend to be a little noisier.  Fit neighboring points
    s1_fb_uA_subset_for_downslope_fit = sq1_fb_uA[max_sq1_safb_servo_uA_max_downslope_idx -
                                                  int(s1_slope_npts_to_fit/2-1):max_sq1_safb_servo_uA_max_downslope_idx+int(s1_slope_npts_to_fit/2+1)]
    max_sq1_safb_servo_uA_subset_for_downslope_fit = max_sq1_safb_servo_uA[max_sq1_safb_servo_uA_max_downslope_idx - int(
        s1_slope_npts_to_fit/2-1):max_sq1_safb_servo_uA_max_downslope_idx+int(s1_slope_npts_to_fit/2+1)]

    m_max_sq1_safb_servo_uA_downslope = None
    try:
        s1_downslope_p1 = np.polyfit(
            s1_fb_uA_subset_for_downslope_fit, max_sq1_safb_servo_uA_subset_for_downslope_fit, 1)
        m_max_sq1_safb_servo_uA_downslope = s1_downslope_p1[0]
        b_max_sq1_safb_servo_uA_downslope = s1_downslope_p1[1]

    except:
        print("Failed to fit SQ1 downslope")
        m_max_sq1_safb_servo_uA_downslope = 0
        b_max_sq1_safb_servo_uA_downslope = 0

    # SQ1 curves tend to be a little noisier.  Fit neighboring points
    s1_fb_uA_subset_for_upslope_fit = sq1_fb_uA[max_sq1_safb_servo_uA_max_upslope_idx -
                                                int(s1_slope_npts_to_fit/2-1):max_sq1_safb_servo_uA_max_upslope_idx+int(s1_slope_npts_to_fit/2+1)]
    max_sq1_safb_servo_uA_subset_for_upslope_fit = max_sq1_safb_servo_uA[max_sq1_safb_servo_uA_max_upslope_idx - int(
        s1_slope_npts_to_fit/2-1):max_sq1_safb_servo_uA_max_upslope_idx+int(s1_slope_npts_to_fit/2+1)]
    s1_upslope_p1 = np.polyfit(
        s1_fb_uA_subset_for_upslope_fit, max_sq1_safb_servo_uA_subset_for_upslope_fit, 1)
    m_max_sq1_safb_servo_uA_upslope = s1_upslope_p1[0]
    b_max_sq1_safb_servo_uA_upslope = s1_upslope_p1[1]

    if m_max_sq1_safb_servo_uA_downslope is not None:
        dI_SSA_IN_pA_dI_SQ1_IN_pA_downslope = m_max_sq1_safb_servo_uA_downslope*(
            # converts SAFB->SAIN
            1.e-6 * 1.e12 * (M_ssa_fb_pH/float(cfg['SSA']['SSA_M_IN_PICOHENRY'])))/(
            1.e-6 * 1.e12 * (M_s1_fb/float(cfg['SQ1']['SQ1_M_IN_PICOHENRY'])))  # converts S1FB->S1IN
    else:
        dI_SSA_IN_pA_dI_SQ1_IN_pA_downslope = None

    dI_SSA_IN_pA_dI_SQ1_IN_pA_upslope = m_max_sq1_safb_servo_uA_upslope*(
        # converts SAFB->SAIN
        1.e-6 * 1.e12 * (M_ssa_fb_pH/float(cfg['SSA']['SSA_M_IN_PICOHENRY'])))/(
            1.e-6 * 1.e12 * (M_s1_fb/float(cfg['SQ1']['SQ1_M_IN_PICOHENRY'])))  # converts S1FB->S1IN

    pd.plot_sq1(col, row, sq1_fb_uA, max_sq1_safb_servo_uA, filter_sq1,
                max_sq1_safb_servo_uA_unfilt,
                s1_max_idx, sq1_fb_uA_midpoints, max_sq1_safb_servo_uA_midpoints,
                max_sq1_safb_servo_uA_max_downslope_idx, s1_fb_uA_subset_for_downslope_fit,
                max_sq1_safb_servo_uA_subset_for_downslope_fit,
                m_max_sq1_safb_servo_uA_downslope, b_max_sq1_safb_servo_uA_downslope,

                max_sq1_safb_servo_uA_max_upslope_idx, s1_fb_uA_subset_for_upslope_fit,
                max_sq1_safb_servo_uA_subset_for_upslope_fit,
                b_max_sq1_safb_servo_uA_upslope, m_max_sq1_safb_servo_uA_upslope,
                s1_phi0_est, M_s1_fb,
                dI_SSA_IN_pA_dI_SQ1_IN_pA_downslope, dI_SSA_IN_pA_dI_SQ1_IN_pA_upslope,
                sa_fb_dac_to_uA, sa_ax, cfg)


def calculate_sq1_parameters(sq1df, sq1_runfile, cfg, col, row, ssa_params,
                             filter_sq1=True, sq1_sgfilter_window_length=5, calc_slopes=False,
                             sq1_sgfilter_poly_deg=2, s1_slope_npts_to_fit=6):
    '''
    Takes in the sq1 data to plot all the parameters for each sq1 row 
    Adapted from David Goldfinger's script
    '''
    sa_fb_dac_to_uA, M_ssa_fb_pH, sa_ax = ssa_params
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

    if(calc_slopes):
        print("Calculating Slopes")
        calculate_sq1_slopes(cfg, sq1df, sq1_fb, sq1_b, max_sq1_safb_servo_span_bias, row, col,
                             sa_fb_dac_to_uA, s1_slope_npts_to_fit, M_ssa_fb_pH, sa_ax,
                             filter_sq1, sq1_sgfilter_window_length=sq1_sgfilter_window_length, sq1_sgfilter_poly_deg=sq1_sgfilter_poly_deg)

    sq1_params = (sq1_safb_servo_curves_dac, sq1_safb_servo_biases_dac,
                  max_sq1_safb_servo_span_bias, max_sq1_safb_servo_span, sq1_fb
                  )

    return sq1_params


def calculate_icminmax(cfg, filter_sq1, row, col,  sq1_params, ssa_params,
                       sq1_sgfilter_window_length, sq1_sgfilter_poly_deg, mod_thresh=0.02):
    '''
    Caluculates the ic min and ic max of each row
    '''
    
    (sq1_safb_servo_curves_dac, sq1_safb_servo_biases_dac,
     max_sq1_safb_servo_span_bias, max_sq1_safb_servo_span, sq1_safb_servo
     ) = sq1_params
    sa_fb_dac_to_uA, M_ssa_fb_pH, sa_ax = ssa_params
    assert len(sq1_safb_servo_biases_dac
               )>1, "Must have more than 1 bias point: bias is currently" + str(sq1_safb_servo_biases_dac)
    # np.max/np.min reversed because of SQUID coil polarity; must flip to get physical current
    sq1_safb_servo_mins_dac = np.array([np.max(sq1_safb_servo_curve_dac)
                                        for sq1_safb_servo_curve_dac in sq1_safb_servo_curves_dac])
    sq1_safb_servo_maxs_dac = np.array([np.min(sq1_safb_servo_curve_dac)
                                        for sq1_safb_servo_curve_dac in sq1_safb_servo_curves_dac])

    # To convert to current, need to flip and zero starting value, then
    # scale to SSA in current units
    sq1_bias_dac_to_uA = 1.e6*(float(cfg['SQ1']['SQ1_BIAS_DAC_MAX_VOLTAGE_VOLTS'])/(  # the 1.e6 converts from A->uA
        (np.power(2, int(cfg['SQ1']['SQ1_BIAS_DAC_NBITS'])))*(
            float(cfg['CRYOCABLE']['CRYOCABLE_ROUNDTRIP_RESISTANCE_OHMS']) +
            float(cfg['SQ1']['SQ1_BIAS_BACKPLANE_RESISTANCE_OHMS']) +
            float(cfg['SQ1']['SQ1_BIAS_BC_RESISTANCE_OHM']))))

    sq1_safb_servo_biases_uA = np.array(
        sq1_safb_servo_biases_dac)*sq1_bias_dac_to_uA

    # sa_fb_servo_maxs_uA =  sa_fb_dac_to_uA*(
    #    np.mean(sq1_safb_servo_curves_dac[0])-sq1_safb_servo_maxs_dac) #The max and min curve should be identical where it is superconducting, but the max curve is superconducting longer.

    #min_sc_branch_current_uA = 1
    #max_sc_branch_current_uA = 5.5
    # indices_sc = np.where(
    #    (sq1_safb_servo_biases_uA > min_sc_branch_current_uA) &
    #    (sq1_safb_servo_biases_uA < max_sc_branch_current_uA))
    #zsc = np.polyfit(sq1_safb_servo_biases_uA[indices_sc], sa_fb_servo_maxs_uA[indices_sc], 1)

    #M_ssa_in_pH = M_ssa_fb_pH/zsc[0]
    # print(M_ssa_in_pH)
    # print(float(cfg['SSA']['SSA_M_IN_PICOHENRY']))

    sa_fb_dac_to_sa_in_uA = (
        sa_fb_dac_to_uA*M_ssa_fb_pH/float(cfg['SSA']['SSA_M_IN_PICOHENRY']))
    sq1_safb_servo_mins_sa_in_uA = sa_fb_dac_to_sa_in_uA*(
        np.mean(sq1_safb_servo_curves_dac[0])-sq1_safb_servo_mins_dac)
    sq1_safb_servo_maxs_sa_in_uA = sa_fb_dac_to_sa_in_uA*(
        np.mean(sq1_safb_servo_curves_dac[0])-sq1_safb_servo_maxs_dac)

    if filter_sq1:
        sq1_safb_servo_mins_sa_in_uA = savgol_filter(
            sq1_safb_servo_mins_sa_in_uA, sq1_sgfilter_window_length, sq1_sgfilter_poly_deg)
        sq1_safb_servo_maxs_uA = savgol_filter(
            sq1_safb_servo_maxs_sa_in_uA, sq1_sgfilter_window_length, sq1_sgfilter_poly_deg)

    # plot max mod
    sq1_safb_servo_span_sa_in_uA = (sq1_safb_servo_maxs_sa_in_uA -
                                    sq1_safb_servo_mins_sa_in_uA)
    max_sq1imod_idx = np.argmax(sq1_safb_servo_span_sa_in_uA)
    max_sq1imod_uA = (sq1_safb_servo_span_sa_in_uA)[max_sq1imod_idx]

    current_threshold = mod_thresh*sq1_safb_servo_maxs_sa_in_uA[max_sq1imod_idx]
    # print(sq1_safb_servo_span_sa_in_uA)
    # print(current_threshold)
    start_sq1imod_idx = np.argmax(
        sq1_safb_servo_span_sa_in_uA > current_threshold)
    if(start_sq1imod_idx == 0):
        start_sq1imod_idx = -1
    # print(start_sq1imod_idx)
    start_sq1imod_uA = (sq1_safb_servo_maxs_sa_in_uA)[start_sq1imod_idx]

    # tasf.write(
    #    f'{col}\t{row}\t{max_sq1_safb_servo_span_bias}\t{sq1_safb_servo_biases_uA[max_sq1imod_idx]: .1f}\t{max_sq1_safb_servo_span}\t{max_sq1imod_uA:.3f}\n')

    # subtract zero offsets

    # Close tune analysis summary file
    # tasf.close()

    ic_params = (sq1_safb_servo_biases_uA, sq1_safb_servo_mins_sa_in_uA,
                 sq1_safb_servo_maxs_sa_in_uA, max_sq1imod_idx, max_sq1imod_uA,
                 start_sq1imod_idx, start_sq1imod_uA)

    return ic_params


def fill_grid_data(value, row, col, grid=None, max_rows=41, max_cols=32):
    '''
    Fills in grid with a value
    '''
    if(grid is None):
        grid = np.zeros((max_rows, max_cols))
    grid[row, col] = value
    return grid


def get_icmaxcolmod(ic_params, ic_params2):
    '''
    Gets ic col, ic max, and modulation given the ic_params
    '''
    (sq1_safb_servo_biases_uA, sq1_safb_servo_mins_sa_in_uA,
     sq1_safb_servo_maxs_sa_in_uA, max_sq1imod_idx, max_sq1imod_uA,
     start_sq1imod_idx, start_sq1imod_uA) = ic_params

    ic_min = sq1_safb_servo_mins_sa_in_uA[start_sq1imod_idx]
    ic_max = sq1_safb_servo_maxs_sa_in_uA[max_sq1imod_idx]
    mod = -(sq1_safb_servo_mins_sa_in_uA[max_sq1imod_idx] -
            sq1_safb_servo_maxs_sa_in_uA[max_sq1imod_idx])
    optimal_bias = sq1_safb_servo_biases_uA[max_sq1imod_idx]

    ic_col = -1
    if(ic_params2 is not None):
        (sq1_safb_servo_biases_uA, sq1_safb_servo_mins_sa_in_uA,
         sq1_safb_servo_maxs_sa_in_uA, max_sq1imod_idx, max_sq1imod_uA,
         start_sq1imod_idx, start_sq1imod_uA) = ic_params2

        ic_col = sq1_safb_servo_maxs_sa_in_uA[start_sq1imod_idx]
        crosstalk_bias = sq1_safb_servo_biases_uA[start_sq1imod_idx]

    return ic_col, ic_min, ic_max, mod, optimal_bias, crosstalk_bias


def ic_driver(cfg, sa_data, sa_runfile, sq1df, sq1_runfile, ctime=None,
              sq1df_off=None,  sq1_runfile_off=None, filter_sq1=True,
              cols=range(0, 16), rows=range(0, 40)):
    sq1_sgfilter_window_length = 5
    sq1_sgfilter_poly_deg = 2
    calc_slopes = False
    show_plot = False
    ic_max_grid = None
    ic_col_grid = None
    mod_grid = None
    optimal_bias_grid = None
    crosstalk_bias_grid = None
    for col in cols:

        print("Analyzing Column: " + str(col))
        try:
            ssa_params = calculate_ssa_parameters(
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
            last_fig = (row == rows[-1])

            sq1_params = calculate_sq1_parameters(sq1df, sq1_runfile, cfg, col, row,
                                                  ssa_params, filter_sq1=filter_sq1, calc_slopes=calc_slopes,
                                                  sq1_sgfilter_window_length=sq1_sgfilter_window_length,
                                                  sq1_sgfilter_poly_deg=sq1_sgfilter_poly_deg)

            ic_params = calculate_icminmax(cfg, filter_sq1, row, col,  sq1_params, ssa_params,
                                           sq1_sgfilter_window_length, sq1_sgfilter_poly_deg)
            if(sq1df_off is not None):
                sq1_params2 = calculate_sq1_parameters(sq1df_off, sq1_runfile_off, cfg, col, row,
                                                       ssa_params, filter_sq1=filter_sq1)

                ic_params2 = calculate_icminmax(cfg, filter_sq1, row, col,  sq1_params2, ssa_params,
                                                sq1_sgfilter_window_length, sq1_sgfilter_poly_deg)
            else:
                ic_params2 = None
            # pd.plot_icminmax(sq1_safb_servo_biases_uA, sq1_safb_servo_mins_sa_in_uA, sq1_safb_servo_maxs_sa_in_uA,
             #     max_sq1imod_idx, max_sq1imod_uA)#, tune_ctime, col, row)
            s1b_minmax_fig, s1b_minmax_ax = pd.plot_icminmax_col(last_fig, col, ic_params,
                                                                 ic_params2=ic_params2, ctime=ctime,
                                                                 s1b_minmax_ax=s1b_minmax_ax,
                                                                 s1b_minmax_fig=s1b_minmax_fig)
            (ic_col, ic_min, ic_max, mod,
             optimal_bias, crosstalk_bias) = get_icmaxcolmod(
                ic_params, ic_params2)
            ic_col_grid = fill_grid_data(ic_col, row, col, grid=ic_col_grid)
            ic_max_grid = fill_grid_data(ic_max, row, col, grid=ic_max_grid)
            mod_grid = fill_grid_data(mod, row, col, grid=mod_grid)
            optimal_bias_grid = fill_grid_data(
                optimal_bias, row, col, grid=optimal_bias_grid)
            crosstalk_bias_grid = fill_grid_data(
                crosstalk_bias, row, col, grid=crosstalk_bias_grid)

    pd.tile_plot(len(rows), len(cols), ic_col_grid,
                 'Ic,col (uA)', str(ctime)+'_Ic_col')
    pd.tile_plot(len(rows), len(cols), ic_max_grid,
                 'Ic,max (uA)', str(ctime)+'_Ic_max')
    pd.tile_plot(len(rows), len(cols), mod_grid,
                 'Modulation (uA)', str(ctime)+'_mod')
    pd.tile_plot(len(rows), len(cols), optimal_bias_grid,
                 'Optimal Bias (uA)', str(ctime)+'_optbias')
    pd.tile_plot(len(rows), len(cols), crosstalk_bias_grid,
                 'Crosstalk Bias Limit (uA)', str(ctime)+'_crosstalk')


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
            ssa_params = calculate_ssa_parameters(
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

            sq1_params = calculate_sq1_parameters(rsdf, rs_runfile, cfg, col, row,
                                                  ssa_params, filter_sq1=filter_sq1, calc_slopes=calc_slopes,
                                                  sq1_sgfilter_window_length=sq1_sgfilter_window_length,
                                                  sq1_sgfilter_poly_deg=sq1_sgfilter_poly_deg)

            # ic_params=calculate_icminmax(cfg, filter_sq1, row, col,  sq1_params, ssa_params,
            #                             sq1_sgfilter_window_length, sq1_sgfilter_poly_deg)
            sq1_params2 = None
            if(rsdf_off is not None):
                sq1_params2 = calculate_sq1_parameters(rsdf_off, rs_runfile_off, cfg, col, row,
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
        numrows = max(sq1df['<row>'].astype(int))+1
        rows = range(0, numrows)
        sq1df_off = None
        sq1_runfile_off = None
        if(args.ctime_off is not None):
            sq1df_off, sq1_runfile_off = rd.get_sq1_tune_data(args.ctime_off)
        time1 = time.time()

        print('Done reading files, time elapsed (s):' + str(time1-time0))

        ic_driver(cfg, sa_data, sa_runfile, sq1df, sq1_runfile,
                  sq1df_off=sq1df_off,  sq1_runfile_off=sq1_runfile_off,
                  filter_sq1=True, ctime=ctime,
                  cols=cols, rows=rows)
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
    main()
