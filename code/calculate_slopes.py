
import numpy as np  
import scipy 
from scipy.signal import savgol_filter
from estimate_phi0 import estimate_phi0


###
# DEPRECATED: NO LONGER USED BY SQUID TUNING ANALYSIS PLOTS
###
def get_midpoints(x):
    # X should be a numpy.ndarray
    return (x[1:] + x[:-1])/2.


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

def calc_ssa_slopes(sa_adu, sa_fb, sa_phi0_est, cfg, d_sa_fb, sa_fb_dac_to_uA,
                    col, sa_bias_ua, M_ssa_fb_pH, show_plot):
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
    return sa_ax
