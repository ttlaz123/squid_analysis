import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import glob
import pandas as pd
import moby2
from moby2.util.mce import MCEFile
from moby2.util.mce import MCERunfile
import configparser
import scipy.constants
from scipy.signal import savgol_filter

# TeX plotting
from matplotlib import rc
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)
rc('text.latex', preamble=r'\usepackage{cmbright}')


def serial_corr(wave, lag=1):
    n = len(wave)
    y1 = wave[lag:]
    y2 = wave[:n-lag]
    corr = np.corrcoef(y1, y2)[0, 1]
    return corr


def autocorr(wave):
    lags = range(len(wave)//2)
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

    if debug:
        plt.figure()
        plt.plot(lags, corrs)
        plt.plot(peaks, corrs[peaks], "x")
        plt.ylabel('Auto-correlation')
        plt.xlabel('Lag')

    sorted_peaks = sorted([pk for _, pk in zip(corrs[peaks], peaks)])

    if debug:
        print(f'-> sorted_peaks[:4]={sorted_peaks[:4]}')

    try:
        phi0_idx = next(pk for pk in sorted_peaks if pk >
                        min_acorr_dist_from_zero)
    except:
        return None

    phi0 = np.abs(phi[phi0_idx]-phi[0])

    if debug:
        print(f'-> phi0={phi0}')

    return phi0


def get_midpoints(x):
    # X should be a numpy.ndarray
    return (x[1:] + x[:-1])/2.


plt.ion()

cfg_file = 'tune_cfg/slac_cd19.cfg'
cfg = configparser.ConfigParser()
cfg.read(cfg_file)

data_cryo = '/mnt/mce/'

# /mnt/mce/20220105/1641449102

filter_sq1 = True
sq1_sgfilter_window_length = 5
sq1_sgfilter_poly_deg = 2
s1_slope_npts_to_fit = 6

tune_ctime = sys.argv[1]
col = int(sys.argv[2])
row = int(sys.argv[3])
tune_dir = glob.glob(f'{data_cryo}/*/{tune_ctime}')[0]

# Tune analysis summary file
tasf = open(f'tas_{tune_ctime}.dat', 'a+')

##
# SSA
##

# Load SSA tune data
sa_tune = glob.glob(f'{tune_dir}/*_ssa')[0]
print(sa_tune)

# Open SSA tune data
sa_data = os.path.join(sa_tune)
mcefile = MCEFile(sa_data)
sa_data = mcefile.Read(field='error', row_col=True)

# Open SSA tune runfile
sa_runfile = MCERunfile(sa_tune+'.run')

sa_offset = sa_runfile.Item('HEADER', f'RB sa offset', type='int')
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


# sa_data.data[row][col] and each row is a copy of the same data.
nrow, ncol, sa_npt = sa_data.data.shape

# col = sa_biased_columns[0]

# Average all the SA curves on this column for better S/N
sa_coadd_adu = np.average([sa_data.data[row][col]
                           for row in range(nrow)], axis=0)

# Co-add is confusing, convert to ADU using sample_num
rc = int(col/8)+1
sample_num = sa_runfile.Item('HEADER', f'RB rc{rc} sample_num', type='int')[0]

sa_adu = sa_coadd_adu/sample_num

sa_bias_dac = sa_bias[col]
sa_bias_ua = 1.e6*sa_bias_dac*float(cfg['SSA']['SSA_BIAS_DAC_MAX_VOLTAGE_VOLTS'])/(
    (2**(float(cfg['SSA']['SSA_BIAS_DAC_NBITS'])))*(float(cfg['SSA']['SSA_BIAS_RC_RESISTANCE_OHMS'])+float(cfg['CRYOCABLE']['CRYOCABLE_ROUNDTRIP_RESISTANCE_OHMS'])))

# Convert sa_fb and sa_adu to physical units.

# Estimate the SSA phi0 from autocorrelation
sa_phi0_est = estimate_phi0(sa_fb, sa_adu, debug=False)

# Find the maximum in the first phi0
sa_max_idx = np.argmax(sa_adu[np.where(sa_fb < (np.min(sa_fb)+sa_phi0_est))])

# Find max and min slope in first full phi0
sa_fb_midpoints = get_midpoints(sa_fb)
sa_adu_midpoints = get_midpoints(sa_adu)
sa_adu_diff = np.diff(sa_adu)
# Max upslope
sa_adu_max_upslope_idx = np.argmax(sa_adu_diff[np.where(sa_fb < sa_phi0_est)])
# Max downslope
sa_adu_max_downslope_idx = np.argmin(
    sa_adu_diff[np.where(sa_fb < sa_phi0_est)])

# Plot
sa_fig, sa_ax = plt.subplots(figsize=(8, 6))
sa_ax.plot(sa_fb, sa_adu)

# Delineate phi0 boundaries
axes = plt.gca()
sa_y_min, sa_y_max = axes.get_ylim()
for phi in np.arange(0, np.max(sa_fb), sa_phi0_est) + sa_fb[sa_max_idx]:
    sa_ax.plot([phi, phi], [sa_y_min, sa_y_max], ls='--',
               color='k', lw=1, alpha=0.5)

# Show plot

# Why is this necessary?
sa_ax.set_ylim([sa_y_min, sa_y_max])
sa_ax.set_xlim([np.min(sa_fb), np.max(sa_fb)])

##
#
# Show max sensitivity points and slopes
# Upslope is green, downslope is green, MCE is magenta
#

# Upslope
sa_upslope_color = 'r'
plt.scatter(sa_fb_midpoints[sa_adu_max_upslope_idx],
            sa_adu_midpoints[sa_adu_max_upslope_idx],
            c=sa_upslope_color)
m_sa_adu_upslope = sa_adu_diff[sa_adu_max_upslope_idx] / \
    (d_sa_fb*sa_fb_dac_to_uA)
b_sa_adu_upslope = (sa_adu_midpoints[sa_adu_max_upslope_idx] -
                    m_sa_adu_upslope*sa_fb_midpoints[sa_adu_max_upslope_idx])

# Plot max upslope trend line
plt.plot([(sa_y_min-b_sa_adu_upslope)/m_sa_adu_upslope,
          (sa_y_max-b_sa_adu_upslope)/m_sa_adu_upslope],
         [sa_y_min, sa_y_max],
         sa_upslope_color+'--')

# Downslope
sa_downslope_color = 'g'
plt.scatter(sa_fb_midpoints[sa_adu_max_downslope_idx],
            sa_adu_midpoints[sa_adu_max_downslope_idx],
            c=sa_downslope_color)
m_sa_adu_downslope = sa_adu_diff[sa_adu_max_downslope_idx] / \
    (d_sa_fb*sa_fb_dac_to_uA)
b_sa_adu_downslope = (sa_adu_midpoints[sa_adu_max_downslope_idx] -
                      m_sa_adu_downslope*sa_fb_midpoints[sa_adu_max_downslope_idx])

# Plot max downslope trend line
plt.plot([(sa_y_min-b_sa_adu_downslope)/m_sa_adu_downslope,
          (sa_y_max-b_sa_adu_downslope)/m_sa_adu_downslope],
         [sa_y_min, sa_y_max],
         sa_downslope_color+'--')

###

###
#
# Anotate
# plt.tight_layout()
# sa_fig.subplots_adjust(right=0.1)
# plt.figtext(0.99, 0.75, ('TEST\n' +
#                         'TEST2\n')
#            )

dVdADU = float(cfg['PREAMPADC']['ADU_TO_VOLTS_AT_PREAMP_INPUT'])

Vmod_adu = np.max(sa_adu)-np.min(sa_adu)
Vmod_mV = 1.e3*Vmod_adu*dVdADU
M_ssa_fb_pH = 1.e12 * \
    scipy.constants.value(u'mag. flux quantum')/(sa_phi0_est*1.e-6)

# Downslopes
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

sa_ax.text(1.05, 0.5, (r"\textbf{Column "+str(col)+"}\n\n" +
                       r'\underline{Measured}' + '\n\n' +
                       r"$I^{SSA}_{bias}$ = "+f"{sa_bias_ua:.1f} $\mu$A\n" +
                       "$V^{SSA}_{mod}$ = " + f"{Vmod_mV:0.3f} mV\n" +
                       "$\Phi^{SSA\,FB}_{0}$ = "+f"{sa_phi0_est:.1f} $\mu$A\n" +
                       r"$M^{SSA}_{FB}$ = " + f"{M_ssa_fb_pH:0.1f} pH" + '\n\n' +
                       r"$dV^{SSA}_{ADU}/dFB^{SSA}_{DAC}$ $\downarrow$ = " + f'{dV_ADU_dSAFB_DAC_downslope:.2f}' + '\n' +
                       r"$dV^{SSA}_{nV}/dFB^{SSA}_{\mu\Phi_0}$ $\downarrow$ = " + f'{dV_nV_dSAFB_uphi0_downslope:.2f}' + '\n' +
                       r"$\Big|dI^{SSA}_{IN,pA}/dV^{SSA}_{nV}\Big|$ $\downarrow$ = " +
                       f'{np.abs(1./dV_nV_dSAIN_pA_downslope):.2f}' + '\n'
                       '\n' +
                       r"$dV^{SSA}_{ADU}/dFB^{SSA}_{DAC}$ $\uparrow$ = " + f'{dV_ADU_dSAFB_DAC_upslope:.2f}' + '\n' +
                       r"$dV^{SSA}_{nV}/dFB^{SSA}_{\mu\Phi_0}$ $\uparrow$ = " + f'{dV_nV_dSAFB_uphi0_upslope:.2f}' + '\n' +
                       r"$\Big|dI^{SSA}_{IN,pA}/dV^{SSA}_{nV}\Big|$ $\uparrow$ = " +
                       f'{np.abs(1./dV_nV_dSAIN_pA_upslope):.2f}' + '\n'
                       "\n" +
                       r'\underline{Assumed}' + '\n\n' +
                       f"nV/ADU = {1.e9*float(cfg['PREAMPADC']['ADU_TO_VOLTS_AT_PREAMP_INPUT']):.1f}" + '\n' +
                       r"$M^{SSA}_{IN}$ = " +
                       f"{cfg['SSA']['SSA_M_IN_PICOHENRY']} pH" + '\n' +
                       f'SSAFB nA/DAC = {1000.*sa_fb_dac_to_uA:.3f}' + '\n' +
                       "$R^{cryo\,cable}_{roundtrip}$ = " +
                       f"{cfg['CRYOCABLE']['CRYOCABLE_ROUNDTRIP_RESISTANCE_OHMS']} $\Omega$"
                       ),
           size=10, weight='bold', linespacing=1.5,
           bbox=dict(edgecolor='k', facecolor='none', pad=10, linewidth=1),
           ha='left',
           va='center',
           transform=sa_ax.transAxes)

plt.subplots_adjust(right=0.75)

plt.xlabel('SSA Feedback Current ($\mu$A)', fontsize=18)
plt.ylabel('SSA Voltage (ADU)', fontsize=18)

plt.savefig(f'{tune_ctime}_sa_c{col}r{row}.png')

# DONE PLOTTING SSA

# START PLOTTING SQ1

# Load SQ1 tune data
sq1_tune = glob.glob(f'{tune_dir}/*_sq1servo_sa.bias')[0]
print(sq1_tune)

# Open SQ1 tune data
sq1_data = os.path.join(sq1_tune)
# sq1_data = mcefile.Read(field='error', row_col=True)
sq1df = pd.read_csv(sq1_data, delim_whitespace=True,
                    error_bad_lines=False, index_col=False)
# In [12]: sq1df.columns
# Out[12]:
# Index(['<bias>', '<flux>', '<row>', '<error00>', '<error01>', '<error02>',
#       '<error03>', '<error04>', '<error05>', '<error06>', '<error07>',
#       '<safb00>', '<safb01>', '<safb02>', '<safb03>', '<safb04>', '<safb05>',
#       '<safb06>', '<safb07>'],
#      dtype='object')

# For requested row, pull the bias


# Open SQ1 tune runfile
sq1_runfile = MCERunfile(sq1_tune.replace('.bias', '.run'))

# What range of SQ1 biases did we sweep?
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
for idx, b in enumerate(sq1_b):
    sq1_safb_servo = sq1df[(sq1df['<row>'] == row) & (
        sq1df['<bias>'] == idx)][f'<safb{col:02}>'].values
    sq1_safb_servo_curves_dac.append(sq1_safb_servo)
    sq1_safb_servo_biases_dac.append(b)
    sq1_safb_servo_span = np.max(sq1_safb_servo)-np.min(sq1_safb_servo)
    if sq1_safb_servo_span > max_sq1_safb_servo_span:
        max_sq1_safb_servo_span = sq1_safb_servo_span
        max_sq1_safb_servo_span_bias = b

print(
    f'max_sq1_safb_servo_span_bias={max_sq1_safb_servo_span_bias}\tmax_sq1_safb_servo_span={max_sq1_safb_servo_span}')

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

###########################
# Plot
s1_fig, s1_ax = plt.subplots(figsize=(8, 6))
sq1_fb_uA = sq1fb_dac_to_nA*sq1_fb/1000.
max_sq1_safb_servo_uA_unfilt = max_sq1_safb_servo*sa_fb_dac_to_uA

max_sq1_safb_servo_uA = max_sq1_safb_servo_uA_unfilt

if filter_sq1:
    max_sq1_safb_servo_uA = savgol_filter(
        max_sq1_safb_servo_uA, sq1_sgfilter_window_length, sq1_sgfilter_poly_deg)

plt.plot(sq1_fb_uA, max_sq1_safb_servo_uA)
plt.scatter(sq1_fb_uA, max_sq1_safb_servo_uA, s=10)

if filter_sq1:
    plt.plot(sq1_fb_uA, max_sq1_safb_servo_uA_unfilt,
             linestyle='--', color='gray')

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

# Delineate phi0 boundaries
axes = plt.gca()
s1_y_min, s1_y_max = axes.get_ylim()
for phi in np.arange(0, np.max(sq1_fb_uA), s1_phi0_est) + sq1_fb_uA[s1_max_idx]:
    s1_ax.plot([phi, phi], [s1_y_min, s1_y_max], ls='--',
               color='k', lw=1, alpha=0.5)

# Why is this necessary?
s1_ax.set_ylim([s1_y_min, s1_y_max])
s1_ax.set_xlim([np.min(sq1_fb_uA), np.max(sq1_fb_uA)])

# Downslope
s1_downslope_color = 'g'
plt.scatter(sq1_fb_uA_midpoints[max_sq1_safb_servo_uA_max_downslope_idx],
            max_sq1_safb_servo_uA_midpoints[max_sq1_safb_servo_uA_max_downslope_idx],
            c=s1_downslope_color)
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

    plt.scatter(s1_fb_uA_subset_for_downslope_fit, max_sq1_safb_servo_uA_subset_for_downslope_fit,
                s=50, facecolors='None', edgecolors=s1_downslope_color)

    # Plot max downslope trend line
    plt.plot([(s1_y_min-b_max_sq1_safb_servo_uA_downslope)/m_max_sq1_safb_servo_uA_downslope,
              (s1_y_max-b_max_sq1_safb_servo_uA_downslope)/m_max_sq1_safb_servo_uA_downslope],
             [s1_y_min, s1_y_max],
             s1_downslope_color+'--')
except:
    print('!! Failed to fit SQ1 downslope!')

# Upslope
s1_upslope_color = 'r'
plt.scatter(sq1_fb_uA_midpoints[max_sq1_safb_servo_uA_max_upslope_idx],
            max_sq1_safb_servo_uA_midpoints[max_sq1_safb_servo_uA_max_upslope_idx],
            c=s1_upslope_color)

# SQ1 curves tend to be a little noisier.  Fit neighboring points
s1_fb_uA_subset_for_upslope_fit = sq1_fb_uA[max_sq1_safb_servo_uA_max_upslope_idx -
                                            int(s1_slope_npts_to_fit/2-1):max_sq1_safb_servo_uA_max_upslope_idx+int(s1_slope_npts_to_fit/2+1)]
max_sq1_safb_servo_uA_subset_for_upslope_fit = max_sq1_safb_servo_uA[max_sq1_safb_servo_uA_max_upslope_idx - int(
    s1_slope_npts_to_fit/2-1):max_sq1_safb_servo_uA_max_upslope_idx+int(s1_slope_npts_to_fit/2+1)]
s1_upslope_p1 = np.polyfit(
    s1_fb_uA_subset_for_upslope_fit, max_sq1_safb_servo_uA_subset_for_upslope_fit, 1)
m_max_sq1_safb_servo_uA_upslope = s1_upslope_p1[0]
b_max_sq1_safb_servo_uA_upslope = s1_upslope_p1[1]

plt.scatter(s1_fb_uA_subset_for_upslope_fit, max_sq1_safb_servo_uA_subset_for_upslope_fit,
            s=50, facecolors='None', edgecolors=s1_upslope_color)

# Plot max upslope trend line
plt.plot([(s1_y_min-b_max_sq1_safb_servo_uA_upslope)/m_max_sq1_safb_servo_uA_upslope,
          (s1_y_max-b_max_sq1_safb_servo_uA_upslope)/m_max_sq1_safb_servo_uA_upslope],
         [s1_y_min, s1_y_max],
         s1_upslope_color+'--')

plt.ylabel('SSA Feedback Current ($\mu$A)', fontsize=18)
plt.xlabel('SQ1 Feedback Current ($\mu$A)', fontsize=18)

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

s1_ax.text(1.05, 0.5, (r"\textbf{SQ1 on " + f'c{col:02}r{row:02}' + "}\n\n" +
                       r'\underline{Measured}' + '\n\n' +
                       "$\Phi^{SQ1\,FB}_{0}$ = " +
                       f"{s1_phi0_est:.1f} $\mu$A\n" + '\n'
                       r"$M^{SQ1}_{FB}$ = " + f"{M_s1_fb:0.1f} pH" + '\n' +
                       '\n' +
                       r"$dI^{SSA}_{FB,\mu A}/dI^{SQ1}_{FB,\mu A}$ $\downarrow$ = " + (f'{m_max_sq1_safb_servo_uA_downslope:.2f}' if m_max_sq1_safb_servo_uA_downslope is not None else 'None') + '\n' +
                       r"$dI^{SSA}_{IN,\mu A}/dI^{SQ1}_{IN,\mu A}$ $\downarrow$ = " + (f'{dI_SSA_IN_pA_dI_SQ1_IN_pA_downslope:.2f}' if dI_SSA_IN_pA_dI_SQ1_IN_pA_downslope is not None else 'None') + '\n' +
                       '\n' +
                       r"$dI^{SSA}_{FB,\mu A}/dI^{SQ1}_{FB,\mu A}$ $\uparrow$ = " + f'{m_max_sq1_safb_servo_uA_upslope:.2f}' + '\n' +
                       r"$dI^{SSA}_{IN,\mu A}/dI^{SQ1}_{IN,\mu A}$ $\uparrow$ = " + f'{dI_SSA_IN_pA_dI_SQ1_IN_pA_upslope:.2f}' + '\n' +
                       "\n" +
                       "\n" +
                       r'\underline{Assumed}' + '\n\n' +
                       f'SSAFB nA/DAC = {1000.*sa_fb_dac_to_uA:.3f}' + '\n' +
                       "$R^{cryo\,cable}_{roundtrip}$ = " +
                       f"{cfg['CRYOCABLE']['CRYOCABLE_ROUNDTRIP_RESISTANCE_OHMS']} $\Omega$" + '\n' +
                       r"$M^{SQ1}_{IN}$ = " +
                       f"{cfg['SQ1']['SQ1_M_IN_PICOHENRY']} pH" + '\n'
                       ),
           size=10, weight='bold', linespacing=1.5,
           bbox=dict(edgecolor='k', facecolor='none', pad=10, linewidth=1),
           ha='left',
           va='center',
           transform=sa_ax.transAxes)

plt.tight_layout()

plt.subplots_adjust(right=0.75)

plt.savefig(f'{tune_ctime}_sq1_c{col}r{row}.png')


# just take first row since should be all the same
# adc_offset = runfile.Item(
#    'HEADER', f'RB rc{rc} adc_offset{col%8}', type='int')[0]

###########################
# Plot SQ1 min & max vs SQ1 bias
s1b_minmax_fig, s1b_minmax_ax = plt.subplots(figsize=(8, 6))
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

#sa_fb_servo_maxs_uA =  sa_fb_dac_to_uA*(
#    np.mean(sq1_safb_servo_curves_dac[0])-sq1_safb_servo_maxs_dac) #The max and min curve should be identical where it is superconducting, but the max curve is superconducting longer.

#min_sc_branch_current_uA = 1
#max_sc_branch_current_uA = 5.5
#indices_sc = np.where(
#    (sq1_safb_servo_biases_uA > min_sc_branch_current_uA) &
#    (sq1_safb_servo_biases_uA < max_sc_branch_current_uA))
#zsc = np.polyfit(sq1_safb_servo_biases_uA[indices_sc], sa_fb_servo_maxs_uA[indices_sc], 1)

#M_ssa_in_pH = M_ssa_fb_pH/zsc[0]
#print(M_ssa_in_pH)
#print(float(cfg['SSA']['SSA_M_IN_PICOHENRY']))

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

plt.plot(sq1_safb_servo_biases_uA, sq1_safb_servo_mins_sa_in_uA,
         lw=2, label='SQ1 min(Imod)')
plt.plot(sq1_safb_servo_biases_uA, sq1_safb_servo_maxs_sa_in_uA,
         lw=2, label='SQ1 max(Imod)')

# plot max mod
sq1_safb_servo_span_sa_in_uA = (sq1_safb_servo_maxs_sa_in_uA -
                                sq1_safb_servo_mins_sa_in_uA)
max_sq1imod_idx = np.argmax(sq1_safb_servo_span_sa_in_uA)
max_sq1imod_uA = (sq1_safb_servo_span_sa_in_uA)[max_sq1imod_idx]

plt.plot([sq1_safb_servo_biases_uA[max_sq1imod_idx], sq1_safb_servo_biases_uA[max_sq1imod_idx]],
         [sq1_safb_servo_mins_sa_in_uA[max_sq1imod_idx],
          sq1_safb_servo_maxs_sa_in_uA[max_sq1imod_idx]],
         'm', lw=3,
         label='$I^{SQ1}_{mod}$ = '+f'{max_sq1imod_uA:.3f} $\mu$A @ '+'$I_{SQ1B,total} = $'+f'{sq1_safb_servo_biases_uA[max_sq1imod_idx]:.1f} $\mu$A')

plt.plot(sq1_safb_servo_biases_uA, sq1_safb_servo_biases_uA)

tasf.write(
    f'{col}\t{row}\t{max_sq1_safb_servo_span_bias}\t{sq1_safb_servo_biases_uA[max_sq1imod_idx]: .1f}\t{max_sq1_safb_servo_span}\t{max_sq1imod_uA:.3f}\n')

plt.legend(loc='upper left', fontsize=8)

plt.ylabel('SSA Input Current ($\mu$A)', fontsize=18)
plt.xlabel('SQ1 Total Bias Current ($\mu$A)', fontsize=18)

plt.ylim(0,40)

plt.tight_layout()
plt.savefig(f'{tune_ctime}_sq1minmax_c{col}r{row}.png')

# subtract zero offsets


###########################
# Plot SQ1 max-min vs SQ1 bias
s1b_span_fig, s1b_span_ax = plt.subplots(figsize=(8, 6))

# Close tune analysis summary file
tasf.close()


###########################
# Calculate R_dyn
step_size = 7
sq1_current_uA = np.array(sq1_safb_servo_curves_dac[0] - sq1_safb_servo_curves_dac)*sa_fb_dac_to_uA * (M_ssa_fb_pH/float(cfg['SSA']['SSA_M_IN_PICOHENRY']))
sq1_current_diff_uA = sq1_current_uA[0:-step_size, :]-sq1_current_uA[step_size-1:-1, :]
sq1_current_min_matrix = np.transpose(np.tile(np.min(sq1_current_uA, axis=1), (np.shape(sq1_current_uA)[1], 1)))
sq1_current_max_matrix = np.transpose(np.tile(np.max(sq1_current_uA, axis=1), (np.shape(sq1_current_uA)[1], 1)))
sq1_current_frac = (sq1_current_uA-sq1_current_min_matrix)/(sq1_current_max_matrix - sq1_current_min_matrix)

sq1_device_uV = (sq1_safb_servo_biases_uA[:, None]-sq1_current_uA)* float(cfg['SQ1']['R_SHUNT'])
sq1_device_diff_uV = sq1_device_uV[:-step_size] - sq1_device_uV[step_size:]
#sq1_bias_diff_uA = sq1_safb_servo_biases_uA[0:-step_size] - sq1_safb_servo_biases_uA[step_size-1:-1]
#sq1_bias_diff_uV = np.transpose(np.tile(sq1_bias_diff_uA* float(cfg['SQ1']['R_SHUNT']), (np.shape(sq1_current_uA)[1], 1)))
r_dyn = sq1_device_diff_uV/sq1_current_diff_uA

###########################
# Plot R_dyn at fixed feedback
r_dyn_fig, r_dyn_ax = plt.subplots(figsize=(8, 6))

# Find the minimum in the first phi0
s1_min_idx = np.argmin(
    max_sq1_safb_servo_uA[np.where(sq1_fb_uA < (np.min(sq1_fb_uA)+s1_phi0_est))])

fb_indices = [s1_max_idx, s1_min_idx]
for idx in fb_indices:
    plt.plot(sq1_current_uA[:-step_size, idx], r_dyn[:, idx])
plt.ylabel('Dynamic SQ1 Resistance ($\Omega$)', fontsize=18)
plt.xlabel('SQ1 Device Current ($\mu$A)', fontsize=18)
plt.ylim(0, 30)
plt.xlim(0, 20)

plt.tight_layout()
plt.savefig(f'{tune_ctime}_r_dyn_vs_i_dev_c{col}r{row}.png')

###########################
# Plot R_dyn at fixed Bias
r_dyn2_fig, r_dyn2_ax = plt.subplots(figsize=(8, 6))
plt.plot(sq1_fb_uA, r_dyn[max_sq1imod_idx, :])

plt.ylabel('Dynamic SQ1 Resistance ($\Omega$)', fontsize=18)
plt.xlabel('SQ1 Feedback Current ($\mu$A)', fontsize=18)
plt.ylim(0, 25)
plt.tight_layout()
plt.savefig(f'{tune_ctime}_r_dyn_vs_fb_c{col}r{row}.png')


r_dyn3_fig, r_dyn3_ax = plt.subplots(figsize=(8, 6))
plt.plot(sq1_current_frac[max_sq1imod_idx], r_dyn[max_sq1imod_idx, :])
    
plt.ylabel('Dynamic SQ1 Resistance ($\Omega$)', fontsize=18)
plt.xlabel('Idev - Idevmin/Idevmax - Idevmin', fontsize=18)
plt.ylim(0, 25)
plt.tight_layout()
plt.savefig(f'{tune_ctime}_r_dyn_vs_ifrac_c{col}r{row}.png')


###########################
# Calculate R_N
#sq1_norm_biases_uV = np.transpose(np.tile(sq1_safb_servo_biases_uA* float(cfg['SQ1']['R_SHUNT']), (np.shape(sq1_current_uA)[1], 1)))
#r_n = sq1_norm_biases_uV[0:-step_size, :]*r_dyn/(sq1_current_uA[0:-step_size, :]*(r_dyn + float(cfg['SQ1']['R_SHUNT'])))
r_n = sq1_device_uV/sq1_current_uA#In case specific values are needed
r_n_min = (sq1_safb_servo_biases_uA - sq1_safb_servo_mins_sa_in_uA)* float(cfg['SQ1']['R_SHUNT'])/sq1_safb_servo_mins_sa_in_uA
r_n_max = (sq1_safb_servo_biases_uA - sq1_safb_servo_maxs_sa_in_uA)* float(cfg['SQ1']['R_SHUNT'])/sq1_safb_servo_maxs_sa_in_uA

sc_region = np.where((sq1_safb_servo_maxs_sa_in_uA > 5) & (sq1_safb_servo_maxs_sa_in_uA < 9))
r_series = np.mean(r_n[sc_region, idx])
print('R_series = ' + f'{r_series}')

r_n_fig, r_n_ax = plt.subplots(figsize=(8, 6))
plt.axhline(y=r_series)

plt.plot(sq1_safb_servo_mins_sa_in_uA, r_n_min)
plt.plot(sq1_safb_servo_maxs_sa_in_uA, r_n_max)

plt.ylabel('SQ1 $R_N$ ($\Omega$)', fontsize=18)
plt.xlabel('SQ1 Device Current ($\mu$A)', fontsize=18)
plt.ylim(0, 8)
plt.xlim(0, 20)
#plt.legend()
plt.grid(visible=True)

plt.tight_layout()
plt.savefig(f'{tune_ctime}_r_N_vs_i_dev_c{col}r{row}.png')

max_sq1_mod_uV = sq1_safb_servo_biases_uA[max_sq1imod_idx]* float(cfg['SQ1']['R_SHUNT'])#Without shunt, this is aka Ic_mod_max in NIST jargon
ic_mod_max = sq1_safb_servo_maxs_sa_in_uA[max_sq1imod_idx]
ic_idx = np.where(sq1_safb_servo_span_sa_in_uA > 0.07)[0][0]
ic_mod_min = np.nanmean(sq1_current_uA[ic_idx, :])
r_n_opt = r_n[max_sq1imod_idx, max_sq1_safb_servo_uA_max_upslope_idx]
print('I_SQ1 = '+ f'{sq1_current_uA[max_sq1imod_idx, max_sq1_safb_servo_uA_max_upslope_idx]}')
print('I_b = '+ f'{sq1_safb_servo_biases_uA[max_sq1_safb_servo_uA_max_upslope_idx]}')
print('R_N = '+ f'{r_n_opt}')
print('R_Dyn = ' + f'{r_dyn[max_sq1imod_idx, max_sq1_safb_servo_uA_max_upslope_idx]}')

print('Ic_mod_min = '+ f'{ic_mod_min}')
print('Ic_mod_max = '+ f'{ic_mod_max}')
print('dV = ' + f'{max_sq1_mod_uV}')
bL_sq1 = (4*ic_mod_max*r_n_opt/(np.pi * max_sq1_mod_uV)) - 1
print('b_L = '+ f'{bL_sq1}')
L_sq1_pH = scipy.constants.value(u'mag. flux quantum') * bL_sq1 * 1e12/(2*ic_mod_max*1e-6)
print('L = ' + f'{L_sq1_pH}')

## Debug step size
r_test_fig, r_test_ax = plt.subplots(figsize=(8, 6))
plt.plot(sq1_safb_servo_biases_uA, sq1_safb_servo_mins_sa_in_uA)
plt.plot(sq1_safb_servo_biases_uA, sq1_safb_servo_maxs_sa_in_uA)
plt.plot(sq1_safb_servo_biases_uA, sq1_safb_servo_biases_uA)
plt.xlim(0, 10)
plt.ylim(0, 5)
plt.savefig(f'{tune_ctime}_testing_c{col}r{row}.png')

step_sizes = range(1, 10)
r_test2_fig, r_test2_ax = plt.subplots(figsize=(8, 6))
for step_size in step_sizes:
    sq1_current_uA = np.array(sq1_safb_servo_curves_dac[0] - sq1_safb_servo_curves_dac)*sa_fb_dac_to_uA * (M_ssa_fb_pH/float(cfg['SSA']['SSA_M_IN_PICOHENRY']))
    sq1_current_diff_uA = sq1_current_uA[0:-step_size, :]-sq1_current_uA[step_size-1:-1, :]
    sq1_current_min_matrix = np.transpose(np.tile(np.min(sq1_current_uA, axis=1), (np.shape(sq1_current_uA)[1], 1)))
    sq1_current_max_matrix = np.transpose(np.tile(np.max(sq1_current_uA, axis=1), (np.shape(sq1_current_uA)[1], 1)))
    sq1_current_frac = (sq1_current_uA-sq1_current_min_matrix)/(sq1_current_max_matrix - sq1_current_min_matrix)

    sq1_bias_diff_uA = sq1_safb_servo_biases_uA[0:-step_size] - sq1_safb_servo_biases_uA[step_size-1:-1]
    sq1_bias_diff_uV = np.transpose(np.tile(sq1_bias_diff_uA* float(cfg['SQ1']['R_SHUNT']), (np.shape(sq1_current_uA)[1], 1)))
    r_dyn = sq1_bias_diff_uV/sq1_current_diff_uA
    plt.plot(sq1_fb_uA, sq1_current_diff_uA[max_sq1imod_idx, :])

#plt.savefig(f'{tune_ctime}_testing2_c{col}r{row}.png')

filename = f'{tune_ctime}_config_c{col}r{row}.txt'
f = open(filename, "w")

f.write("#Note all values given in SI units unless labeled otherwise\n")
f.write("\n[PREAMPADC]\n")
f.write(f"RS On Measurement time: {tune_ctime}\n")
#f.write("Chip ID: {}")
f.write(f"Row: {row}\n")
f.write("ADU_TO_VOLTS_AT_PREAMP_INPUT=7.075e-7\n")

f.write("\n[CRYOCABLE]\n")
f.write(f"CRYOCABLE_ROUNDTRIP_RESISTANCE_OHMS = {cfg['CRYOCABLE']['CRYOCABLE_ROUNDTRIP_RESISTANCE_OHMS']}\n")
f.write("#From cryostat config file\n")

f.write("\n[SSA]\n")
f.write(f"I_SSA_BIAS_OPT = {sa_bias_ua*1e-6:.2e}\n")
f.write(f"V_SSA_MOD_OPT = {Vmod_mV*1e-3:.2e}\n")
f.write(f"PHI_SSA_FB_0 = {sa_phi0_est*1e-6:.2e}\n")
f.write(f"M_SSA_FB = {M_ssa_fb_pH*1e-12:.2e}\n")
f.write(f"dV_SSA_dI_SSA_IN_DOWNSLOPE = {dV_nV_dSAIN_pA_downslope*1e3:.2e}\n")
f.write(f"dV_SSA_dI_SSA_IN_UPSLOPE = {dV_nV_dSAIN_pA_upslope*1e3:.2e}\n")
        
f.write("\n#Assumed from previous testing\n")
f.write(f"M_SSA_IN = {float(cfg['SSA']['SSA_M_IN_PICOHENRY'])*1e-12:.2e}\n")
f.write(f"SSA_FB_CAL = {1e-6*sa_fb_dac_to_uA:.2e}\n")
f.write("#MCE Calibration in Amps per DAC unit\n\n")

f.write("\n[SQ1]\n")
f.write(f"PHI_SQ1_FB_0 = {s1_phi0_est*1e-6:.2e}\n")
f.write(f"M_SQ1_FB = {M_s1_fb*1e-12:.2e}\n")
f.write(f"Ic_MOD_MIN = {ic_mod_min*1e-6:.2e}\n")
f.write(f"Ic_MOD_MAX = {ic_mod_max*1e-6:.2e}\n")
f.write(f"dI_SSA_IN_dI_SQ1_IN_downslope = {dI_SSA_IN_pA_dI_SQ1_IN_pA_downslope:.2e}\n")
f.write(f"dI_SSA_IN_dI_SQ1_IN_upslope = {dI_SSA_IN_pA_dI_SQ1_IN_pA_upslope:.2e}\n")
f.write(f"I_SQ1_OPERATING_UPSLOPE = {sq1_current_uA[max_sq1imod_idx, max_sq1_safb_servo_uA_max_upslope_idx]:.2e}\n")
f.write(f"I_BIAS_OPERATING_UPSLOPE = {sq1_safb_servo_biases_uA[max_sq1_safb_servo_uA_max_upslope_idx]:.2e}\n")
f.write(f"R_N_OPERATING_UPSLOPE = {r_n[max_sq1imod_idx, max_sq1_safb_servo_uA_max_upslope_idx]:.2e}\n")
f.write(f"R_DYN_OPERATING_UPSLOPE = {r_dyn[max_sq1imod_idx, max_sq1_safb_servo_uA_max_upslope_idx]:.2e}\n")
f.write(f"I_SQ1_OPERATING_DOWNSLOPE = {sq1_current_uA[max_sq1imod_idx, max_sq1_safb_servo_uA_max_downslope_idx]:.2e}\n")
f.write(f"I_BIAS_OPERATING_DOWNSLOPE = {sq1_safb_servo_biases_uA[max_sq1_safb_servo_uA_max_downslope_idx]:.2e}\n")
f.write(f"R_N_OPERATING_DOWNSLOPE = {r_n[max_sq1imod_idx, max_sq1_safb_servo_uA_max_downslope_idx]:.2e}\n")
f.write(f"R_DYN_OPERATING_DOWNSLOPE = {r_dyn[max_sq1imod_idx, max_sq1_safb_servo_uA_max_downslope_idx]:.2e}\n")
f.write(f"V_MODULATION = {max_sq1_mod_uV:.2e}\n")
f.write(f"b_L = {bL_sq1:.2e}\n")
f.write(f"L_SQ1 = {L_sq1_pH*1e-12:.2e}\n")
f.write(f"R_SERIES = {r_series}\n")

f.write("\n#Assumed from previous testing\n")
f.write(f"M_SQ1_IN = {float(cfg['SQ1']['SQ1_M_IN_PICOHENRY'])*1e-12:.2e}\n")

f.close()
