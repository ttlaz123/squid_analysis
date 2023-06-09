'''
Written by Tom Liu, Documentation last updated 2023 May 31
Handles all the plotting done in Squid tuning analysis
'''

from matplotlib.ticker import AutoMinorLocator
from collections import OrderedDict
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
# Reduces memory leakage or something idk
matplotlib.use('Agg')


def plot_ssa_debug(sa_fb, sa_adu):
    """
    DEPRECATED
    """
    sa_fig, sa_ax = plt.subplots(figsize=(8, 6))
    sa_ax.plot(sa_fb, sa_adu)
    # plt.show()
    plt.close()


def plot_ssa(cfg, sa_fb, sa_adu, sa_phi0_est, sa_max_idx,
             sa_fb_midpoints, sa_adu_midpoints,
             sa_adu_max_upslope_idx, sa_adu_max_downslope_idx,
             sa_fb_dac_to_uA, b_sa_adu_upslope, m_sa_adu_upslope,
             b_sa_adu_downslope, m_sa_adu_downslope,
             dV_ADU_dSAFB_DAC_downslope, dV_nV_dSAFB_uphi0_downslope, dV_nV_dSAIN_pA_downslope,
             dV_ADU_dSAFB_DAC_upslope, dV_nV_dSAFB_uphi0_upslope, dV_nV_dSAIN_pA_upslope,
             col, sa_bias_ua, Vmod_mV, M_ssa_fb_pH, show_plot=True):
    '''
    DEPRECATED
    Plots SSA tuning results along with the optimal slopes
    '''
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

    # Plot max downslope trend line
    plt.plot([(sa_y_min-b_sa_adu_downslope)/m_sa_adu_downslope,
              (sa_y_max-b_sa_adu_downslope)/m_sa_adu_downslope],
             [sa_y_min, sa_y_max],
             sa_downslope_color+'--')

    sa_ax.text(1.05, 0.5, (r"\textbf{Column "+str(col)+"}\n\n" +
                           r'\underline{Measured}' + '\n\n' +
                           r"$I^{SSA}_{bias}$ = "+f"{sa_bias_ua:.1f} $\mu$A\n" +
                           "$V^{SSA}_{mod}$ = " + f"{Vmod_mV:0.3f} mV\n" +
                           "$\Phi^{SSA\,FB}_{0}$ = "+f"{sa_phi0_est:.1f} $\mu$A\n" +
                           r"$M^{SSA}_{FB}$ = " + f"{M_ssa_fb_pH:0.1f} pH" + '\n\n' +
                           r"$dV^{SSA}_{ADU}/dFB^{SSA}_{DAC}$ $\downarrow$ = " + f'{dV_ADU_dSAFB_DAC_downslope:.2f}' + '\n' +
                           r"$dV^{SSA}_{nV}/dFB^{SSA}_{\mu\Phi_0}$ $\downarrow$ = " + f'{dV_nV_dSAFB_uphi0_downslope:.2f}' + '\n' +
                           r"$|dI^{SSA}_{IN,pA}/dV^{SSA}_{nV}|$ $\downarrow$ = " +
                           f'{np.abs(1./dV_nV_dSAIN_pA_downslope):.2f}' + '\n'
                           '\n' +
                           r"$dV^{SSA}_{ADU}/dFB^{SSA}_{DAC}$ $\uparrow$ = " + f'{dV_ADU_dSAFB_DAC_upslope:.2f}' + '\n' +
                           r"$dV^{SSA}_{nV}/dFB^{SSA}_{\mu\Phi_0}$ $\uparrow$ = " + f'{dV_nV_dSAFB_uphi0_upslope:.2f}' + '\n' +
                           r"$|dI^{SSA}_{IN,pA}/dV^{SSA}_{nV}|$ $\uparrow$ = " +
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
    if(show_plot):
        plt.show()
    # plt.savefig(f'{tune_ctime}_sa_c{col}r{row}.png')
    plt.close()
    return sa_ax


def plot_sq1(col, row, sq1_fb_uA, max_sq1_safb_servo_uA, filter_sq1,
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
             sa_fb_dac_to_uA, cfg, show_plot=False):
    '''
    DEPRECATED
    Plots sq1 tuning curves with optimal slopes
    '''
    s1_fig, s1_ax = plt.subplots(figsize=(8, 6))

    plt.plot(sq1_fb_uA, max_sq1_safb_servo_uA)
    plt.scatter(sq1_fb_uA, max_sq1_safb_servo_uA, s=10)

    if filter_sq1:
        plt.plot(sq1_fb_uA, max_sq1_safb_servo_uA_unfilt,
                 linestyle='--', color='gray')

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

    plt.scatter(s1_fb_uA_subset_for_downslope_fit, max_sq1_safb_servo_uA_subset_for_downslope_fit,
                s=50, facecolors='None', edgecolors=s1_downslope_color)

    # Plot max downslope trend line
    plt.plot([(s1_y_min-b_max_sq1_safb_servo_uA_downslope)/m_max_sq1_safb_servo_uA_downslope,
              (s1_y_max-b_max_sq1_safb_servo_uA_downslope)/m_max_sq1_safb_servo_uA_downslope],
             [s1_y_min, s1_y_max],
             s1_downslope_color+'--')

    # Upslope
    s1_upslope_color = 'r'
    plt.scatter(sq1_fb_uA_midpoints[max_sq1_safb_servo_uA_max_upslope_idx],
                max_sq1_safb_servo_uA_midpoints[max_sq1_safb_servo_uA_max_upslope_idx],
                c=s1_upslope_color)

    plt.scatter(s1_fb_uA_subset_for_upslope_fit, max_sq1_safb_servo_uA_subset_for_upslope_fit,
                s=50, facecolors='None', edgecolors=s1_upslope_color)

    # Plot max upslope trend line
    plt.plot([(s1_y_min-b_max_sq1_safb_servo_uA_upslope)/m_max_sq1_safb_servo_uA_upslope,
              (s1_y_max-b_max_sq1_safb_servo_uA_upslope)/m_max_sq1_safb_servo_uA_upslope],
             [s1_y_min, s1_y_max],
             s1_upslope_color+'--')

    plt.ylabel('SSA Feedback Current ($\mu$A)', fontsize=18)
    plt.xlabel('SQ1 Feedback Current ($\mu$A)', fontsize=18)

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
               transform=s1_ax.transAxes)

    plt.tight_layout()

    plt.subplots_adjust(right=0.75)
    if(show_plot):
        plt.show()
    # plt.savefig(f'{tune_ctime}_sq1_c{col}r{row}.png')
    plt.close()


def plot_iccurve(ic_params, ax, alpha=1,
                 min_label=None, min_color=None,
                 max_label=None, max_color=None,
                 mod_label=None, mod_color=None,
                 thresh_label=None, thresh_color=None,
                 picked_label=None, picked_color=None, picked_bias_idx=None):
    """
    Plots an I-V curve with specified parameters and annotations.

    Args:
        ic_params (dict): Dictionary containing the following keys:
            - 'bias': List or array of bias currents.
            - 'fb_min': List or array of minimum feedback currents.
            - 'fb_max': List or array of maximum feedback currents.
            - 'bias_max_idx': Index of the maximum bias current.
            - 'bias_min_idx': Index of the bias current at which modulation begins.

        ax (AxesSubplot): The matplotlib axes on which to plot the I-V curve.

        alpha (float, optional): Transparency value for the plotted curves. Defaults to 1.

        min_label (str, optional): Label for the minimum feedback current curve. Defaults to None.

        min_color (str or tuple, optional): Color for the minimum feedback current curve. Defaults to None.

        max_label (str, optional): Label for the maximum feedback current curve. Defaults to None.

        max_color (str or tuple, optional): Color for the maximum feedback current curve. Defaults to None.

        mod_label (str, optional): Label for the modified feedback current curve. Defaults to None.

        mod_color (str or tuple, optional): Color for the modified feedback current curve. Defaults to None.

        thresh_label (str, optional): Label for the threshold lines. Defaults to None.

        thresh_color (str or tuple, optional): Color for the threshold lines. Defaults to None.

        picked_label (str, optional): Label for the picked bias current annotation. Defaults to None.

        picked_color (str or tuple, optional): Color for the picked bias current annotation. Defaults to None.

        picked_bias_idx (int, optional): Index of the picked bias current. Defaults to None.

    Returns:
        AxesSubplot: The matplotlib axes with the I-V curve plot.
    """
    bias_current = ic_params['bias']
    fb_mins_current = ic_params['fb_min']
    fb_maxes_current = ic_params['fb_max']
    max_sq1imod_idx = ic_params['bias_max_idx']
    start_sq1imod_idx = ic_params['bias_min_idx']

    ax.plot(bias_current, fb_mins_current, lw=2,
            label=min_label, color=min_color, alpha=alpha)
    ax.plot(bias_current, fb_maxes_current, lw=2,
            label=max_label, color=max_color, alpha=alpha)

    if(mod_label is not None):
        xs = [bias_current[max_sq1imod_idx], bias_current[max_sq1imod_idx]]
        ys = [fb_mins_current[max_sq1imod_idx],
              fb_maxes_current[max_sq1imod_idx]]
        ax.plot(xs, ys, lw=3,
                label=mod_label, color=mod_color, alpha=alpha)
    if(thresh_label is not None):
        xh = [0, bias_current[-1]]
        yh = [fb_mins_current[start_sq1imod_idx],
              fb_mins_current[start_sq1imod_idx]]
        xv = [bias_current[max_sq1imod_idx], bias_current[max_sq1imod_idx]]
        yv = [0, max(fb_maxes_current)]
        ax.plot(xh, yh, linestyle="dotted", lw=3,
                label=thresh_label, color=thresh_color, alpha=alpha)
        ax.plot(xv, yv, linestyle="dotted", lw=3,
                color=thresh_color, alpha=alpha)
    if(picked_bias_idx is not None):
        xh = [0, bias_current[picked_bias_idx]]
        yh = [fb_maxes_current[picked_bias_idx],
              fb_maxes_current[picked_bias_idx]]
        xv = [bias_current[picked_bias_idx], bias_current[picked_bias_idx]]
        yv = [0, fb_maxes_current[picked_bias_idx]]
        ax.plot(xh, yh, linestyle="dotted", lw=2,
                label=picked_label, color=picked_color, alpha=alpha)
        ax.plot(xv, yv, linestyle="dotted", lw=2,
                color=picked_color, alpha=alpha)

    return ax


def plot_icminmax(col, row, ic_params_rson, ic_params_rsoff=None,
                  ctime=None, convert_units=False,
                  savedir='output_data',
                  fig=None, ax=None,
                  show_plot=False):
    '''
    Plots the SQ1 Min and Max against bias current for a single row on a column
    '''
    alpha = 1
    figsize = (8, 6)
    sq1min_rson_label = 'SQ1 min, rs on'
    sq1min_rson_color = 'blue'
    sq1max_rson_label = 'SQ1 max, rs on'
    sq1max_rson_color = 'red'

    max_mod_color = 'purple'
    max_mod_label = ('$I^{SQ1}_{mod}$')

    if fig is None:
        fig, ax = plt.subplots(figsize=figsize)

    plot_iccurve(ic_params_rson, ax, alpha=alpha,
                 min_label=sq1min_rson_label, min_color=sq1min_rson_color,
                 max_label=sq1max_rson_label, max_color=sq1max_rson_color,
                 mod_label=max_mod_label, mod_color=max_mod_color)

    if ic_params_rsoff is not None:
        plot_iccurve(ic_params_rsoff, ax, alpha=alpha,
                     min_label='SQ1 min, rs off', min_color='aqua',
                     max_label='SQ1 max, rs off', max_color='lime',
                     thresh_label='Bias Limit', thresh_color='deeppink',
                     picked_bias_idx=ic_params_rsoff['bias_min_idx'],
                     picked_label=None, picked_color=None)

    if convert_units:
        ax.set_ylabel('SSA Input Current ($\mu$A)', fontsize=18)
        ax.set_xlabel('SQ1 Total Bias Current ($\mu$A)', fontsize=18)
        ax.set_ylim(0, 40)
    else:
        ax.set_ylabel('SSA FB (DAC)', fontsize=18)
        ax.set_xlabel('SQ1 Total Bias Current (DAC)', fontsize=18)
        ax.set_ylim(0, 10000)

    ax.set_title(str(ctime) + ' Ic Check Column ' +
                 str(col) + ' Row ' + str(row))
    fig.tight_layout()

    savename = str(ctime) + '_icminmax_units' + uname + \
        '_row' + str(row) + '_col' + str(col) + '.png'
    savename = os.path.join(savedir, savename)
    print('saving to: ' + savename)
    fig.savefig(savename)

    if show_plot:
        plt.show()

    return fig, ax


def plot_icminmax_col(last_fig, col, ic_params, ic_params2=None, ctime=None,
                      s1b_minmax_ax=None, s1b_minmax_fig=None, manual_bias_idx=None, convert_units=False,
                      show_plot=False, savedir='../output_data'):
    '''
    plots the ic col, ic min, and ic max given the ic_params. ic_params2 is assumed to be when row select is turned off
    Returns the axes so they can be continually plotted over for the entire column
    '''

    alpha = 0.1
    sq1_safb_servo_biases_uA = ic_params['bias']
    sq1_safb_servo_mins_sa_in_uA = ic_params['fb_min']
    sq1_safb_servo_maxs_sa_in_uA = ic_params['fb_max']
    max_sq1imod_idx = ic_params['bias_max_idx']
    start_sq1imod_idx = ic_params['bias_min_idx']
    start_sq1imod_uA = sq1_safb_servo_biases_uA[max_sq1imod_idx]
    if(s1b_minmax_ax is None):
        s1b_minmax_fig, s1b_minmax_ax = plt.subplots(figsize=(8, 6))
    if(last_fig):

        s1b_minmax_ax.plot(sq1_safb_servo_biases_uA, sq1_safb_servo_mins_sa_in_uA,
                           lw=2, label='SQ1 min, rs on', color='blue', alpha=alpha)
        s1b_minmax_ax.plot(sq1_safb_servo_biases_uA, sq1_safb_servo_maxs_sa_in_uA,
                           lw=2, label='SQ1 max, rs on', color='red', alpha=alpha)
        s1b_minmax_ax.plot([sq1_safb_servo_biases_uA[max_sq1imod_idx], sq1_safb_servo_biases_uA[max_sq1imod_idx]],
                           [sq1_safb_servo_mins_sa_in_uA[max_sq1imod_idx],
                            sq1_safb_servo_maxs_sa_in_uA[max_sq1imod_idx]], lw=3, color='purple', alpha=alpha)
        if(manual_bias_idx is not None):
            s1b_minmax_ax.plot([sq1_safb_servo_biases_uA[manual_bias_idx], sq1_safb_servo_biases_uA[manual_bias_idx]],
                               [0, sq1_safb_servo_maxs_sa_in_uA[max_sq1imod_idx]], lw=2, color='pink', alpha=1,
                               label='Manually Chosen Bias', linestyle='dotted')

        if(ic_params2 is not None):
            sq1_safb_servo_biases_uA = ic_params2['bias']
            sq1_safb_servo_mins_sa_in_uA = ic_params2['fb_min']
            sq1_safb_servo_maxs_sa_in_uA = ic_params2['fb_max']
            max_sq1imod_idx = ic_params2['bias_max_idx']
            start_sq1imod_idx = ic_params2['bias_min_idx']
            start_sq1imod_uA = sq1_safb_servo_mins_sa_in_uA[max_sq1imod_idx]
            s1b_minmax_ax.plot(sq1_safb_servo_biases_uA, sq1_safb_servo_mins_sa_in_uA,
                               lw=2, label='SQ1 min, rs off', color='aqua', alpha=alpha)
            s1b_minmax_ax.plot(sq1_safb_servo_biases_uA, sq1_safb_servo_maxs_sa_in_uA,
                               lw=2, label='SQ1 max, rs off', color='lime', alpha=alpha)

            bias_limit = sq1_safb_servo_biases_uA[start_sq1imod_idx]
            s1b_minmax_ax.plot([bias_limit, bias_limit],
                               [0, sq1_safb_servo_biases_uA[-1]], label='Bias Limit', color='deeppink', lw=3, linestyle="dotted")
            s1b_minmax_ax.plot([0, sq1_safb_servo_biases_uA[-1]],
                               [start_sq1imod_uA, start_sq1imod_uA],  color='deeppink', lw=3, linestyle="dotted")

        leg = s1b_minmax_ax.legend(loc='upper left', fontsize=8)
        for lh in leg.legendHandles:
            lh.set_alpha(1)
        s1b_minmax_ax.xaxis.set_minor_locator(AutoMinorLocator(5))
        if(convert_units):
            uname = 'ua'
            s1b_minmax_ax.set_ylabel('SSA Input Current ($\mu$A)', fontsize=18)
            s1b_minmax_ax.set_xlabel(
                'SQ1 Total Bias Current ($\mu$A)', fontsize=18)
            s1b_minmax_ax.set_ylim(0, 40)
        else:
            s1b_minmax_ax.set_ylabel('SSA FB (DAC)', fontsize=18)
            s1b_minmax_ax.set_xlabel(
                'SQ1 Total Bias Current (DAC)', fontsize=18)
            s1b_minmax_ax.set_ylim(0, 10000)
            uname = 'dac'
        s1b_minmax_fig.suptitle(
            str(ctime) + ' Ic Check Column ' + str(col) + ' All Rows')
        s1b_minmax_fig.tight_layout()
        savename = str(ctime) + '_icminmax_units'+uname + \
            '_summary_col' + str(col) + '.png'
        print('saving to: ' + os.path.join(savedir, savename))
        s1b_minmax_fig.savefig(os.path.join(savedir, savename))
        if(show_plot):
            s1b_minmax_fig.show()
        print("Figures open: " + str(plt.get_fignums()))
        plt.close('all')
        return None, None

    else:
        s1b_minmax_ax.plot(sq1_safb_servo_biases_uA, sq1_safb_servo_mins_sa_in_uA,
                           lw=2, color='blue', alpha=alpha)
        s1b_minmax_ax.plot(sq1_safb_servo_biases_uA, sq1_safb_servo_maxs_sa_in_uA,
                           lw=2,  color='red', alpha=alpha)

        s1b_minmax_ax.plot([sq1_safb_servo_biases_uA[max_sq1imod_idx],
                            sq1_safb_servo_biases_uA[max_sq1imod_idx]],
                           [sq1_safb_servo_mins_sa_in_uA[max_sq1imod_idx],
                            sq1_safb_servo_maxs_sa_in_uA[max_sq1imod_idx]],
                           lw=3, color='purple', alpha=alpha,
                           )
        bias_limit = sq1_safb_servo_biases_uA[start_sq1imod_idx]
        # s1b_minmax_ax.plot([bias_limit, bias_limit], [0,sq1_safb_servo_biases_uA[-1] ], alpha=alpha, color='orange')
        # s1b_minmax_ax.plot([0, sq1_safb_servo_biases_uA[-1]], [start_sq1imod_uA, start_sq1imod_uA ], alpha=alpha, color='orange')
        if(ic_params2 is not None):
            sq1_safb_servo_biases_uA = ic_params2['bias']
            sq1_safb_servo_mins_sa_in_uA = ic_params2['fb_min']
            sq1_safb_servo_maxs_sa_in_uA = ic_params2['fb_max']
            max_sq1imod_idx = ic_params2['bias_max_idx']
            start_sq1imod_idx = ic_params2['bias_min_idx']
            start_sq1imod_uA = sq1_safb_servo_biases_uA[max_sq1imod_idx]
            s1b_minmax_ax.plot(sq1_safb_servo_biases_uA, sq1_safb_servo_mins_sa_in_uA,
                               lw=2, color='aqua', alpha=alpha)
            s1b_minmax_ax.plot(sq1_safb_servo_biases_uA, sq1_safb_servo_maxs_sa_in_uA,
                               lw=2,  color='lime', alpha=alpha)

        return s1b_minmax_fig, s1b_minmax_ax


def plot_rsservo_col(last_fig, col, chip_num, sq1_params, sq1_params2=None, ctime=None,
                     s1b_minmax_ax=None, s1b_minmax_fig=None, show_plot=False):
    '''
    plots the rs servo at the given bias point. 
    '''
    colors = ['deepskyblue', 'red', 'green', 'purple']
    color = colors[chip_num]
    savedir = '../output_data/'
    alpha = 0.2
    (sq1_safb_servo_curves_dac, sq1_safb_servo_biases_dac,
     max_sq1_safb_servo_span_bias, max_sq1_safb_servo_span, sq1_safb_servo
     ) = sq1_params
    # print(sq1_safb_servo_curves_dac)
    # print(sq1_safb_servo)
    if(s1b_minmax_ax is None):
        s1b_minmax_fig, s1b_minmax_ax = plt.subplots(figsize=(8, 6))
    if(last_fig):

        # s1b_minmax_ax.plot(sq1_safb_servo_biases_uA, sq1_safb_servo_mins_sa_in_uA,
        #        lw=2, label='SQ1 min(Imod)', color='blue', alpha=alpha)
        s1b_minmax_ax.plot(
            sq1_safb_servo, sq1_safb_servo_curves_dac[0], alpha=alpha, color=color, label='Chip Number: ' + str(chip_num))
        if(sq1_params2 is not None):
            (sq1_safb_servo_curves_dac, sq1_safb_servo_biases_dac,
             max_sq1_safb_servo_span_bias, max_sq1_safb_servo_span, sq1_safb_servo
             ) = sq1_params2
            # s1b_minmax_ax.plot(sq1_safb_servo_biases_uA, sq1_safb_servo_mins_sa_in_uA,
            #    lw=2, label='Ic,col', color='aqua', alpha=alpha)
            s1b_minmax_ax.plot(
                sq1_safb_servo, sq1_safb_servo_curves_dac[0], alpha=alpha, color='aqua')
        handles, labels = plt.gca().get_legend_handles_labels()
        by_label = OrderedDict(zip(labels, handles))
        leg = s1b_minmax_ax.legend(
            by_label.values(), by_label.keys(), loc='upper left', fontsize=8)

        for lh in leg.legendHandles:
            lh.set_alpha(1)
        s1b_minmax_ax.set_ylabel('SSA Input Current (DAC units)', fontsize=18)
        s1b_minmax_ax.set_xlabel(
            'SQ1 Total Bias Current (DAC units)', fontsize=18)

        # s1b_minmax_ax.set_ylim(0, 40)
        s1b_minmax_fig.suptitle('RS Check Column ' + str(col))
        s1b_minmax_fig.tight_layout()
        savename = str(ctime) + '_rs_summary_col' + str(col) + '.png'
        print('saving to: ' + os.path.join(savedir, savename))
        s1b_minmax_fig.savefig(os.path.join(savedir, savename))
        if(show_plot):
            plt.show()
        plt.close()

    else:

        s1b_minmax_ax.plot(
            sq1_safb_servo, sq1_safb_servo_curves_dac[0],  alpha=alpha, color=color, label='Chip Number: ' + str(chip_num))

        # s1b_minmax_ax.plot([bias_limit, bias_limit], [0,sq1_safb_servo_biases_uA[-1] ], alpha=alpha, color='orange')
        # s1b_minmax_ax.plot([0, sq1_safb_servo_biases_uA[-1]], [start_sq1imod_uA, start_sq1imod_uA ], alpha=alpha, color='orange')
        if(sq1_params2 is not None):
            (sq1_safb_servo_curves_dac, sq1_safb_servo_biases_dac,
             max_sq1_safb_servo_span_bias, max_sq1_safb_servo_span, sq1_safb_servo
             ) = sq1_params2
            s1b_minmax_ax.plot(
                sq1_safb_servo, sq1_safb_servo_curves_dac[0], alpha=alpha, color='aqua')
        # plt.show()

    return s1b_minmax_fig, s1b_minmax_ax


def tile_plot(num_rows, num_columns, data, label, title, vmin=0, vmax=20,
              savedir='../output_data', show_plot=False):
    '''
    Assumes data to be plotted is accessed by data[row][col]
    '''
    fig, ax = plt.subplots()
    im = plt.imshow(data,
                    interpolation='none', aspect='equal',
                    vmin=vmin, vmax=vmax, cmap='bwr')

    ax = plt.gca()

    # Major ticks
    spacing = 5
    ax.set_xticks(np.arange(0, num_columns, spacing))
    ax.set_yticks(np.arange(0, num_rows, spacing))

    # Labels for major ticks
    ax.set_xticklabels(np.arange(0, num_columns, spacing))
    ax.set_yticklabels(np.arange(0, num_rows, spacing))

    # Minor ticks
    ax.set_xticks(np.arange(-.5, num_columns, 1), minor=True)
    ax.set_yticks(np.arange(-.5, num_rows, 1), minor=True)

    ax.set_xlabel('Column')
    ax.set_ylabel('Row')
    # Gridlines based on minor ticks
    ax.grid(which='minor', color='w', linestyle='-', linewidth=2)

    # Remove minor ticks
    ax.tick_params(which='minor', bottom=False, left=False)
    fig.suptitle(title)
    cbar = fig.colorbar(im)
    cbar.set_label(label)
    savename = os.path.join(savedir, title + '.png')
    print('saving: ' + savename)
    plt.savefig(savename)
    if(show_plot):
        plt.show()
    plt.close('all')
