import numpy as np  

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