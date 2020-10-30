import os
from astropy.io import fits
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.interpolate import RegularGridInterpolator
import numpy as np

# sort out paths
LOCALPATH = os.path.dirname(os.path.realpath(__file__))

# directory for fits files
FITS_DIR = os.path.join(LOCALPATH, 'grids')

#----------------------------------------
def _findel(num, arr):

    '''
    Finds index of nearest element in array to a given number

    Parameters
    ----------
    num : int/float
        number to search for nearest element in array

    arr : :obj:`ndarray` of :obj:`float`
        array to search for 'num'

    Returns
    -------
    idx : int
        index of array element with value closes to 'num'

    S Halverson - JPL - 29-Sep-2019
    '''
    arr = np.array(arr)
    idx = (np.abs(arr - num)).argmin()
    return idx
#----------------------------------------

def kpf_photon_noise_estimate(teff, vmag, exp_time,
                              quiet=False):
    '''
    Estimates the approximate photon-limited radial velocity uncertainty
    for a given KPF target and exposure.

	This function is a simple interpolator for a pre-calculated grid of
    photon-limited RV uncertainties. The four dimensions to the 'master'
    grid are:

	- Stellar effective temperature [K]
	- V-band magnitude
	- Exposure time [s]
	- KPF echelle order

	The precomputed sigma_rv grid was calculated assuming parameter bounds of:
	- Teff: 2700 - 6600 Kelvin
	- Vmag: 5 - 19
	- Texp: 10 - 3600 seconds
	Interpolating beyond these bounds will result in errors in the
    returned RV values

	The interpolator loads a set of pre-defined .fits files with the
    relevant grid, basis vectors (Teff, Vmag, Texp)

    Parameters
    ----------
    teff : :obj:`float`
        Target stellar effective temperature [K]

    vmag : :obj:`float`
        Target V magnitude

	exp_time : :obj:`float`
    	Target exposure time [seconds]

	Keywords:
    ----------
    quiet : boolean
        Suppresses text output in terminal

	Outputs
    ----------
	sigma_rv : :obj:`float`
        Total photon-limited RV uncertainty over KPF bandpass

    wvl_arr : :obj:`arr` of :obj:`float`
        Echelle order central wavelengths

    snr_rv_ord : :obj:`arr` of :obj:`float`
        Median SNR values for each echelle order

    sigma_rv_ord : :obj:`arr` of :obj:`float`
        Photon-limited uncertainties for each order

    '''

    # Grid files for teff, vmag, exp_time
    teff_grid_file = os.path.join(FITS_DIR, 'photon_grid_teff.fits')
    vmag_grid_file = os.path.join(FITS_DIR, 'photon_grid_vmag.fits')
    exp_time_grid_file = os.path.join(FITS_DIR, 'photon_grid_exptime.fits')

    teff_grid = fits.getdata(teff_grid_file)
    vmag_grid = fits.getdata(vmag_grid_file)
    exp_time_grid = fits.getdata(exp_time_grid_file)

    # Master grid files for interpolation
    snr_grid_file = os.path.join(FITS_DIR, 'snr_master_order.fits')
    sigma_rv_grid_file = os.path.join(FITS_DIR, 'dv_uncertainty_master_order.fits')
    sigma_rv_total_file = os.path.join(FITS_DIR, 'dv_uncertainty_master.fits')
    wvl_ord_file = os.path.join(FITS_DIR, 'order_wvl_centers.fits')

    # read in data cubes and arrays
    snr_grid_ord = fits.getdata(snr_grid_file)
    sigma_rv_grid_ord = fits.getdata(sigma_rv_grid_file)
    sigma_rv_grid = fits.getdata(sigma_rv_total_file)
    wvl_ords = fits.getdata(wvl_ord_file)

    # check that inputs are within the grid boundaries
    flag_bound = True
    if teff < np.min(teff_grid) or teff > np.max(teff_grid):
        print("Temperature out of bounds (%d K to %d K)" % (np.amin(teff_grid), np.amax(teff_grid)))
        flag_bound = False
    if vmag < np.min(vmag_grid) or vmag > np.max(vmag_grid):
        print("Magnitude out of bounds (V = %d to V = %d)" % (np.amin(vmag_grid), np.amax(vmag_grid)))
        flag_bound = False
    if not flag_bound:
        return np.nan

    # Separate order and wavelength vectors
    order_arr = wvl_ords[0]
    wvl_arr = wvl_ords[1]

    # Get fractional indices for input parameters
    teff_index_spline = InterpolatedUnivariateSpline(teff_grid,
                                                     np.arange(len(teff_grid),
                                                               dtype=np.double))
    teff_location = teff_index_spline(teff)

    vmag_index_spline = InterpolatedUnivariateSpline(vmag_grid,
                                                     np.arange(len(vmag_grid),
                                                               dtype=np.double))

    vmag_location = vmag_index_spline(vmag)

    # Exposure time in log space
    exp_time_index_spline = InterpolatedUnivariateSpline(np.log10(exp_time_grid),
                                                         np.arange(len(exp_time_grid),
                                                                   dtype=np.double))
    exp_time_location = exp_time_index_spline(np.log10(exp_time))

    # Interpolate sigma_rv
    interpolation_grid = (np.arange(len(exp_time_grid)),
                          np.arange(len(vmag_grid)),
                          np.arange(len(teff_grid)))

    sigma_rv_interpolator = RegularGridInterpolator(interpolation_grid,
                                                    sigma_rv_grid)

    fit_point = [exp_time_location, vmag_location, teff_location]
    sigma_rv = sigma_rv_interpolator(fit_point)[0]

    # Compute uncertainties per order
    sigma_rv_ord = np.empty(len(order_arr))
    snr_rv_ord = np.empty(len(order_arr))

    if not quiet:
        print('')
        print('--------------------------')
        print('Order #   Wavelength [nm]   Mean SNR     sigma_rv [m/s]')

    # Interpolate grids for each order
    for ind, ord_val in enumerate(order_arr):
        # Sigma RV interpolation
        sigma_rv_grid_ord_l = sigma_rv_grid_ord[ind]
        sigma_rv_interpolator_ord_l = RegularGridInterpolator(interpolation_grid,
                                                              sigma_rv_grid_ord_l)
        sigma_rv_val_ord_l = sigma_rv_interpolator_ord_l(fit_point)[0]
        sigma_rv_ord[ind] = sigma_rv_val_ord_l

        # SNR grid interpolation
        snr_grid_ord_l = snr_grid_ord[ind]
        snr_interpolator_ord_l = RegularGridInterpolator(interpolation_grid,
                                                         snr_grid_ord_l)
        snr_val_ord_l = snr_interpolator_ord_l(fit_point)[0]
        snr_rv_ord[ind] = snr_val_ord_l
        if not quiet:
            print('%6.1lf   %8.1lf    %12.2lf %15.3lf' %
                  (ord_val, wvl_arr[ind], snr_rv_ord[ind], sigma_rv_ord[ind]))
    if not quiet:
        print('--------------------------')
        print('')
        print('Total velocity uncertainty: %lf m/s' % sigma_rv)
        print('')

    return sigma_rv, wvl_arr, snr_rv_ord, sigma_rv_ord

def kpf_etc_rv(teff, vmag, sigma_rv):
    '''
    Estimates the exposure time required to reach a specified RV uncertainty
    value (sigma_rv) for a given stellar target. The target is defined by the
    stellar effective temperature (teff) and V magnitude (vmag)

    This function interpolates over a pre-computed grid of RV uncertainty values
    to estimate the optimum exposure length. This is done by looping through
    'trial' exposure time guesses to see which yields sigma_rv closest to the
    desired value.

    Parameters
    ----------
    teff : :obj:`float`
        Target effective temperature

    vmag : :obj:`float`
        Target V magnitude

    sigma_rv : :obj:`float`
        Desired integrated RV uncertainty (photons only)

    Returns
    -------
    exptime : :obj:`float`
        Estimated exposure time for reaching specified sigma_rv

    S Halverson - JPL - 29-Oct-2020
    '''

    # Grid files for teff, vmag, exp_time
    teff_grid_file = os.path.join(FITS_DIR, 'photon_grid_teff.fits')
    vmag_grid_file = os.path.join(FITS_DIR, 'photon_grid_vmag.fits')
    exp_time_grid_file = os.path.join(FITS_DIR, 'photon_grid_exptime.fits')

    # Master grid files for interpolation
    sigma_rv_total_file = os.path.join(FITS_DIR, 'dv_uncertainty_master.fits')

    # read in data cubes and arrays
    sigma_rv_grid = fits.getdata(sigma_rv_total_file)
    teff_grid = fits.getdata(teff_grid_file)
    vmag_grid = fits.getdata(vmag_grid_file)
    exptime_grid = fits.getdata(exp_time_grid_file)
    logexp = np.log10(exptime_grid)

    # check that inputs are within the grid boundaries
    flag_bound = True
    if teff < np.min(teff_grid) or teff > np.max(teff_grid):
        print("Temperature out of bounds (%d K to %d K)" % (np.amin(teff_grid), np.amax(teff_grid)))
        flag_bound = False
    if vmag < np.min(vmag_grid) or vmag > np.max(vmag_grid):
        print("Magnitude out of bounds (V = %d to V = %d)" % (np.amin(vmag_grid), np.amax(vmag_grid)))
        flag_bound = False
    if not flag_bound:
        return np.nan

    # Get fractional indices for relevant input parameters
    teff_index_spline = InterpolatedUnivariateSpline(teff_grid,
                                                     np.arange(len(teff_grid),
                                                               dtype=np.double))
    teff_location = teff_index_spline(teff)

    vmag_index_spline = InterpolatedUnivariateSpline(vmag_grid,
                                                     np.arange(len(vmag_grid),
                                                               dtype=np.double))
    vmag_location = vmag_index_spline(vmag)

    # dummy criterea variables for loop
    ind = 2
    maxout = 1e10

    # while trial exposure time yields worse precision, keep increasing until
    # you reach specified sigma_rv
    while maxout > sigma_rv:

        # dummy guess trial exposure
        trial_exp = min(exptime_grid) + ind

        # fractional index for trial_exp in exptime_grid
        exptime_index = InterpolatedUnivariateSpline(logexp,
                                                     np.arange(len(exptime_grid),
                                                               dtype=np.double))(np.log10(trial_exp))

        # recompute expected sigma_rv based on trial exposure time
        sigma_rv_interpolator = RegularGridInterpolator((np.arange(len(exptime_grid)),
                                                         np.arange(len(vmag_grid)),
                                                         np.arange(len(teff_grid))),
                                                        sigma_rv_grid)

        inputs = [exptime_index, vmag_location, teff_location]

        # store as new maximum
        maxout = sigma_rv_interpolator(inputs)[0]

        # increase exposure time by 1 second
        ind += 1

    # last 'trial' exposure time is correct answer
    exptime = trial_exp

    return exptime


def kpf_etc_snr(teff, vmag, snr, wavelength):
    '''
    Estimates the exposure time required to reach a specified signal-to-noise
    value (snr) at a specified wavelength for a given stellar target. The
    target is defined by the stellar effective temperature (teff) and V
    magnitude (vmag)

    This function interpolates over a pre-computed grid of RV uncertainty values
    to estimate the optimum exposure length. This is done by looping through
    'trial' exposure time guesses to see which yields sigma_rv closest to the
    desired value.

    Parameters
    ----------
    teff : :obj:`float`
        Target effective temperature

    vmag : :obj:`float`
        Target V magnitude

    snr : :obj:`float`
        Desired spectral SNR

    wavelength : :obj:'float'
        Wavelength of desired SNR

    Returns
    -------
    exptime : :obj:`float`
        Estimated exposure time for reaching specified snr

    S Halverson - JPL - 29-Oct-2020
    '''

    # Grid files for teff, vmag, exp_time
    teff_grid_file = os.path.join(FITS_DIR, 'photon_grid_teff.fits')
    vmag_grid_file = os.path.join(FITS_DIR, 'photon_grid_vmag.fits')
    exp_time_grid_file = os.path.join(FITS_DIR, 'photon_grid_exptime.fits')
    wvl_ord_file = os.path.join(FITS_DIR, 'order_wvl_centers.fits')

    # Master grid files for interpolation
    snr_grid_file = os.path.join(FITS_DIR, 'snr_master_order.fits')
    snr_grid_all = fits.getdata(snr_grid_file)

    # find closest order to specified wavelength
    wvl_ords = fits.getdata(wvl_ord_file)[1]
    ord_ind = _findel(wavelength, wvl_ords)
    snr_grid = snr_grid_all[ord_ind]

    # read in data cubes and arrays
    teff_grid = fits.getdata(teff_grid_file)
    vmag_grid = fits.getdata(vmag_grid_file)
    exptime_grid = fits.getdata(exp_time_grid_file)
    logexp = np.log10(exptime_grid)

    # check that inputs are within the grid boundaries
    flag_bound = True
    if teff < np.min(teff_grid) or teff > np.max(teff_grid):
        print("Temperature out of bounds (%d K to %d K)" % (np.amin(teff_grid), np.amax(teff_grid)))
        flag_bound = False
    if vmag < np.min(vmag_grid) or vmag > np.max(vmag_grid):
        print("Magnitude out of bounds (V = %d to V = %d)" % (np.amin(vmag_grid), np.amax(vmag_grid)))
        flag_bound = False
    if not flag_bound:
        return np.nan

    # Get fractional indices for relevant input parameters
    teff_index_spline = InterpolatedUnivariateSpline(teff_grid,
                                                     np.arange(len(teff_grid),
                                                               dtype=np.double))
    teff_location = teff_index_spline(teff)

    vmag_index_spline = InterpolatedUnivariateSpline(vmag_grid,
                                                     np.arange(len(vmag_grid),
                                                               dtype=np.double))
    vmag_location = vmag_index_spline(vmag)

    # dummy criterea variables for loop
    ind = 2
    minout = 0.

    # while trial exposure time yields worse precision, keep increasing until
    # you reach specified sigma_rv
    while minout < snr:

        # dummy guess trial exposure
        trial_exp = min(exptime_grid) + ind

        # fractional index for trial exposure time in exptime_grid
        exptime_index = InterpolatedUnivariateSpline(logexp, np.arange(len(exptime_grid),
                                                                       dtype=np.double))(np.log10(trial_exp))

        # recompute expected SNR based on trial exposure time
        snr_interpolator = RegularGridInterpolator((np.arange(len(exptime_grid)),
                                                    np.arange(len(vmag_grid)),
                                                    np.arange(len(teff_grid))),
                                                   snr_grid)
        inputs = [exptime_index, vmag_location, teff_location]

        # store as new maximum
        minout = snr_interpolator(inputs)[0]

        # increase exposure time by 1 second
        ind += 1

    # last 'trial' exposure time is correct answer
    exptime = trial_exp
    return exptime
