#*************ESTIMATES PHOTON-LIMITED VELOCITY UNCETAINTIES FOR KPF*************
#
#	This function is a simple interpolator for a pre-calculated grid of photon-limited RV uncertainties.
#	The four dimensions to the 'master' grid are: 
#
#	- Stellar effective temperature [K]
#	- V-band magnitude
#	- Exposure time [s]
#	- KPF echelle order
#
#	The precomputed sigma_rv grid was calculated assuming parameter bounds of:
#	- Teff: 2700 - 6600 Kelvin
#	- Vmag: 5 - 19
#	- Texp: 10 - 3600 seconds
#	Interpolating beyond these bounds will result in errors in the returned RV values
#
#	The interpolator loads a set of pre-defined .fits files with the relevant grid, basis vectors (Teff, Vmag, Texp)
#
#*************DETAILS ON FITS FILES*************
#
#	Grid files:  
#	'dv_uncertainty_master_order.fits':  order-by-order velocity uncertainites [m/s] (4D grid)
#   'dv_uncertainty_master.fits': 		 integrated KPF velocity uncertainty values [m/s] (3D grid)
#	'snr_master_order.fits':		  	 mean order-by-order SNR values (4d grid)
#
#	Basis vectors:
#	'order_wvl_centers.fits': 			 mean wavelengths of echelle orders [nm]
#	'photon_grid_teff.fits':			 stellar effective temperature array [K]
#	'photon_grid_vmag.fits':			 stellar V magnitude
#	'photon_grid_exptime.fits':			 exposure time array [s]
#
#**************************
#
#	Function inputs:
#	- Teff:	Stellar effective temperature [K]
#	- Vmag:	V magnitude of target
#	- exp_time:	Exposure time [seconds]
#	
#	Keywords:
#	- fits_dir: Directory location of fits files, including master grid + basis arrays
#	- save_file: Flag keyword for saving order-by-order sigma_rv values to a text file
#	- quiet: Suppresses text output in terminal
#
#	Outputs:
#	- dv: Total sigma_rv over the KPF bandpass
#       - orderinfo: tuple containing 
#	        Mean wavelength for each echelle order	[nm] (array)
#               Mean SNR for each order in the KPF bandpass (array)
#	        Estimated velocity uncertainty for each echelle order [m/s] (array)
#
#Example call: dv, orderinfo = KPF_photon_noise_estimate(5400.,8.,600.,fits_dir='grids/')
#

# We use astropy for reading fits files, since PyFITS has a planned deprecation
#       see http://docs.astropy.org/en/stable/io/fits/appendix/faq.html
from astropy.io import fits
from scipy import interpolate
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.interpolate import RegularGridInterpolator
import numpy as np


def KPF_photon_noise_estimate(Teff, Vmag, exp_time,
                              fits_dir='./', quiet=False, save_file=True):

    # Grid files for teff, vmag, exp_time
    teff_grid_file = fits_dir + 'photon_grid_teff.fits'
    vmag_grid_file = fits_dir + 'photon_grid_vmag.fits'
    exp_time_grid_file = fits_dir + 'photon_grid_exptime.fits'

    teff_grid = fits.open(teff_grid_file)[0].data
    vmag_grid = fits.open(vmag_grid_file)[0].data
    exp_time_grid = fits.open(exp_time_grid_file)[0].data

    # Master grid files for interpolation
    SNR_grid_file = fits_dir + 'snr_master_order.fits'
    sigma_rv_grid_file = fits_dir + 'dv_uncertainty_master_order.fits'
    sigma_rv_total_file = fits_dir + 'dv_uncertainty_master.fits'
    wvl_ord_file = fits_dir + 'order_wvl_centers.fits'
    
    #read in data cubes and arrays
    #    sigma_rv_grid_ord ( Teff, Vmag, Exposure time, Echelle order)	- order-binned velocity uncertainty
    #    SNR_grid_ord ( Teff, Vmag, Exposure time, Echelle order)	- mean spectral SNR within order
    #	 sigma_rv_grid ( Teff, Vmag, Exposure time )	- total
    snr_grid_ord = fits.open(SNR_grid_file)[0].data
    sigma_rv_grid_ord = fits.open(sigma_rv_grid_file)[0].data
    sigma_rv_grid = fits.open(sigma_rv_total_file)[0].data
    wvl_ords = fits.open(wvl_ord_file)[0].data
    
    # Separate order and wavelength vectors
    order_arr = wvl_ords[0]
    wvl_arr = wvl_ords[1]

## This is a direct translation of th IDL routine
    # Get fractional indices for input parameters 
    teff_index_spline = InterpolatedUnivariateSpline(teff_grid, 
                            np.arange(len(teff_grid), dtype=np.double))
    teff_location = teff_index_spline(Teff)
    vmag_index_spline = InterpolatedUnivariateSpline(vmag_grid, 
                            np.arange(len(vmag_grid), dtype=np.double))
    vmag_location = vmag_index_spline(Vmag)
    # Exposure time in log space
    exp_time_index_spline = InterpolatedUnivariateSpline(np.log10(exp_time_grid), 
                                np.arange(len(exp_time_grid), dtype=np.double))
    exp_time_location = exp_time_index_spline(np.log10(exp_time))

    # Interpolate sigma_rv
    interpolation_grid = (np.arange(len(exp_time_grid)),
                          np.arange(len(vmag_grid)),
                          np.arange(len(teff_grid)))
    sigma_rv_interpolator = RegularGridInterpolator( interpolation_grid,
                                sigma_rv_grid)
    fit_point = [exp_time_location, vmag_location, teff_location]
    sigma_rv_val = sigma_rv_interpolator(fit_point)[0]
### This is more straightforward, but doesn't account for the non-linear spacing of 
###      log10(exp_time_grid) 
#    sigma_rv_interpolator = RegularGridInterpolator(
#                                (np.log10(exp_time_grid),vmag_grid,teff_grid),
#                                sigma_rv_grid)
#    # Total sigma_RV value
#    sigma_rv_val = sigma_rv_interpolator([np.log10(exp_time), Vmag, Teff])[0]

    # Compute uncertainties per order
    Norders = len(order_arr)
    sigma_rv_ord = np.empty(Norders)
    snr_rv_ord = np.empty(Norders)

    if not quiet:
        print('')
        print('--------------------------')
        print('Order #   Wavelength [nm]   Mean SNR     sigma_rv [m/s]')
    # Interpolate grids for each order
    for l in range(Norders):
        # Sigma RV interpolation
        sigma_rv_grid_ord_l = sigma_rv_grid_ord[l] 
        sigma_rv_interpolator_ord_l = RegularGridInterpolator( interpolation_grid,
                                          sigma_rv_grid_ord_l)
        sigma_rv_val_ord_l = sigma_rv_interpolator_ord_l(fit_point)[0]
        sigma_rv_ord[l] = sigma_rv_val_ord_l
        
        # SNR grid interpolation
        snr_grid_ord_l = snr_grid_ord[l] 
        snr_interpolator_ord_l = RegularGridInterpolator( interpolation_grid,
                                          snr_grid_ord_l)
        snr_val_ord_l = snr_interpolator_ord_l(fit_point)[0]
        snr_rv_ord[l] = snr_val_ord_l
        if not quiet:
            print('%6.1lf   %8.1lf    %12.2lf %15.3lf' % 
                     (order_arr[l], wvl_arr[l], snr_rv_ord[l], sigma_rv_ord[l]))
    if not quiet:
        print('--------------------------')
        print('')
        print('Total velocity uncertainty: %lf m/s' % sigma_rv_val)
        print('')

    snr_ord = snr_rv_ord
    dv_ord = sigma_rv_ord


    if save_file:

        # Parameter strings
        exp_time_str = str(np.int(np.round(exp_time))) + 's'
        teff_str = str(np.int(np.round(Teff))) + 'K'
        vmag_str = str(np.round(Vmag,2))
        sigma_rv_str = str(np.round(sigma_rv_val, 3))
        vsini_str = '2kms'

        save_file_name = 'dv_photon_' + teff_str + '_' + vmag_str + '_' + exp_time_str + '_' + vsini_str + '.txt'

        header = 'Teff: ' + teff_str + '\n'
        header += 'Vmag: ' + vmag_str + '\n'
        header += 't_exposure: ' + exp_time_str + '\n'
        header += 'vsini: ' + vsini_str + '\n'
        header += ' ' + '\n'
        header += '--------------------------' + '\n'
        header += '' + '\n'
        header += 'Total velocity uncertainty: ' + str(np.round(sigma_rv_val, 2)) + '\n'
        header += '' + '\n'
        header += 'Wavelength [nm]        Mean SNR per pixel   Sigma_rv [m/s]' + '\n'
        # Save an ASCII file with wavelength, SNR, dv estimate
        np.savetxt(save_file_name, 
                      np.transpose([wvl_arr, snr_rv_ord, sigma_rv_ord]),
                      fmt=('%20.3lf %20.3lf %20.3lf'), header=header)

    return (sigma_rv_val, (wvl_arr, snr_ord, dv_ord))


#dv, orderinfo = KPF_photon_noise_estimate(5400.,8.,600.,fits_dir='grids/')

















