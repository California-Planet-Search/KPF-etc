# KPF-etc
README.md

This repository contains a set of pre-calculated photon-limited velocity uncertainties for the Keck Planet Finder project, as well as a simple Python lookup function for quickly estimating photon-limited RV precision for a user-specified target observation. 

In addition to a photon noise estimator (etc.kpf_photon_noise_estimate), there are also functions that calculate the required exposure time to reach a desired photon-limited RV uncertainty (etc.kpf_etc_rv), and to reach a desired spectral signal-to-noise value at a specific wavelength (etc.kpf_etc_rv).

To install, simply download the package, change directories into the downloaded folder, and run:

	pip install .

Example uses below:

	from kpf_etc.etc import kpf_photon_noise_estimate, kpf_etc_rv, kpf_etc_snr

	# calculate photon noise estimate for a known exposure length
	exp_time = 46. # s
	vmag = 8.
	teff = 5000. # K
	sigma_rv_val, wvl_arr, snr_ord, dv_ord = kpf_photon_noise_estimate(teff,vmag,exp_time)  

	# exposure time estimate to reach a desired RV precision
	sigma_rv_desired = 0.5 # m/s
	exposure_time_sigma_rv = kpf_etc_rv(teff, vmag, sigma_rv_desired)

	# exposure time estimate to reach a desired spectral SNR
	snr_desired = 500.
	wavelength_desired = 550. # nm
	exposure_time_snr = kpf_etc_snr(teff, vmag, snr_desired, wavelength_desired)


The model grid descriptions (located in kpf_etc/grids) are as follows:

	'dv_uncertainty_master_order.fits':  order-by-order velocity uncertainites [m/s] (4D grid, Teff x Vmag x Exposure time x Echelle order)
    'dv_uncertainty_master.fits': 		 integrated KPF velocity uncertainty values [m/s] (3D grid, Teff x Vmag x Exposure time)
	'snr_master_order.fits':		  	 mean order-by-order SNR values (4d grid, Teff x Vmag x Exposure time x Echelle order)

The basis arrays used to generate the grids are contained in these files:

	'order_wvl_centers.fits': 			 mean wavelengths of echelle orders [nm]
	'photon_grid_teff.fits':			 stellar effective temperature array [K]
	'photon_grid_vmag.fits':			 stellar V magnitude
	'photon_grid_exptime.fits':			 exposure time array [s]

The precomputed grids were calculated assuming parameter bounds of:

	- Teff: 2700 - 6600 Kelvin
	- Vmag: 2 - 19
	- Texp: 10 - 3600 seconds
	
Interpolating beyond these bounds will result in errors in the returned RV values.
	
All calculations are based on BT-Settl synthetic stellar spectra, available at: http://phoenix.astro.physik.uni-goettingen.de/?page_id=15

The current uncertainty grid uses the latest KPF PDR throughput model (30 micron STA detector, green & red camera), a median WMKO atmospheric transmission curve (Buton et al 2012), and excludes regions of the spectrum with 1% or deeper telluric features (based on telfit model). A vsini value of 2 km/s is assumed for all stellar targets, and the observation airmass is fixed at 1.2.

Written by Sam Halverson (sphalverson@gmail.com) and Sean Mills