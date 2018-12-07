;*************ESTIMATES PHOTON-LIMITED VELOCITY UNCETAINTIES FOR KPF*************
;
;	This function is a simple interpolator for a pre-calculated grid of photon-limited RV uncertainties.
;	The four dimensions to the 'master' grid are: 
;
;	- Stellar effective temperature [K]
;	- V-band magnitude
;	- Exposure time [s]
;	- KPF echelle order
;
	;	The precomputed sigma_rv grid was calculated assuming parameter bounds of:
	;	- Teff: 2700 - 6600 Kelvin
	;	- Vmag: 5 - 19
	;	- Texp: 10 - 3600 seconds
	;	Interpolating beyond these bounds will result in errors in the returned RV values
;
;	The interpolator loads a set of pre-defined .fits files with the relevant grid, basis vectors (Teff, Vmag, Texp)
;
;*************DETAILS ON FITS FILES*************
;
;	Grid files:  
;	'dv_uncertainty_master_order.fits':  order-by-order velocity uncertainites [m/s] (4D grid)
;   'dv_uncertainty_master.fits': 		 integrated KPF velocity uncertainty values [m/s] (3D grid)
;	'snr_master_order.fits':		  	 mean order-by-order SNR values (4d grid)
;
;	Basis vectors:
;	'order_wvl_centers.fits': 			 mean wavelengths of echelle orders [nm]
;	'photon_grid_teff.fits':			 stellar effective temperature array [K]
;	'photon_grid_vmag.fits':			 stellar V magnitude
;	'photon_grid_exptime.fits':			 exposure time array [s]
;
;**************************
;
;	Function inputs:
;	- Teff:	Stellar effective temperature [K]
;	- Vmag:	V magnitude of target
;	- Exp_time:	Exposure time [seconds]
;	
;	Keywords:
;	- FITS_DIR: Directory location of fits files, including master grid + basis arrays
;	- SAVE_FILE: Flag keyword for saving order-by-order sigma_rv values to a text file
;
;	Outputs:
;	- SNR_ORD:	Mean SNR for each order in the KPF bandpass (array)
;	- DV_ORD:	Estimated velocity uncertainty for each echelle order [m/s] (array)
;	- WVL_ARR:	Mean wavelength for each echelle order	[nm] (array)
;
;Example call: dv = KPF_PHOTON_NOISE_ESTIMATE(5400d,8d,600d,SNR_ORD,DV_ORD,FITS_DIR='grids/')
;
FUNCTION KPF_PHOTON_NOISE_ESTIMATE,TEFF,VMAG,EXP_TIME, $
	SNR_ORD,DV_ORD,WVL_ARR, $
	FITS_DIR=FITS_DIR,SAVE_FILE=SAVE_FILE

;fits files directory containing all master data cubes, basis vectors
IF NOT KEYWORD_SET(FITS_DIR) THEN FITS_DIR = './'

;grid files for teff, vmag, exp_time
teff_grid_file = FITS_DIR + 'photon_grid_teff.fits'
vmag_grid_file = FITS_DIR + 'photon_grid_vmag.fits'
exptime_grid_file = FITS_DIR + 'photon_grid_exptime.fits'

;read in the various data cube axis
teff_grid = READFITS(teff_grid_file,/SILENT)
vmag_grid = READFITS(vmag_grid_file,/SILENT)
exp_time_grid = READFITS(exptime_grid_file,/SILENT)

;master grid files (fits) -- contains master data cubes to interpolate over
SNR_grid_file = FITS_DIR + 'snr_master_order.fits'
sigma_rv_grid_file = FITS_DIR + 'dv_uncertainty_master_order.fits'
sigma_rv_total_file = FITS_DIR + 'dv_uncertainty_master.fits'
wvl_ord_file = FITS_DIR + 'order_wvl_centers.fits'

;read in data cubes and arrays
;    sigma_rv_grid_ord ( Teff, Vmag, Exposure time, Echelle order)	- order-binned velocity uncertainty
;    SNR_grid_ord ( Teff, Vmag, Exposure time, Echelle order)	- mean spectral SNR within order
;	 sigma_rv_grid ( Teff, Vmag, Exposure time )	- total
snr_grid_ord = READFITS(SNR_grid_file,/SILENT)
sigma_rv_grid_ord = READFITS(sigma_rv_grid_file,/SILENT)
sigma_rv_grid = READFITS(sigma_rv_total_file,/SILENT)
wvl_ords = READFITS(wvl_ord_file,/SILENT)

;separate order and wavelength vectors
order_arr = wvl_ords[*,0]
wvl_arr = wvl_ords[*,1]

;get fractional indecies for input parameters (Teff, vmag, exp_time)
;teff_location = INTERPOL(FINDGEN(N_ELEMENTS(teff_arr_master)), ALOG10(teff_arr_master), ALOG10(teff),/SPLINE)
teff_location = INTERPOL(FINDGEN(N_ELEMENTS(teff_grid)), teff_grid, teff, /SPLINE)
vmag_location = INTERPOL(FINDGEN(N_ELEMENTS(vmag_grid)), vmag_grid, VMAG, /SPLINE)

;get fractional index of exposure time in log space, since the sampling is logarithmic.
exptime_location = INTERPOL(FINDGEN(N_ELEMENTS(exp_time_grid)), ALOG10(exp_time_grid), ALOG10(exp_time), /SPLINE)

;total sigma_rv value
sigma_rv_val = INTERPOLATE(sigma_rv_grid, teff_location, vmag_location, exptime_location, /GRID, /DOUBLE)

;make vectors for specific input parameters
sigma_rv_ord = DBLARR(N_ELEMENTS(order_arr))
snr_rv_ord = DBLARR(N_ELEMENTS(order_arr))

;for each order, interpolate grids
FOR l = 0, N_ELEMENTS(order_arr) - 1 DO BEGIN
	;sigma_rv interpolation
	sigma_rv_grid_ord_l = sigma_rv_grid_ord[*,*,*,l]
	sigma_rv_val_l = INTERPOLATE(sigma_rv_grid_ord_l, teff_location, vmag_location, exptime_location, /GRID, /DOUBLE)
	sigma_rv_ord[l] = sigma_rv_val_l
	
	;snr grid interpolation
	snr_grid_ord_l = snr_grid_ord[*,*,*,l]
	snr_rv_val_l = INTERPOLATE(snr_grid_ord_l, teff_location, vmag_location, exptime_location, /GRID, /DOUBLE)
	snr_rv_ord[l] = snr_rv_val_l
;	PRINT,STRCOMPRESS(STRING(order_arr[l],format='(f20.1)'),/REMOVE),'		' $
;		,STRCOMPRESS(STRING(wvl_arr[l],format='(f20.1)'),/REMOVE),'		' $
;		,STRCOMPRESS(STRING(snr_rv_ord[l],format='(f20.2)'),/REMOVE),'		' $
;		,STRCOMPRESS(STRING(sigma_rv_ord[l],format='(f20.3)'),/REMOVE)

ENDFOR
	
;output arrays for mean SNR and velocity uncertainty across orders
SNR_ORD = snr_rv_ord
DV_ORD = sigma_rv_ord

;-------------PARAMETER STRINGS-----------------
;setup strings for exposure time, stellar parameters
exp_str_for = '(i0' + STRCOMPRESS(STRING(FLOOR(ALOG10(EXP_TIME)+1d),format='(i01)'),/remove) + ')'
exp_time_str = STRCOMPRESS(STRING(EXP_TIME,format=exp_str_for),/remove) + 's'

;teff string
teff_str_for = '(i0' + STRCOMPRESS(STRING(FLOOR(ALOG10(teff)+1d),format='(i01)'),/remove) + ')'
teff_str = STRCOMPRESS(STRING(teff,format=teff_str_for),/remove) + 'K'

;v_mag string
v_mag_str_for = '(f20.1)'
v_mag_str = STRCOMPRESS(STRING(VMAG,format=v_mag_str_for),/remove)

;vsini string
vsini_str = '2 kms'

;sigma_rv total string
sigma_rv_str_fot = '(f10.3)'
sigma_rv_str = STRCOMPRESS(STRING(sigma_rv_val,format=sigma_rv_str_fot),/remove) + ' m/s'
;----------------------------------

;save ascii file containing sigma_rv, snr, if desired
IF KEYWORD_SET(SAVE_FILE) THEN BEGIN
	save_file_name_tags = 'dv_photon_' + teff_str + '_' + v_mag_str + '_' + exp_time_str
	
	;save ascii file with mean order SNR, dv estimate
	OPENW,UnitW,save_file_name_tags + '.txt',/Get_LUN
	PRINTF,UnitW,'Teff: ' + teff_str
	PRINTF,UnitW,'Vmag: ' + v_mag_str
	PRINTF,UnitW,'t_exposure: ' + exp_time_str
	PRINTF,UnitW,'vsini: ' + vsini_str
	PRINTF,UnitW,' '
	PRINTF,UnitW,'--------------------------'
	PRINTF,UnitW,''
	PRINTF,UnitW,'Total velocity uncertainty: ' + STRCOMPRESS(STRING(sigma_rv_val,format='(f20.2)'),/REMOVE) + ' m/s'
	PRINTF,UnitW,''

	PRINTF,UnitW,'Wavelength [nm]','	Mean SNR per pixel','	Sigma_rv [m/s]'
	FOR i=0,N_ELEMENTS(order_arr)-1 DO BEGIN
	     PRINTF,UnitW,wvl_arr[i],snr_rv_ord[i],sigma_rv_ord[i],format='(f20.3,f20.3,f20.3)'
	ENDFOR
	FREE_LUN,UnitW
ENDIF

RETURN,sigma_rv_val

END