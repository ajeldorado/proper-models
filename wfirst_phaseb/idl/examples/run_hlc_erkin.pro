;   Copyright 2019 California Institute of Technology
; ------------------------------------------------------------------

nlam = 7
bandwidth = 0.1d	;-- bandpass fractional width
lam0 = 0.575d		;-- central wavelength (microns)
minlam = lam0 * (1 - bandwidth/2)
maxlam = lam0 * (1 + bandwidth/2)
lam_array = dindgen(nlam) / (nlam-1) * (maxlam - minlam) + minlam

final_sampling = 0.1
npsf = 512

print, 'Computing unaberrated coronagraphic field using default DM wavefront maps'

optval = {cor_type:'hlc_erkin', use_hlc_dm_patterns:1, zindex:4, zval_m:0.19e-9, use_errors:0, $
	  final_sampling_lam0:final_sampling}
prop_run_multi, 'wfirst_phaseb', fields, lam_array, npsf, /quiet, PASSVALUE=optval
images = abs(fields)^2
image = total(images,3) / nlam

print, 'Computing 7 lam/D offset source using compact model and default DM wavefront maps'

optval = {cor_type:'hlc_erkin', use_hlc_dm_patterns:1, final_sampling_lam0:final_sampling, source_x_offset:7.0}
prop_run_multi, 'wfirst_phaseb_compact', fields, lam_array, npsf, /quiet, PASSVALUE=optval
images = abs(fields)^2
psf = total(images,3) / nlam
max_psf = max(psf)

ni = image / max_psf

window, xs=400, ys=400

showcontrast, ni, final_sampling, 3, 9, /circ, /mat, grid=400, min=1e-10, max=1e-7, mag=1.3

end
