;   Copyright 2019 California Institute of Technology
; ------------------------------------------------------------------

nlam = 7
bandwidth = 0.1d	;-- bandpass fractional width
lam0 = 0.575d		;-- central wavelength (microns)
minlam = lam0 * (1 - bandwidth/2)
maxlam = lam0 * (1 + bandwidth/2)
lam_array = dindgen(nlam) / (nlam-1) * (maxlam - minlam) + minlam

final_sampling = 0.1	;-- image sampling is 0.1 lam0/D
npsf = 512

print, 'Computing field before flattening' 

optval = {cor_type:'hlc', use_errors:1, polaxis:10, use_hlc_dm_patterns:1, $
	  final_sampling_lam0:final_sampling}

prop_run_multi, 'wfirst_phaseb', fields, lam_array, npsf, /quiet, PASSVALUE=optval 
images = abs(fields)^2
image_before = total(images,3) / nlam

print, 'Computing field after flattening' 

fits_read, 'errors_polaxis10_dm.fits', dm1

optval = {cor_type:'hlc', use_errors:1, polaxis:10, use_hlc_dm_patterns:1, $
	  final_sampling_lam0:final_sampling, use_dm1:1, dm1_m:dm1}

prop_run_multi, 'wfirst_phaseb', fields, lam_array, npsf, /quiet, PASSVALUE=optval
images = abs(fields)^2
image_after = total(images,3) / nlam

print, 'Computing 7 lam/D offset source using compact model'

optval = {cor_type:'hlc', use_hlc_dm_patterns:1, final_sampling_lam0:final_sampling, $
	  source_x_offset:7.0}

prop_run_multi, 'wfirst_phaseb_compact', fields, lam_array, npsf, /quiet, PASSVALUE=optval
images = abs(fields)^2
psf = total(images,3) / nlam
max_psf = max(psf)

ni_before = image_before / max_psf
ni_after = image_after / max_psf

window, xs=800, ys=400
showcontrast, ni_before, /mat, grid=400, min=1e-7, max=1e-2, mag=1.3
xyouts, 200, 380, 'Before', charsize=1.5, align=0.5, /dev
showcontrast, ni_after, /mat, grid=400, min=1e-7, max=1e-2, xoff=400, mag=1.3
xyouts, 600, 380, 'After', charsize=1.5, align=0.5, /dev

end
