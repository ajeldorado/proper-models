nlam = 7
bandwidth = 0.1d	;-- bandpass fractional width
lam0 = 0.575d		;-- central wavelength (microns)
minlam = lam0 * (1 - bandwidth/2)
maxlam = lam0 * (1 + bandwidth/2)
lam_array = dindgen(nlam) / (nlam-1) * (maxlam - minlam) + minlam

final_sampling = 0.1
npsf = 512

print, 'Computing aberrated coronagraphic image with full prescription'

optval = {cor_type:'hlc', use_errors:1, polaxis:-2, zindex:4, zval_m:0.19e-9, $ 
	  use_hlc_dm_patterns:1, final_sampling_lam0:final_sampling} 

t1 = systime(/sec)
prop_run_multi, 'wfirst_phaseb', fields, lam_array, npsf, /quiet, PASSVALUE=optval
print, 'Time for full prescription = ', systime(/sec) - t1
images = abs(fields)^2
image_full = total(images,3) / nlam

print, 'Computing FPM exit pupil fields with full prescription'

optval = {cor_type:'hlc', use_errors:1, polaxis:-2, zindex:4, zval_m:0.19e-9, $
	  use_hlc_dm_patterns:0, use_fpm:0, final_sampling_lam0:final_sampling, $
	  end_at_fpm_exit_pupil:1,output_field_rootname:'hlc_input_field'} 

prop_run_multi, 'wfirst_phaseb', fields, lam_array, npsf, /quiet, PASSVALUE=optval

print, 'Computing field with compact prescription and pre-computed input fields'

optval = {cor_type:'hlc', polaxis:-2, use_hlc_dm_patterns:1, $
	  final_sampling_lam0:final_sampling, input_field_rootname:'hlc_input_field'} 

t1 = systime(/sec)
prop_run_multi, 'wfirst_phaseb_compact', fields, lam_array, npsf, /quiet, PASSVALUE=optval
print, 'Time for compact prescription with input fields = ', systime(/sec) - t1
images = abs(fields)^2
image = total(images,3) / nlam

print, 'Computing offset source to derive NI'

optval = {cor_type:'hlc', use_hlc_dm_patterns:1, source_x_offset:7.0, $
	  final_sampling_lam0:final_sampling }
prop_run_multi, 'wfirst_phaseb_compact', fields, lam_array, npsf, /quiet, PASSVALUE=optval 
psfs = abs(fields)^2
psf = total(psfs,3) / nlam
max_psf = max(psf)

ni_full = image_full / max_psf
ni = image / max_psf 

window, xs=800, ys=400

showcontrast, ni_full, final_sampling, 3, 9, /circ, /mat, grid=400, min=1e-7, max=1e-2, mag=1.3
xyouts, 200, 380, 'Full model', charsize=1.5, align=0.5, /dev
showcontrast, ni, final_sampling, 3, 9, /circ, /mat, grid=400, min=1e-7, max=1e-2, xoff=400, mag=1.3
xyouts, 600, 380, 'Compact model + input fields', charsize=1.5, align=0.5, /dev

end
