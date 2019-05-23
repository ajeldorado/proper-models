lam_array = [ 0.54625d, 0.5534375d, 0.560625d, 0.5678125d, 0.575d, 0.5821875d, 0.589375d, 0.5965625d, 0.60375d ]
lam0 = 0.575d
nlam = n_elements(lam_array)

final_sampling = 0.1

;-- write FPM exit pupil fields with full prescription

prop_run_multi, 'wfirst_phaseb', fields, lam_array, 512, /quiet, $
	PASSVALUE={cor_type:'hlc',polaxis:-2,zindex:4,zval_m:0.19e-9,use_errors:1,$
		   final_sampling_lam0:final_sampling,use_hlc_dm_patterns:0,use_fpm:0,$
		   end_at_fpm_exit_pupil:1,output_field_rootname:'hlc_input_field'}

;-- compute coronagraphic image with full prescription
t1 = systime(/sec)
prop_run_multi, 'wfirst_phaseb', fields, lam_array, 512, /quiet, $
	PASSVALUE={cor_type:'hlc',polaxis:-2,zindex:4,zval_m:0.19e-9,use_errors:1,$
		   final_sampling_lam0:final_sampling}
print, systime(/sec) - t1
images = abs(fields)^2
image_full = total(images,3) / nlam

;-- compute coronagraphic image with compact prescription and pre-computed input fields

t1 = systime(/sec)
prop_run_multi, 'wfirst_phaseb_compact', fields, lam_array, 512, /quiet, $
	PASSVALUE={cor_type:'hlc',polaxis:-2,final_sampling_lam0:final_sampling,$
		   input_field_rootname:'hlc_input_field'}
print, systime(/sec) - t1
images = abs(fields)^2
image = total(images,3) / nlam

;-- move source by 7.0 lam/D

prop_run_multi, 'wfirst_phaseb_compact', fields, lam_array, 512, /quiet, $
	PASSVALUE={cor_type:'hlc',source_x_offset:7.0,final_sampling_lam0:final_sampling}
psfs = abs(fields)^2
psf = total(psfs,3) / nlam

ni_full = image_full / max(psf)
ni = image / max(psf)

window, xs=800, ys=400
showcontrast, ni_full, final_sampling, 3, 9, /circ, /mat, grid=400, min=1e-7, max=1e-2
showcontrast, ni, final_sampling, 3, 9, /circ, /mat, grid=400, min=1e-7, max=1e-2, xoff=400

end
