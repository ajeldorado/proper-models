lam_array = [ 0.54625d, 0.5534375d, 0.560625d, 0.5678125d, 0.575d, 0.5821875d, 0.589375d, 0.5965625d, 0.60375d ]
nlam = n_elements(lam_array)

final_sampling = 0.1

;-- compute field before flattening 

prop_run_multi, 'wfirst_phaseb', fields, lam_array, 512, /quiet, $
	PASSVALUE={cor_type:'hlc',final_sampling_lam0:final_sampling}
images = abs(fields)^2
image_before = total(images,3) / nlam

;-- compute field after flattening 

fits_read, 'dm1_flatten.fits', dm1

prop_run_multi, 'wfirst_phaseb', fields, lam_array, 512, /quiet, $
	PASSVALUE={cor_type:'hlc',use_dm1:1,dm1_m:dm1,final_sampling_lam0:final_sampling}
images = abs(fields)^2
image_after = total(images,3) / nlam

;-- move source by 7.0 lam/D

prop_run_multi, 'wfirst_phaseb_compact', fields, lam_array, 512, /quiet, $
	PASSVALUE={cor_type:'hlc',source_x_offset:7.0,final_sampling_lam0:final_sampling}
psfs = abs(fields)^2
psf = total(psfs,3) / nlam

ni_before = image_before / max(psf)
ni_after = image_after / max(psf)

window, xs=800, ys=400
showcontrast, ni_before, /mat, grid=400, min=1e-7, max=1e-2
xyouts, 200, 380, 'Before', charsize=1.5, align=0.5, /dev
showcontrast, ni_after, /mat, grid=400, min=1e-7, max=1e-2, xoff=400
xyouts, 600, 380, 'After', charsize=1.5, align=0.5, /dev

end
