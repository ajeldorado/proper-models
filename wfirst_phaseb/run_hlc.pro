lam_array = [ 0.54625d, 0.5534375d, 0.560625d, 0.5678125d, 0.575d, 0.5821875d, 0.589375d, 0.5965625d, 0.60375d ]
lam0 = 0.575d
nlam = n_elements(lam_array)

final_sampling = 0.1

;-- compute coronagraphic image

prop_run_multi, 'wfirst_phaseb', fields, lam_array, 512, /quiet, PASSVALUE={cor_type:'hlc',zindex:4,zval_m:0.19e-9,use_errors:0,final_sampling_lam0:final_sampling}
images = abs(fields)^2
image = total(images,3) / nlam

;-- move source by 7.0 lam/D

prop_run_multi, 'wfirst_phaseb_compact', fields, lam_array, 512, /quiet, PASSVALUE={cor_type:'hlc',source_x_offset:7.0,final_sampling_lam0:final_sampling}
psfs = abs(fields)^2
psf = total(psfs,3) / nlam

ni = image / max(psf)

window, xs=400, ys=400
showcontrast, ni, final_sampling, 3, 9, /circ, /mat, grid=400, min=1e-10, max=1e-7

end
