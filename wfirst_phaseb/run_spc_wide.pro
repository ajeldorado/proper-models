lam0 = 0.825d
lam_min = lam0 * 0.95d
lam_max = lam0 * 1.05d
nlam = 9
lam_array = dindgen(nlam) / (nlam - 1) * (lam_max - lam_min) + lam_min

final_sampling = 0.1

;-- compute coronagraphic image

prop_run_multi, 'wfirst_phaseb', fields, lam_array, 512, /quiet, PASSVALUE={cor_type:'spc-wide',lam0:lam0,use_errors:0,zindex:4,zval_m:0.19e-9,final_sampling_lam0:final_sampling}
images = abs(fields)^2
image = total(images,3) / nlam

;-- move source by 7.0 lam/D

prop_run_multi, 'wfirst_phaseb_compact', fields, lam_array, 512, /quiet, PASSVALUE={cor_type:'spc-wide',lam0:lam0,source_x_offset:12.0,final_sampling_lam0:final_sampling}
psfs = abs(fields)^2
psf = total(psfs,3) / nlam

ni = image / max(psf)

showcontrast, ni, final_sampling, 5.4, 20, /circ, /mat, min=1e-10, max=1e-7

end
