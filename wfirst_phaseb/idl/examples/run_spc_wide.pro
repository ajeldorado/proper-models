lam0 = 0.825d
lam_min = lam0 * 0.95d
lam_max = lam0 * 1.05d
nlam = 9
lam_array = dindgen(nlam) / (nlam - 1) * (lam_max - lam_min) + lam_min

final_sampling = 0.15
npsf = 512

print, 'Computing coronagraphic image with full model'

optval = {cor_type:'spc-wide', use_errors:0, zindex:4, zval_m:0.19e-9, final_sampling_lam0:final_sampling}

prop_run_multi, 'wfirst_phaseb', fields, lam_array, npsf, /quiet, PASSVALUE=optval
images = abs(fields)^2
image = total(images,3) / nlam

print, 'Computing offset source with compact model'

optval = {cor_type:'spc-wide', source_x_offset:7.0, final_sampling_lam0:final_sampling}

prop_run_multi, 'wfirst_phaseb_compact', fields, lam_array, npsf, /quiet, PASSVALUE=optval
psfs = abs(fields)^2
psf = total(psfs,3) / nlam

ni = image / max(psf)

window, xs=512, ys=512
showcontrast, ni, final_sampling, 5.4, 20, /circ, mag=1.3, /mat, min=1e-10, max=1e-7

end
