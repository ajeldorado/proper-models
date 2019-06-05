clear all
delete(gcp('nocreate'))

npsf = 256;

lambda0 = 0.730; 
nlam = 9;  
bw = 0.18;
lam_array = [(1-bw/2):bw/(nlam-mod(nlam,2)):(1+bw/2)]*lambda0;

% compute unberrated coronagraphic field 

optval.cor_type ='spc-ifs_long'; 
optval.final_sampling_lam0 = 0.1;
optval.use_errors = 0;
optval.zindex = 4;
optval.zval_m = 0.19e-9;

fields = prop_run_multi( ['wfirst_phaseb'], lam_array, npsf, 'quiet', 'passvalue',optval );
image  = sum(abs(fields).^2,3) / nlam;

% offaxis compact psf pk for normalizing factor (no errors) 

optval.zindex = 0;
optval.zval_m = 0;
optval.source_x_offset = 7.0;

fields = prop_run_multi( ['wfirst_phaseb_compact'], lam_array, npsf, 'quiet', 'passvalue',optval );
psf  = sum(abs(fields).^2,3) / nlam;
max_psf = max(max(psf));

ni = image / max_psf;

figure(1), clf, imagesc([log10(ni)]);
colorbar, caxis([-10 -7]), axis image
colormap jet

return
