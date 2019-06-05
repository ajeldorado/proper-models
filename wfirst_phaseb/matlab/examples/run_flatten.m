clear all
delete(gcp('nocreate'))

npsf = 256;

lambda0 = 0.575; 
nlam = 7;  
bw = 0.10;
lam_array = [(1-bw/2):bw/(nlam-mod(nlam,2)):(1+bw/2)]*lambda0;

optval.final_sampling_lam0 = 0.1;
optval.cor_type ='hlc'; 
optval.use_errors = 1;
optval.polaxis = 10;
optval.use_hlc_dm_patterns = 1;

% before flattening 

EE = prop_run_multi(['wfirst_phaseb' ], lam_array, npsf, 'quiet', 'passvalue',optval );

% after flattening

optval.dm1_m = fitsread(['errors_polaxis10_dm.fits']);
optval.use_dm1 = 1;
Ef = prop_run_multi(['wfirst_phaseb' ], lam_array, npsf, 'quiet', 'passvalue',optval );

% offaxis compact psf pk for normalizing factor (no errors) 

optval.use_errors = 0;
optval.polaxis = 0;
optval.source_x_offset = 7.0;
EE0 = prop_run_multi(['wfirst_phaseb_compact'  ], lam_array, npsf, 'quiet', 'passvalue',optval );

max_psf0 = max(max(sum(abs(EE0).^2,3))) / nlam;
psf_bfF  = sum(abs(EE).^2,3) / nlam / max_psf0;
psf_afF  = sum(abs(Ef).^2,3) / nlam / max_psf0;

figure(1), clf,imagesc([log10(abs(psf_bfF)) log10(abs(psf_afF))]);%
colorbar, caxis([-7 -2]), axis image%
title([ strrep(optval.cor_type,'_',', ') ', ' num2str(nlam) 'lam, before and after flattening' ])
colormap jet

return
