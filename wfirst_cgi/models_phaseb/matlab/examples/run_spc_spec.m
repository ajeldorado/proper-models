%   Copyright 2019 California Institute of Technology
% ------------------------------------------------------------------

clear all
delete(gcp('nocreate'))


npsf = 256;

lambda0 = 0.730; 
nlam = 9;  
bw = 0.15;
lam_array = [(1-bw/2):bw/(nlam-mod(nlam,2)):(1+bw/2)]*lambda0;

optval.cor_type ='spc-spec_long'; 
optval.final_sampling_lam0 = 0.1;

% compute unaberrated coronagraphic field using compact model

fields = prop_run_multi( ['wfirst_phaseb_compact'], lam_array, npsf, 'quiet', 'passvalue',optval );
image_noab  = sum(abs(fields).^2,3) / nlam;

% compute aberrated coronagraphic field with DM actuator settings

optval.dm1_m = fitsread(['spc-spec_long_with_aberrations_dm1.fits']);
optval.dm2_m = fitsread(['spc-spec_long_with_aberrations_dm2.fits']);
optval.use_dm1 = 1;
optval.use_dm2 = 1;
optval.use_errors = 1;
optval.polaxis = 10;
fields = prop_run_multi( ['wfirst_phaseb'], lam_array, npsf, 'quiet', 'passvalue',optval );
image_ab  = sum(abs(fields).^2,3) / nlam;

% offaxis compact psf pk for normalizing factor (no errors) 

clear optval
optval.cor_type ='spc-spec_long'; 
optval.final_sampling_lam0 = 0.1;
optval.source_x_offset = 7.0;

fields = prop_run_multi( ['wfirst_phaseb_compact'], lam_array, npsf, 'quiet', 'passvalue',optval );
psf  = sum(abs(fields).^2,3) / nlam;
max_psf = max(max(psf));

ni_noab = image_noab / max_psf;
ni_ab = image_ab / max_psf;

figure(1), clf, imagesc([log10(ni_noab) log10(ni_ab)]);
colorbar, caxis([-10 -7]), axis image
title([ 'Unaberrated, aberrated+DM pistons' ])
colormap jet

return
