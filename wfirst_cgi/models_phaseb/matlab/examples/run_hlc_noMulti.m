%   Copyright 2019 California Institute of Technology
% ------------------------------------------------------------------

clear all
% delete(gcp('nocreate'))

data_dir = '/Users/ajriggs/Repos/proper-models/wfirst_cgi/data_phaseb/';

addpath /Users/ajriggs/Repos/proper-models/wfirst_cgi/models_phaseb/matlab/

optval.data_dir = data_dir;

npsf = 256;

lambda0 = 0.575; 
nlam = 1;%7;  
bw = 0.10;
lam_array = lambda0;%[(1-bw/2):bw/(nlam-mod(nlam,2)):(1+bw/2)]*lambda0;

% compute unberrated coronagraphic field using Dwight's DM wavefront mapsg 

optval.cor_type ='hlc'; 
optval.final_sampling_lam0 = 0.1;
optval.use_errors = 0;
optval.use_hlc_dm_patterns = 1;
optval.zindex = 4;
optval.zval_m = 0.19e-9;


fields = prop_run( ['model_full_wfirst_phaseb'], lam_array, npsf, 'quiet', 'passvalue',optval );
% fields = prop_run( ['wfirst_phaseb'], lam_array, npsf, 'quiet', 'passvalue',optval );
image  = sum(abs(fields).^2,3) / nlam;

% compute unberrated coronagraphic field using DM actuator settings

optval.dm1_m = fitsread(['hlc_dm1.fits']);
optval.dm2_m = fitsread(['hlc_dm2.fits']);
optval.use_dm1 = 1;
optval.use_dm2 = 1;
optval.use_hlc_dm_patterns = 0;

fields = prop_run( ['model_full_wfirst_phaseb'], lam_array, npsf, 'quiet', 'passvalue',optval );
% fields = prop_run( ['wfirst_phaseb'], lam_array, npsf, 'quiet', 'passvalue',optval );
image_dm  = sum(abs(fields).^2,3) / nlam;

% compute aberrated coronagraphic field using DM actuator settings

optval.dm1_m = fitsread(['hlc_with_aberrations_dm1.fits']);
optval.dm2_m = fitsread(['hlc_with_aberrations_dm2.fits']);
optval.use_errors = 1;
optval.polaxis = 10;


fields = prop_run( ['model_full_wfirst_phaseb'], lam_array, npsf, 'quiet', 'passvalue',optval );
% fields = prop_run( ['wfirst_phaseb'], lam_array, npsf, 'quiet', 'passvalue',optval );
image_ab  = sum(abs(fields).^2,3) / nlam;

% offaxis compact psf pk for normalizing factor (no errors) 

clear optval
optval.data_dir = data_dir;
optval.cor_type ='hlc'; 
optval.final_sampling_lam0 = 0.1;
optval.use_hlc_dm_patterns = 1;
optval.source_x_offset = 7.0;


fields = prop_run( ['model_compact_wfirst_phaseb'], lam_array, npsf, 'quiet', 'passvalue',optval );
% fields = prop_run( ['wfirst_phaseb_compact'], lam_array, npsf, 'quiet', 'passvalue',optval );
psf  = sum(abs(fields).^2,3) / nlam;
max_psf = max(max(psf));

ni = image / max_psf;
ni_dm = image_dm / max_psf;
ni_ab = image_ab / max_psf;

figure(1), clf, imagesc([log10(ni) log10(ni_dm), log10(ni_ab)]);
colorbar, caxis([-10 -7]), axis image
title([ 'Unaberrated+DM map, unaberrated+DM pistons, aberrated+DM pistons' ])
colormap jet

return
