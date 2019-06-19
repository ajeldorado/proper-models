%   Copyright 2019 California Institute of Technology
% ------------------------------------------------------------------

clear all
delete(gcp('nocreate'))

addpath(genpath('~/Downloads/wfirst_phaseb_v1.2/matlab/'))
optval.data_dir = '/Users/ajriggs/Downloads/phaseb_data/';

n = 256;

lambda0 = 0.575; 
nlam = 1;
bw = 0.01;
% lam_array = [(1-bw/2):bw/(nlam-mod(nlam,2)):(1+bw/2)]*lambda0; 
lam_array = lambda0;
% compute coronagraphic field with full prescription

optval.final_sampling_lam0 = 0.1;
optval.zindex = 4;
optval.zval_m = 0.19e-9;
optval.use_errors = 1;
optval.polaxis = -2;
optval.use_hlc_dm_patterns = 1;
EE = prop_run(['wfirst_phaseb'], lam_array, n, 'quiet', 'passvalue',optval );
image_full = sum(abs(EE).^2,3) / nlam;

% write FPM exit pupil fields with full prescription

optval.end_at_fpm_exit_pupil = 1;
optval.output_field_rootname = ['hlc_input_field'];
optval.use_hlc_dm_patterns = 0;
optval.use_fpm = 0;
fld = prop_run(['wfirst_phaseb'], lam_array, n, 'quiet', 'passvalue',optval );

% compute coronagraphic field with compact prescription and input fields

clear optval
optval.data_dir = '/Users/ajriggs/Downloads/phaseb_data/';


optval.final_sampling_lam0 = 0.1;
optval.polaxis = -2;
optval.use_hlc_dm_patterns = 1;
optval.input_field_rootname = ['hlc_input_field'];
Ee = prop_run(['wfirst_phaseb_compact'], lam_array, n, 'quiet', 'passvalue',optval ); 
image = sum(abs(Ee).^2,3) / nlam;

%  move source to 7 lam/D to get PSF peak for normalization

clear optval
optval.data_dir = '/Users/ajriggs/Downloads/phaseb_data/';

optval.final_sampling_lam0 = 0.1;
optval.use_hlc_dm_patterns = 1;
optval.source_x_offset = 7;         
EE0 = prop_run(['wfirst_phaseb_compact'], lam_array, n, 'quiet', 'passvalue',optval );
max_psf0 = max(max(sum(abs(EE0).^2,3))) / nlam;

ni_full = image_full / max_psf0;
ni = image / max_psf0;
   
figure(1), clf, imagesc( [log10(ni_full) log10(ni)] );
colorbar, caxis([-7 -2]),axis image 
colormap jet
title([ num2str(nlam) 'lam, full vs compact w/ field input from full'])

return

