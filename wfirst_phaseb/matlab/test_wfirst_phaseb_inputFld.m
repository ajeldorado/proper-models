%  3/14/2019, Hanying Zhou  
clear all
% addpath('/home/hanying/afta_im/properMatlab2/proper_v3.0c_matlab_4jun18')
addpath('~/Documents/MATLAB/PROPER');
addpath(genpath('/Users/ajriggs/Documents/Sim/cgi/wfirst_phaseb/'))
cd /Users/ajriggs/Documents/Sim/cgi/wfirst_phaseb/

testcase = 3;%1;  % 1 -hlc; 2- hlc_erkin; 3 -spc-ifs (long); 4- ifs_short;  5-spc_wide; 

%-------------------------------------------
nlam = 1;%9;  
npsf = 256;
bw = 0.10;                          % 0.18 for ifs
optval.final_sampling_lam0 = 0.1;   % 0.2 for disc
optval.source_x_offset = 7;         % 12 for disc

switch testcase
    case 1      
        lambda0 =.575; 
        optval.cor_type ='hlc'; 
    case 2      
        lambda0 =.575; 
        optval.cor_type ='hlc_erkin'; 
    case 3
        lambda0 =.730;  
        bw = 0.18;
        optval.cor_type ='spc-ifs_long';
    case 4
        lambda0 =.660;  
        bw = 0.18;
        optval.cor_type ='spc-ifs_short';
    case 5
        lambda0 =.825;
        optval.cor_type ='spc-wide'; 
        optval.final_sampling_lam0 = 0.2; 
        optval.source_x_offset = 12;
end
        
lam_array = [(1-bw/2):bw/(nlam-mod(nlam,2)):(1+bw/2)]*lambda0; 
if nlam ==1  bw = 0; lam_array =lambda0; end

%%
% % 1. offaxis, compact  psf pk for normalizing factor; 
% % EE0 = prop_run_multi(['wfirst_phaseb_v2_compact' ], lam_array, npsf, 'quiet', 'passvalue',optval );
% EE0 = prop_run(['wfirst_phaseb_v2_compact' ], lam_array, npsf, 'quiet', 'passvalue',optval );

optval.source_x_offset =0;
optval.zindex = 4;
optval.zval_m = 0.19e-9;
optval.use_errors =1;
optval.polaxis = 10; %-2; 
% 
% % 2. full model, for regular psf 
% % EE = prop_run_multi(['wfirst_phaseb_v2'], lam_array, npsf, 'quiet', 'passvalue',optval );
% EE = prop_run(['wfirst_phaseb_v2'], lam_array, npsf, 'quiet', 'passvalue',optval );


% 3. full model, for field

optval.dm1_m = fitsread (['dm1_flatten.fits']);
optval.use_dm1 =1;

optval.end_at_fpm_exit_pupil =1;
optval.output_field_rootname = ['fld_at_xtPup'];
optval.use_fpm=0;
optval.use_hlc_dm_patterns=0;
nout = 512; 			% nout > pupil_daim_pix
if testcase >2; nout =1024; end % >= pupil_daim_pix; 
% fld = prop_run_multi(['wfirst_phaseb_v2'], lam_array, nout, 'quiet', 'passvalue',optval );
fld = prop_run(['wfirst_phaseb_v2'], lam_array, nout, 'quiet', 'passvalue',optval );

figure(601); imagesc(angle(fld)); axis xy equal tight; colorbar; colormap hsv;
figure(602); imagesc(abs(fld)); axis xy equal tight; colorbar; colormap parula;
% figure(603); imagesc(angle(fld)); axis xy equal tight; colorbar; colormap hsv;
% figure(604); imagesc(abs(fld)); axis xy equal tight; colorbar; colormap parula;

% figure(605); imagesc(angle(fld)); axis xy equal tight; colorbar; colormap hsv;
% figure(606); imagesc(abs(fld)); axis xy equal tight; colorbar; colormap parula;

%%
return
%%

% 4. compact, for psf; use field input from full

optval.use_fpm= 1;
optval.use_hlc_dm_patterns= 1;
optval.input_field_rootname = ['fld_at_xtPup'];
optval.zindex = 0;
optval.zval_m = 0;
    
Ee= prop_run_multi(['wfirst_phaseb_v2_compact'], lam_array, npsf, 'quiet', 'passvalue',optval ); 
    

max_psf0 = max(max(sum(abs(EE0).^2,3)))/nlam;
psf_full = sum(abs(EE).^2,3)/nlam/max_psf0;
psf_cp   = sum(abs(Ee).^2,3)/nlam/max_psf0;
    
figure(1), clf,imagesc([log10(abs(psf_full)) log10(abs(psf_cp))]);
colorbar, caxis([-7 -2]),axis image %
title([ strrep(optval.cor_type,'_',', ') ', ' num2str(nlam) 'lam, full vs compact w/ field input from full'])

    
print('-dpng', [optval.cor_type  '_full_vs_cp_fldFrmFll_NI']);

