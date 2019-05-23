%  3/14/2019, Hanying Zhou  
clear all
% addpath('/home/hanying/afta_im/properMatlab2/proper_v3.0c_matlab_4jun18')
addpath('~/Documents/MATLAB/PROPER');
addpath(genpath('/Users/ajriggs/Documents/Sim/cgi/wfirst_phaseb/'))
cd /Users/ajriggs/Documents/Sim/cgi/wfirst_phaseb/

testcase = 1;  % 1 -hlc; 2- hlc_erkin;

%--------------------------------------
nlam = 9;  	
npsf = 256; 
lambda0 =.575;
bw = 0.10;
optval.final_sampling_lam0 = 0.1;   
optval.source_x_offset = 7;         

switch testcase
    case 1      
        optval.cor_type ='hlc'; 
    case 2      
        optval.cor_type ='hlc_erkin'; 
end

lam_array = [(1-bw/2):bw/(nlam-mod(nlam,2)):(1+bw/2)]*lambda0;

    
% 1. offaxis compact psf pk for normalizing factor
EE0 = prop_run_multi(['wfirst_phaseb_v2_compact'  ], lam_array, npsf, 'quiet', 'passvalue',optval );

optval.source_x_offset =0;
optval.fpm_axis ='p';
optval.zindex = 4;
optval.zval_m = 0.19e-9;

% 2. pol p
Ep = prop_run_multi(['wfirst_phaseb_v2' ], lam_array, npsf, 'quiet', 'passvalue',optval );


% 3. pol s
optval.fpm_axis ='s';
Es = prop_run_multi(['wfirst_phaseb_v2' ], lam_array, npsf, 'quiet', 'passvalue',optval );


max_psf0 = max(max(sum(abs(EE0).^2,3)))/nlam;
psf_p  = sum(abs(Ep).^2,3)/nlam/max_psf0;
psf_s  = sum(abs(Es).^2,3)/nlam/max_psf0;

figure(1), clf,imagesc([log10(abs(psf_p)) log10(abs(psf_s))]);%
colorbar, caxis([-10 -7]), axis image%
title([ strrep(optval.cor_type,'_',', ') ', ' num2str(nlam) 'lam, p- vs s- polarization'])

print('-dpng', [optval.cor_type  '_p_vs_s_NI']);
