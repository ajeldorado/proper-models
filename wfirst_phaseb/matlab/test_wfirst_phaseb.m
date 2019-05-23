%  3/14/2019, Hanying Zhou  

clear all
% addpath('/home/hanying/afta_im/properMatlab2/proper_v3.0c_matlab_4jun18')
addpath('~/Documents/MATLAB/PROPER');
addpath(genpath('/Users/ajriggs/Documents/Sim/cgi/wfirst_phaseb/'))
cd /Users/ajriggs/Documents/Sim/cgi/wfirst_phaseb/

delete(gcp('nocreate'))

testcase = 5;%1;  % 1 -hlc; 2- hlc_erkin; 3 -spc-ifs (short); 4- ifs_short;  5-spc_wide; 
isFull   = 1;  % 1 -full model; 0 -compact model

%---------------------------------------------
% nlam = 4;%9;  
npsf = 256;
% bw = 0.10;
optval.final_sampling_lam0 = 0.1;   % 0.2 for disc
optval.source_x_offset = 7;         % 12 for disc

switch testcase
    case 1      
        lambda0 =.575; 
        bw = 0.10; 
        nlam = 9;
        optval.cor_type ='hlc'; 
    case 2      
        lambda0 =.575; 
        bw = 0.10;
        nlam = 9;
        optval.cor_type ='hlc_erkin'; 
    case 3
        lambda0 =.730;  
        bw = 0.18;
        %nlam = 9;
        optval.cor_type ='spc-ifs_long';
    case 4
        lambda0 =.660;  
        bw = 0.18;
        %nlam = 9;
        optval.cor_type ='spc-ifs_short';
    case 5
        lambda0 =.825;
        bw = 0.10; 
        nlam = 4;
        optval.cor_type ='spc-wide'; 
        optval.final_sampling_lam0 = 0.2; 
        optval.source_x_offset = 12;
end
        
lam_array = ((1-bw/2):bw/(nlam-mod(nlam,2)):(1+bw/2))*lambda0; 
if(nlam ==1);  bw = 0; lam_array =lambda0; end
if isFull mdl_str ='';    else; mdl_str = '_compact'; end


% 1. offaxis compact model for normalizing factor 
EE0 = prop_run_multi(['wfirst_phaseb_v2_compact'  ], lam_array, npsf, 'quiet', 'passvalue',optval );

if isFull                   % full & compact diff
    optval.zindex = 4;
    optval.zval_m = 0.19e-9;
end

optval.source_x_offset =0;

% 2. regular psf full or compact
EE = prop_run_multi(['wfirst_phaseb_v2' mdl_str], lam_array, npsf, 'quiet', 'passvalue',optval );


max_psf0  = max(max(sum(abs(EE0).^2,3)))/nlam;
psf = sum(abs(EE).^2,3)/nlam/max_psf0;
    
figure(1), clf,imagesc(log10(abs(psf))); colorbar, caxis([-10 -7]), axis image%
title([ strrep(optval.cor_type,'_',', ') ', ' num2str(nlam) 'lam' ])
    
% print('-dpng', [optval.cor_type mdl_str '_NI']);
    

