%  3/14/2019, Hanying Zhou  
clear all
% addpath('/home/hanying/afta_im/properMatlab2/proper_v3.0c_matlab_4jun18')
addpath('~/Documents/MATLAB/PROPER');
addpath(genpath('/Users/ajriggs/Documents/Sim/cgi/wfirst_phaseb/'))
cd /Users/ajriggs/Documents/Sim/cgi/wfirst_phaseb/

delete(gcp('nocreate'))

testcase = 1;%3;%1;  % 1 -hlc; 2- hlc_erkin; 3 -spc-ifs (long); 4- ifs_short;  5-spc_wide; 

%-------------------------------------------
nlam = 1;%9;  
npsf = 256;
bw = 0.10;                          % 0.18 for ifs
optval.final_sampling_lam0 = 0.1;   % 0.2 for disc
optval.source_x_offset = 7;         % 12 for disc

switch testcase
    case 1      
        lambda0 =.575; 
%         lambda0 = .730; %--For generating DM flat maps
%         lambda0 = .825; %--For generating DM flat maps
%         lambda0 = .660; %--For generating DM flat maps
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

% 1. offaxis compact psf pk for normalizing factor; 
% EE0 = prop_run_multi(['wfirst_phaseb_v2_compact'  ], lam_array, npsf, 'quiet', 'passvalue',optval );
EE0 = prop_run(['wfirst_phaseb_v2_compact'  ], lambda0, npsf, 'quiet', 'passvalue',optval );

optval.source_x_offset =0;
optval.use_errors =1;

%--DEBUGGING
optval.polaxis = 10; %--DEBUGGING -----------------------------------------

% 2. before flatten, full
% EE = prop_run_multi(['wfirst_phaseb_v2' ], lam_array, npsf, 'quiet', 'passvalue',optval );
EE = prop_run(['wfirst_phaseb_v2' ], lambda0, npsf, 'quiet', 'passvalue',optval );


% 3. post flatten, full
% optval.dm1_m = fitsread(['dm1_flatten_pol10.fits']);
optval.dm1_m = 0;
% optval.dm1_m = fitsread(['dm1_flatten.fits']);
% optval.dm1_m = fitsread (['/home/hanying/afta_im/idl_etc/dm1_flatten.fits']);
optval.use_dm1 = 1;
% Ef = prop_run_multi(['wfirst_phaseb_v2' ], lam_array, npsf, 'quiet', 'passvalue',optval );
Ef = prop_run(['wfirst_phaseb_v2' ], lambda0, npsf, 'quiet', 'passvalue',optval );

max_psf0 = max(max(sum(abs(EE0).^2,3)))/nlam;
psf_bfF  = sum(abs(EE).^2,3)/nlam/max_psf0;
psf_afF  = sum(abs(Ef).^2,3)/nlam/max_psf0;

figure(1), clf,imagesc([log10(abs(psf_bfF)) log10(abs(psf_afF))]);%
colorbar, caxis([-7 -2]), axis image%
title([ strrep(optval.cor_type,'_',', ') ', ' num2str(nlam) 'lam, before and after flattening' ])

% print('-dpng', [optval.cor_type '_bf_vs_aftr_flatten_NI']);

%%
% return
%%
%--DEBUGGING

Nbeam = 309;
Narray = 310;

optval.output_dim = Narray;
optval.use_fpm = 0;
optval.use_hlc_dm_patterns= 0;
optval.end_at_fpm_exit_pupil =1;
optval.output_field_rootname = ['fld_at_xtPup'];
Epup = prop_run(['wfirst_phaseb_v2b' ], lam_array, npsf, 'quiet', 'passvalue',optval );
figure(110); imagesc(angle(Epup)); axis xy equal tight; colorbar;

mask = abs(Epup);
mask = mask/max(mask(:));
mask(mask<1/4) = 0;
mask(mask~=0) = 1;
mask = logical(mask);
wfe = angle(Epup).*mask;

% lambda0 = 575; % nm
wfe = wfe/(2*pi/(lambda0*1e3));

figure(111); imagesc(wfe); axis xy equal tight; colorbar;
% figure(111); imagesc(angle(Epup),0.3*[-1 1]); axis xy equal tight; colorbar;



dm.inf0 = fitsread('influence_dm5v2.fits');
dm.dm_spacing = 0.9906e-3;
dm.dx_inf0 = dm.dm_spacing/10;

xs = (-Narray/2:Narray/2-1)/Nbeam;
[X,Y] = meshgrid(xs);
RS = sqrt(X.^2+Y.^2);
% figure(301); imagesc(RS*Nbeam); axis xy equal tight; colorbar;


%%
Nbeam2 = 46.3;
Narray2 = 48;
xs2 = cosd(5.7)*(-(Narray2-1)/2:(Narray2-1)/2)/Nbeam2;
ys2 = (-(Narray2-1)/2:(Narray2-1)/2)/Nbeam2;
[X2,Y2] = meshgrid(xs2,ys2);

% figure(121); imagesc(X); axis xy equal tight; colorbar;
% figure(122); imagesc(Y); axis xy equal tight; colorbar;
% figure(123); imagesc(X2); axis xy equal tight; colorbar;
% figure(124); imagesc(Y2); axis xy equal tight; colorbar;

% unmasked = ones(size(mask));
% unmasked(logical(mask)) = 0;
% figure(201); imagesc(unmasked); axis xy equal tight; colorbar;


% %--Boxcar Smoothing
% Npad = 330;
% wfePad = padOrCropEven(wfe,Npad);
% maskPad = logical(padOrCropEven(mask,Npad));
% wfePad0 = wfePad;
% Nboxcar = 7;
% for ii=1:100
%     for jj=1:Npad
%         wfePad(:,jj) = smooth(wfePad(:,jj),Nboxcar);
%     end
%     wfePad(maskPad) = wfePad0(maskPad);
%     for jj=1:Npad
%         wfePad(jj,:) = smooth(wfePad(jj,:),Nboxcar);
%     end
%     
%     wfePad(maskPad) = wfePad0(maskPad);
% end
% wfe1 = padOrCropEven(wfePad,310);



%--Convolution smoothing
Npad = 350;
xspad = (-Npad/2:Npad/2-1)/Nbeam;
[Xpad,Ypad] = meshgrid(xspad);
RSpad = sqrt(Xpad.^2+Ypad.^2);
THpad = atan2(Ypad,Xpad);

Rconv = 2;%1.5;%3.5;%3.5; 
conkern = zeros(size(RS));
conkern(RS*Nbeam<=Rconv) = 1;
conkern = conkern/sum(conkern(:));
figure(301); imagesc(conkern); axis xy equal tight; colorbar; drawnow;

wfePad = padOrCropEven(wfe,Npad);
maskPad = logical(padOrCropEven(mask,Npad));
wfePad0 = wfePad;
figure(303); imagesc(maskPad); axis xy equal tight; colorbar; drawnow;

% unmaskPad = ones(size(maskPad));
% unmaskPad(logical(maskPad)) = 0;
% figure(201); imagesc(unmaskPad); axis xy equal tight; colorbar;


% wfe1 = ifftshift( ifft2(fft2(fftshift(conkern)).*fft2(fftshift(wfe))) );
Nitr = 1000;
for ii=1:Nitr
    if(ii==Nitr)
        Rconv = 3.5;
        conkern = zeros(size(RS));
        conkern(RS*Nbeam<=Rconv) = 1;
        conkern = conkern/sum(conkern(:));
    end
    wfePad(maskPad) = wfePad0(maskPad);
    wfePad = ifftshift( ifft2(fft2(fftshift(padOrCropEven(conkern,Npad))).*fft2(fftshift(padOrCropEven(wfePad,Npad)))) );

end
wfe1 = padOrCropEven(wfePad,310); 


% %--Sort
% % inds = find(wfePad0==0 & RSpad>=0.49);
% 
% wfe1Pad = wfePad;
% wfe1Pad(RSpad>=0.50) = 0;
% figure(221); imagesc(wfe1Pad); axis xy equal tight; colorbar;  drawnow;
% 
% inds = find(wfe1Pad==0);
% 
% 
% RSsubset = RSpad(inds);
% [RSsubsetSort,inds2] = sort(RSsubset);
% THsubset = THpad(inds);
% THsubsetSort = THsubset(inds2);
% 


%




% wfe1 = fftshift( fft2(fftshift(padOrCropEven(conkern,Npad))).*fft2(fftshift(padOrCropEven(wfe,Npad))) );
% wfe1 = real(padOrCropEven(wfe1,310));
% figure(302); imagesc(imag(wfe1)); axis xy equal tight; colorbar;  drawnow;

% wfeLowRes = interp2(X,Y,wfe1,X2,Y2,'cubic',0);
% wfeLowRes = padOrCropEven(wfeLowRes,48);


wfeLowRes = interp2(Xpad,Ypad,wfePad,X2,Y2,'cubic',0);



figure(201); imagesc(wfe); axis xy equal tight; colorbar;  drawnow;
figure(202); imagesc(wfe1); axis xy equal tight; colorbar;  drawnow;

figure(112); imagesc(wfeLowRes); axis xy equal tight; colorbar;  drawnow;

% %--Influence function resampled to actuator map resolution
% actres2 = 1; %--pixels per actuator width
% N2 = ceil_even(N1*actres2/actres1)+1; %--Make odd to have peak of 1
% xq = (-(N2-1)/2:(N2-1)/2)/actres2;
% [Xq,Yq] = meshgrid(xq);
% inf2 = interp2(X,Y,inf1,Xq,Yq,'cubic',0);



Vout = -1e-9*1/2*falco_fit_dm_surf(dm,wfeLowRes);
fitswrite(Vout,sprintf('dm1_flatten_pol10_%dnm.fits',round(1e3*lambda0)));

figure(113); imagesc(Vout); axis xy equal tight; colorbar; drawnow;



%--Try out the flat map
optval.dm1_m = Vout;
optval.use_dm1 =1;

% optval.output_dim = Narray;
% optval.use_fpm = 0;
% optval.use_hlc_dm_patterns= 0;
% optval.end_at_fpm_exit_pupil =1;
% optval.output_field_rootname = ['fld_at_xtPup'];
Epup2 = prop_run(['wfirst_phaseb_v2b' ], lam_array, npsf, 'quiet', 'passvalue',optval );
figure(120); imagesc(angle(Epup2)); axis xy equal tight; colorbar; drawnow;
figure(121); imagesc(angle(Epup2).*mask); axis xy equal tight; colorbar; drawnow;
figure(122); imagesc(angle(Epup2).*mask,0.2*[-1 1]); axis xy equal tight; colorbar; drawnow;

% lambda0 = 575; % nm
wfe2 = angle(Epup2).*mask;
wfe2 = wfe2/(2*pi/(lambda0*1e3));

figure(211); imagesc(wfe2,20*[-1 1]); axis xy equal tight; colorbar; drawnow;

% function Vout = falco_fit_dm_surf(dm,Vin)
% 
% %--Starting influence function
% inf1 = dm.inf0;
% N1 = length(inf1);
% actres1 = dm.dm_spacing/dm.dx_inf0;
% x = (-(N1-1)/2:(N1-1)/2)/actres1;
% [X,Y] = meshgrid(x);
% 
% %--Influence function resampled to actuator map resolution
% actres2 = 1; %--pixels per actuator width
% N2 = ceil_even(N1*actres2/actres1)+1; %--Make odd to have peak of 1
% xq = (-(N2-1)/2:(N2-1)/2)/actres2;
% [Xq,Yq] = meshgrid(xq);
% inf2 = interp2(X,Y,inf1,Xq,Yq,'cubic',0);
% 
% %--Perform the fit
% [Vout, ~] = prop_fit_dm(Vin, inf2);
% 
% end %--END OF FUNCTION
%%

% % 3. post flatten, full
% % clear optval
% % optval = rmfield(optval,'output_field_rootname');
% optval.output_dim = 256;
% optval.use_fpm = 1;
% optval.use_hlc_dm_patterns= 1;
% optval.end_at_fpm_exit_pupil = 0;
% % optval.dm1_m = fitsread(['dm1_flatten.fits']);
% optval.final_sampling_lam0 = 0.1;   % 0.2 for disc
% optval.source_x_offset = 7;         % 12 for disc
% 
% % optval.dm1_m = fitsread (['/home/hanying/afta_im/idl_etc/dm1_flatten.fits']);
% optval.use_dm1 = 1;
% % Ef = prop_run_multi(['wfirst_phaseb_v2' ], lam_array, npsf, 'quiet', 'passvalue',optval );
% Ef = prop_run(['wfirst_phaseb_v2' ], lambda0, npsf, 'quiet', 'passvalue',optval );
% 
% max_psf0 = max(max(sum(abs(EE0).^2,3)))/nlam;
% psf_bfF  = sum(abs(EE).^2,3)/nlam/max_psf0;
% psf_afF  = sum(abs(Ef).^2,3)/nlam/max_psf0;
% 
% 
% figure(90); imagesc(log10(abs(psf_afF))); colorbar; axis xy equal tight;
% 
% figure(91), clf,imagesc([log10(abs(psf_bfF)) log10(abs(psf_afF))]);%
% colorbar, caxis([-7 -2]), axis image%
% title([ strrep(optval.cor_type,'_',', ') ', ' num2str(nlam) 'lam, before and after flattening' ])
% 
% 
