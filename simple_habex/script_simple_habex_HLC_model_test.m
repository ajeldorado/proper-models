% Copyright 2020, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
% Script to test the Habex PROPER prescription in Matlab
%
% Written by A.J. Riggs (JPL, CIT) in February, 2020
% 
% Matlab's PATH must know the locations of 
% - the Habex prescription 'habex_multi_coro'
% - the PROPER library
% - the FALCO package

clear

prescription = 'habex_multi_coro';

optval.cor_type = 'hlc';

lambda0 = 575e-9;    %--Central wavelength of the whole spectral bandpass [meters]
lambda_um = lambda0 * 1e6;
optval.lambda0_um = lambda0 * 1e6;

%--Focal planes
res = 3;
FOV = 14;
Rout = 12;
c_range = [-10, -2.5];
optval.nout = ceil_even(1 + res*(2*FOV)); %  dimensions of output in pixels (overrides output_dim0)
optval.final_sampling_lam0 = 1/res;	%   final sampling in lambda0/D
optval.use_field_stop = true;	%-- use field stop (0 = no stop)
optval.field_stop_radius = Rout;   %-- field stop radius in lam0/D

% %--Pupil Plane Resolutions
NbeamFull = 314.581;
optval.pupil_diam_pix = NbeamFull;
gridsize = 1024;%2^ceil(log2(NbeamFull)); 

optval.map_dir = '/Users/ajriggs/Documents/habex/maps/';	%-- directory containing optical surface error maps
optval.mask_dir = '/Users/ajriggs/Documents/habex/run819/';	%-- directory containing Habex HLC files
% optval.normLyotDiam = 0.95;
% optval.vortexCharge = 6;
optval.pupil_fn = 'run819_roman_pupil_rotated.fits';
optval.lyot_stop_fn = 'run819_roman_lyot_rotated.fits';
% optval.pupil_fn = 'run819_luvoir_pupil_aj.fits';
% optval.lyot_stop_fn = 'run819_luvoir_lyot_aj.fits';
optval.nlams = 9;
optval.bw = 0.10;
optval.fpm_aoi = '5.5'; % string, AOI in degrees.
optval.fpm_pol = 's';


optval.xoffset = 0;
optval.yoffset = 0;

%% No Correction
%--Phase Retrieval Pupil

optval.use_pr = true;
optval.pr_pupil_diam_pix = NbeamFull;%248;

optval.use_dm1 = false;
optval.use_errors = 1;		%-- 1 = use optical surface errors, 0 = none
optval.use_fpm = 0;		%-- use focal plane mask (0 = no FPM)
optval.use_lyot_stop = 1;	%-- use Lyot stop (0 = no stop)
optval.use_field_stop = 0;	%-- use field stop (0 = no stop)
% [Epup, sampling_m]  = habex_multi_coro(lambda_m, gridsize, optval);
[Epup, sampling_m] = prop_run(prescription, lambda_um, gridsize, 'quiet', 'passvalue', optval);

mask = 0*Epup;
mask(abs(Epup) > 0.1*max(abs(Epup(:)))) = 1;
% mask = ones(size(Epup));

figure(1); imagesc(abs(Epup)); axis xy equal tight; colorbar; title('abs(E$_{pupil}$)', 'Interpreter','Latex'); set(gca,'Fontsize',20); drawnow;
figure(2); imagesc(mask.*angle(Epup),[-1, 1]); axis xy equal tight; colorbar; title('angle(E$_{pupil}$)', 'Interpreter','Latex'); set(gca,'Fontsize',20); drawnow;

% PSF for normalization
optval.use_pr = false;
optval.use_errors = 1;		%-- 1 = use optical surface errors, 0 = none
optval.use_fpm = 0;		%-- use focal plane mask (0 = no FPM)
optval.use_lyot_stop = 1;	%-- use Lyot stop (0 = no stop)
optval.use_field_stop = 0;	%-- use field stop (0 = no stop)
% [EforNorm, sampling_m]  = habex_multi_coro(lambda_m, gridsize, optval);
[EforNorm, sampling_m] = prop_run(prescription, lambda_um, gridsize, 'quiet', 'passvalue', optval);

IforNorm = abs(EforNorm).^2;
I00 = max(IforNorm(:));
figure(3); imagesc(log10(IforNorm/I00),[-7 0]); axis xy equal tight; colorbar; title('PSF for Normalization)', 'Interpreter','Latex'); set(gca,'Fontsize',20); drawnow;

% Coronagraphic PSF
optval.use_errors = 1;		%-- 1 = use optical surface errors, 0 = none
optval.use_fpm = 1;		%-- use focal plane mask (0 = no FPM)
optval.use_lyot_stop = 1;	%-- use Lyot stop (0 = no stop)
optval.use_field_stop = 1;	%-- use field stop (0 = no stop)
% [Ecoro, sampling_m]  = habex_multi_coro(lambda_m, gridsize, optval);
[Ecoro, sampling_m] = prop_run(prescription, lambda_um, gridsize, 'quiet', 'passvalue', optval);
Inorm = abs(Ecoro).^2 / I00;
figure(4); imagesc(log10(Inorm), c_range); axis xy equal tight; colorbar;
title('Coronagraphic PSF before Flattening', 'Fontsize', 16); drawnow;


%% Corrected

optval.use_dm1 = true;
optval.dm1 = fitsread([optval.map_dir, 'flat_map_hlc.fits']);

%--Phase Retrieval Pupil
optval.use_pr = true;
optval.pr_pupil_diam_pix = NbeamFull;%248;
optval.use_errors = 1;		%-- 1 = use optical surface errors, 0 = none
optval.use_fpm = 0;		%-- use focal plane mask (0 = no FPM)
optval.use_lyot_stop = 1;	%-- use Lyot stop (0 = no stop)
optval.use_field_stop = 0;	%-- use field stop (0 = no stop)
% [Epup, sampling_m]  = habex_multi_coro(lambda_um, gridsize, optval);
[Epup, sampling_m] = prop_run(prescription, lambda_um, gridsize, 'quiet', 'passvalue', optval);

mask = 0*Epup;
mask(abs(Epup) > 0.1*max(abs(Epup(:)))) = 1;
% mask = ones(size(Epup));

figure(11); imagesc(abs(Epup)); axis xy equal tight; colorbar; title('abs(E$_{pupil}$)', 'Interpreter','Latex'); set(gca,'Fontsize',20); drawnow;
figure(12); imagesc(mask.*angle(Epup),[-1, 1]); axis xy equal tight; colorbar; title('angle(E$_{pupil}$)', 'Interpreter','Latex'); set(gca,'Fontsize',20);drawnow;

% PSF for normalization
optval.use_pr = false;
optval.use_errors = 1;		%-- 1 = use optical surface errors, 0 = none
optval.use_fpm = 0;		%-- use focal plane mask (0 = no FPM)
optval.use_lyot_stop = 1;	%-- use Lyot stop (0 = no stop)
optval.use_field_stop = 0;	%-- use field stop (0 = no stop)
% [EforNorm, sampling_m]  = habex_multi_coro(lambda_m, gridsize, optval);
[EforNorm, sampling_m] = prop_run(prescription, lambda_um, gridsize, 'quiet', 'passvalue', optval);
IforNorm = abs(EforNorm).^2;
I00 = max(IforNorm(:));
figure(13); imagesc(log10(IforNorm/I00),[-7 0]); axis xy equal tight; colorbar; title('PSF for Normalization)', 'Interpreter','Latex'); set(gca,'Fontsize',20); drawnow;

% Coronagraphic PSF
optval.use_errors = 1;		%-- 1 = use optical surface errors, 0 = none
optval.use_fpm = 1;		%-- use focal plane mask (0 = no FPM)
optval.use_lyot_stop = 1;	%-- use Lyot stop (0 = no stop)
optval.use_field_stop = 1;	%-- use field stop (0 = no stop)
% [Ecoro, sampling_m]  = habex_multi_coro(lambda_m, gridsize, optval);
[Ecoro, sampling_m] = prop_run(prescription, lambda_um, gridsize, 'quiet', 'passvalue', optval);
Inorm = abs(Ecoro).^2 / I00;
figure(14); imagesc(log10(Inorm), c_range); axis xy equal tight; colorbar;
title('Coronagraphic PSF after Flattening', 'Fontsize', 16);
drawnow;

% Coronagraphic PSF
optval.use_hlc_dm_patterns = true;
optval.dm1wfe_fn = 'run819_roman_dm1wfe.fits';
optval.dm2wfe_fn = 'run819_roman_dm2wfe.fits';
optval.use_errors = 1;		%-- 1 = use optical surface errors, 0 = none
optval.use_fpm = 1;		%-- use focal plane mask (0 = no FPM)
optval.use_lyot_stop = 1;	%-- use Lyot stop (0 = no stop)
optval.use_field_stop = 1;	%-- use field stop (0 = no stop)
% [Ecoro, sampling_m]  = habex_multi_coro(lambda_m, gridsize, optval);
[Ecoro, sampling_m] = prop_run(prescription, lambda_um, gridsize, 'quiet', 'passvalue', optval);
Inorm = abs(Ecoro).^2 / I00;
figure(15); imagesc(log10(Inorm), c_range); axis xy equal tight; colorbar;
title("Coronagraphic PSF after Flattening and with Dwight's DM Patterns", 'Fontsize', 16);
drawnow;

%%
% Coronagraphic PSF without aberrations
optval.use_dm1 = false;
optval.use_errors = false;		%-- 1 = use optical surface errors, 0 = none
optval.use_hlc_dm_patterns = true;
optval.dm1wfe_fn = 'run819_roman_dm1wfe.fits';
optval.dm2wfe_fn = 'run819_roman_dm2wfe.fits';
% optval.use_errors = 1;		%-- 1 = use optical surface errors, 0 = none
optval.use_fpm = 1;		%-- use focal plane mask (0 = no FPM)
optval.use_lyot_stop = 1;	%-- use Lyot stop (0 = no stop)
optval.use_field_stop = 1;	%-- use field stop (0 = no stop)
% [Ecoro, sampling_m]  = habex_multi_coro(lambda_m, gridsize, optval);
[Ecoro, sampling_m] = prop_run(prescription, lambda_um, gridsize, 'quiet', 'passvalue', optval);
Inorm = abs(Ecoro).^2 / I00;
figure(16); imagesc(log10(Inorm), c_range); axis xy equal tight; colorbar;
title("Coronagraphic PSF with Dwight's DM Patterns and No Errors", 'Fontsize', 16);
drawnow;
