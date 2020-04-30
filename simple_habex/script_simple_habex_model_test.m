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
% - the Habex prescription 'habex_vortex'
% - the PROPER library
% - the FALCO package

clear all

prescription = 'habex';

lambda0 = 550e-9;    %--Central wavelength of the whole spectral bandpass [meters]
lambda_um = lambda0 * 1e6;
optval.lambda0_um = 0.550;

%--Focal planes
res = 3;
FOV = 30;
Rout = 26; 
optval.nout = ceil_even(1 + res*(2*FOV)); %  dimensions of output in pixels (overrides output_dim0)
optval.final_sampling_lam0 = 1/res;	%   final sampling in lambda0/D
optval.use_field_stop = true;	%-- use field stop (0 = no stop)
optval.field_stop_radius = Rout;   %-- field stop radius in lam0/D

% %--Pupil Plane Resolutions
NbeamFull = 62*7;
optval.pupil_diam_pix = NbeamFull;
gridsize = 2^ceil(log2(NbeamFull)); 

optval.map_dir = '/Users/ajriggs/Documents/habex/maps/';	%-- directory containing optical surface error maps
optval.normLyotDiam = 0.95;
optval.vortexCharge = 6;

optval.xoffset = 0;
optval.yoffset = 0;

%% No Correction
%--Phase Retrieval Pupil

optval.use_pr = true;
optval.pr_pupil_diam_pix = 248;

optval.use_dm1 = false;
optval.use_errors = 1;		%-- 1 = use optical surface errors, 0 = none
optval.use_fpm = 0;		%-- use focal plane mask (0 = no FPM)
optval.use_lyot_stop = 1;	%-- use Lyot stop (0 = no stop)
optval.use_field_stop = 0;	%-- use field stop (0 = no stop)
% [Epup, sampling_m]  = habex_vortex(lambda_m, gridsize, optval);
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
% [EforNorm, sampling_m]  = habex_vortex(lambda_m, gridsize, optval);
[EforNorm, sampling_m] = prop_run(prescription, lambda_um, gridsize, 'quiet', 'passvalue', optval);

IforNorm = abs(EforNorm).^2;
I00 = max(IforNorm(:));
figure(3); imagesc(log10(IforNorm/I00),[-7 0]); axis xy equal tight; colorbar; title('PSF for Normalization)', 'Interpreter','Latex'); set(gca,'Fontsize',20); drawnow;

% Coronagraphic PSF
optval.use_errors = 1;		%-- 1 = use optical surface errors, 0 = none
optval.use_fpm = 1;		%-- use focal plane mask (0 = no FPM)
optval.use_lyot_stop = 1;	%-- use Lyot stop (0 = no stop)
optval.use_field_stop = 1;	%-- use field stop (0 = no stop)
% [Ecoro, sampling_m]  = habex_vortex(lambda_m, gridsize, optval);
[Ecoro, sampling_m] = prop_run(prescription, lambda_um, gridsize, 'quiet', 'passvalue', optval);
Inorm = abs(Ecoro).^2 / I00;
figure(4); imagesc(log10(Inorm), [-7 -2]); axis xy equal tight; colorbar;
title('Coronagraphic PSF before Flattening', 'Fontsize', 16); drawnow;


%% Corrected

optval.use_dm1 = true;
optval.dm1 = fitsread([optval.map_dir, 'flat_map.fits']);

%--Phase Retrieval Pupil
optval.use_pr = true;
optval.pr_pupil_diam_pix = 248;
optval.use_errors = 1;		%-- 1 = use optical surface errors, 0 = none
optval.use_fpm = 0;		%-- use focal plane mask (0 = no FPM)
optval.use_lyot_stop = 1;	%-- use Lyot stop (0 = no stop)
optval.use_field_stop = 0;	%-- use field stop (0 = no stop)
% [Epup, sampling_m]  = habex_vortex(lambda_um, gridsize, optval);
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
% [EforNorm, sampling_m]  = habex_vortex(lambda_m, gridsize, optval);
[EforNorm, sampling_m] = prop_run(prescription, lambda_um, gridsize, 'quiet', 'passvalue', optval);
IforNorm = abs(EforNorm).^2;
I00 = max(IforNorm(:));
figure(13); imagesc(log10(IforNorm/I00),[-7 0]); axis xy equal tight; colorbar; title('PSF for Normalization)', 'Interpreter','Latex'); set(gca,'Fontsize',20); drawnow;

% Coronagraphic PSF
optval.use_errors = 1;		%-- 1 = use optical surface errors, 0 = none
optval.use_fpm = 1;		%-- use focal plane mask (0 = no FPM)
optval.use_lyot_stop = 1;	%-- use Lyot stop (0 = no stop)
optval.use_field_stop = 1;	%-- use field stop (0 = no stop)
% [Ecoro, sampling_m]  = habex_vortex(lambda_m, gridsize, optval);
[Ecoro, sampling_m] = prop_run(prescription, lambda_um, gridsize, 'quiet', 'passvalue', optval);
Inorm = abs(Ecoro).^2 / I00;
figure(14); imagesc(log10(Inorm), [-7 -2]); axis xy equal tight; colorbar;
title('Coronagraphic PSF after Flattening', 'Fontsize', 16);
drawnow;
