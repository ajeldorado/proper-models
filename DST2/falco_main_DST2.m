% Copyright 2021, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Script to perform WFSC for the c vortex model.
%
% FALCO Written by A.J. Riggs (JPL, CIT) in February 2020.
    % DST 2 Model written by Matthew Noyes (JPL) January 2021
    % FALCO adapted for DST 2 Model by A.J. Riggs and Matt Noyes in
    % February 2021
%
% Matlab's path must know the locations of 
% - The DST2 Prescription Code:   'DST2.m'
% - The DST2 Falco Defaults Code: 'falco_defaults_DST2.m'
% - Polar Transform Code:         'polarTransform.m' (optional)
% - The PROPER Library Directory: (\falco-matlab-master)
% - The FALCO Package Directory:  (\lib_external\proper)
% - The OAP and DM surface files
%       - 'flat_map.fits'
%       - 'OAP#_#-WFE.fits'
%
% In falco_defaults_DST2.m, change the value of mp.full.map_dir to be
% correct for your computer.
% -------------------------------------------------------------------------

clear, close all; clc



%% Step 1: Define Necessary Paths on Your Computer System

%--Required packages are FALCO and PROPER. 
% Add FALCO to the MATLAB path with the command:  addpath(genpath(full_path_to_falco)); savepath;
% Add PROPER to the MATLAB path with the command: addpath(full_path_to_proper); savepath;

addpath(genpath('C:\Users\matnoyes\Documents\GitHub\proper-models\'));
addpath(genpath('C:\Users\matnoyes\Documents\GitHub\falco-matlab\'));
addpath(genpath('C:\Users\matnoyes\Documents\Projects\HCIT_Camilo\DST2\Maps\OAP_FITS\'))

%% Step 2: Load default model parameters

falco_defaults_DST2

%% Step 3a: Overwrite default values as desired

%%--Special Computational Settings
mp.flagParfor = true; %--whether to use parfor for Jacobian calculation
mp.flagPlot = true;

%--Record Keeping
mp.SeriesNum = 1;
mp.TrialNum = 1;

%--DEBUGGING
mp.fracBW = 0.01;       %--fractional bandwidth of the whole bandpass (Delta lambda / lambda0)
mp.Nsbp = 1;            %--Number of sub-bandpasses to divide the whole bandpass into for estimation and control
mp.flagParfor = false; %--whether to use parfor for Jacobian calculation

%% Step 3b: Obtain the phase retrieval phase.

optval = mp.full;
optval.xoffset = 0; 
optval.yoffset = 0;
optval.use_dm1 = true;
optval.dm1 = fitsread([optval.map_dir, 'flat_map.fits']); % Be sure the location to this is set properly in defaults code.
optval.dm2 = zeros(50,50);

optval.use_pr = true;
optval.use_fpm = 0;
nout = mp.P1.full.Narr;     % nout > pupil_diam_pix
optval.output_dim = nout;   % Get the Input Pupil's E-field

if(mp.Nsbp==1)
    lambdaFacs = 1;
else
    lambdaFacs = linspace(1-mp.fracBW/2,1+mp.fracBW/2,mp.Nsbp);
end

prescription = 'DST2';

%--Get the Input Pupil's E-field
Nf = mp.full.gridsize; %--N full
Nc = mp.full.gridsize; %--N compact
mp.P1.compact.E = zeros(mp.full.gridsize, mp.full.gridsize, mp.Nsbp);
for si=1:mp.Nsbp
    lambda_um = 1e6*mp.lambda0*lambdaFacs(si);

    [fldFull, ~] = prop_run(prescription, lambda_um, nout, 'quiet', 'passvalue', optval);
    if(mp.flagPlot)
        figure(605); imagesc(angle(fldFull)); axis xy equal tight; colorbar; colormap hsv; drawnow;
        figure(606); imagesc(abs(fldFull)); axis xy equal tight; colorbar; colormap parula; drawnow;
    end

    %%--Downsampling for the compact model
    dxF = 1
    dxC = mp.P1.full.Nbeam/mp.P1.compact.Nbeam

    xF = (-Nf/2:(Nf/2-1))*dxF;
    xC = (-Nc/2:(Nc/2-1))*dxC;
    [Xf, Yf] = meshgrid(xF);
    [Xc, Yc] = meshgrid(xC);
    
    if mp.P1.full.Nbeam ~= mp.P1.compact.Nbeam
        fldC = interp2(Xf, Yf, fldFull, Xc, Yc, 'cubic', 0); %--Downsample by interpolation
    else
        fldC = fldFull;
    end
    
    if(mp.flagPlot)
        figure(607+si-1); imagesc(angle(fldC)); axis xy equal tight; colorbar; colormap hsv; drawnow;
        figure(608); imagesc(abs(fldC)); axis xy equal tight; colorbar; colormap parula; drawnow;
    end
    
    mp.P1.compact.E(:,:,si) = propcustom_relay(fldC, 1);
    
end

%%

%% %% (Optional) Step 3c: Obtain DM1 Commands to Flatten the Wavefront Prior to WFSC
%  %% Run this section, on its own, if the OAP/DM surface data changes at
%  %% all. Only needs to be run once to generate the correct flats_map.fits
% 
% optval.use_pr = true;
% optval.pr_pupil_diam_pix = mp.P1.compact.Nbeam;
%  
% optval.use_dm1 = false;
% optval.use_dm2 = false;
%  
% optval.use_errors = 1;   %-- 1 = use optical surface errors, 0 = none
% % optval.use_fpm = 0;    %-- use focal plane mask (0 = no FPM)
% % optval.use_lyot_stop = 0; %-- use Lyot stop (0 = no stop)
% % optval.use_field_stop = 0;  %-- use field stop (0 = no stop)
% [Epup, sampling_m] = DST2(mp.lambda0, mp.P1.full.Narr, optval);
%  
% mask = 0*Epup;
% mask(abs(Epup) > 1e-1*max(abs(Epup(:)))) = 1;
% % mask = ones(size(Epup));
%  
% surfaceToFit = -0.5*mask.*angle(Epup)*(mp.lambda0/(2*pi));
%  
% figure(1); imagesc(abs(Epup)); axis xy equal tight; colorbar; 
% % title('', 'Fontsize', 16);
% drawnow;
% figure(2); imagesc(surfaceToFit); axis xy equal tight; colorbar; 
%  
%  
% mp.dm1.inf0 = fitsread(mp.dm1.inf_fn);
% mp.dm1.dm_spacing = 400e-6;
% mp.dm1.dx_inf0 = mp.dm1.dm_spacing/10;
% mp.dm1.dx = mp.P2.D/mp.P2.compact.Nbeam;
% mp.dm1.centering = 'pixel';
% V0 = falco_fit_dm_surf(mp.dm1,surfaceToFit);
% figure(3); imagesc(V0); axis xy equal tight; colorbar; drawnow;
%  
% fitswrite(V0, [mp.full.map_dir, 'flat_map.fits'])


%% Step 4: Generate the label associated with this trial

mp.runLabel = ['Series',num2str(mp.SeriesNum,'%04d'),'_Trial',num2str(mp.TrialNum,'%04d_'),...
    mp.coro,'_',mp.whichPupil,'_',num2str(numel(mp.dm_ind)),'DM',num2str(mp.dm1.Nact),'_z',num2str(mp.d_dm1_dm2),...
    '_IWA',num2str(mp.Fend.corr.Rin),'_OWA',num2str(mp.Fend.corr.Rout),...
    '_',num2str(mp.Nsbp),'lams',num2str(round(1e9*mp.lambda0)),'nm_BW',num2str(mp.fracBW*100),...
    '_',mp.controller];


%% Step 5: Perform the Wavefront Sensing and Control

% mp.P1.compact.mask = ones(mp.full.gridsize,mp.full.gridsize);  %Matt created this line as a means to fixx a bug. It didn't work

[mp, out] = falco_flesh_out_workspace(mp);

mp.P1.compact.mask = ones(size(mp.P1.compact.mask));

[mp, out] = falco_wfsc_loop(mp, out);

%% %% Step 6: Contrast Curve Code
%  %% This messy section takes the output of the entire simulation and
%  %% generates a contrast curve. It takes a little while, so be patient. 

gridsize = mp.full.gridsize;
beam_ratio = mp.full.beam_ratio;

Ndim = 1;

conversionMap = falco_compute_NI_to_contrast_map(mp, Ndim);
 
Ifinal = falco_get_summed_image(mp);
Cfinal = Ifinal ./ conversionMap;
 
figure(320); imagesc(conversionMap); axis xy equal tight; colorbar;
figure(321); imagesc(log10(Ifinal), [-10, -4]); axis xy equal tight; colorbar;

x = ((gridsize/2):1:((gridsize-1)/2))*beam_ratio;

figure(322); imagesc(mp.Fend.xisDL,mp.Fend.etasDL,real(log10(Cfinal)), [-10, -4]); axis xy equal tight; colorbar;
xlabel('\lambda/d'); ylabel('\lambda/d'); title('DST 2 Coronagraph Contrast')

N=74;
[pt,rvec,qvec] = polarTransform(Ifinal, [N/2+1,N/2+1], N/2, N/2+1, 360, 'bilinear');
rad_pro = mean(pt,2);

ldd = mp.Fend.xisDL;

figure(103); % plotting normalized intesity as a function of the mean distance from the center
semilogy(ldd(74/2:74),rad_pro); title('log_1_0(Normalized Intensity)');
xlabel('\rho [\lambda/D]', 'FontSize', 20); ylabel('log_1_0(Normalized Intensity)', 'FontSize', 20); set(gca,'FontSize',12); 

N=74;
[pt,rvec,qvec] = polarTransform(Cfinal, [N/2+1,N/2+1], N/2, N/2+1, 360, 'bilinear');
rad_pro = mean(pt,2);

ldd = mp.Fend.xisDL;

figure(104); % plotting normalized intesity as a function of the mean distance from the center
semilogy(ldd(74/2:74),rad_pro); title('log_1_0(Contrast)');

xlabel('\rho [\lambda/D]', 'FontSize', 20); ylabel('log_1_0(Contrast)', 'FontSize', 20); set(gca,'FontSize',12); 