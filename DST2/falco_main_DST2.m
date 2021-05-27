% Copyright 2021, by the California Institute of Technology.
% ALL RIGHTS RESERVED. United States Government Sponsorship acknowledged. 
% Any commercial use must be negotiated with the Office of Technology
% Transfer at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Script to perform WFSC for the DST2 model.
%
% Matlab's path must know the locations of 
% - The DST2 Prescription Code:   'DST2.m'
% - The DST2 Falco Defaults Code: 'falco_defaults_DST2.m'
% - Polar Transform Code:         'polarTransform.m' (optional)
% - The FALCO Package Directory:
% - The OAP and DM surface files
%       - 'OAP#_#-WFE.fits'
%
% In falco_defaults_DST2.m, change the value of mp.full.map_dir to be
% correct for your computer.
% -------------------------------------------------------------------------

clear

%% Define Necessary Paths on Your Computer System

%--Required packages are FALCO and PROPER. 
% Add FALCO to the MATLAB path with the command:  addpath(genpath(full_path_to_falco)); savepath;

% Add path to DM quilting maps
path_to_this_file = fileparts(mfilename('fullpath'));
addpath([path_to_this_file filesep 'quilting'])

%% Load default model parameters

falco_defaults_DST2

%% Overwrite default values as desired

%%--Special Computational Settings
mp.flagParfor = true; %--whether to use parfor for Jacobian calculation
mp.flagPlot = true;

%--Record Keeping
mp.SeriesNum = 1;
mp.TrialNum = 0;

%--DEBUGGING (with just 1 wavelength)
mp.fracBW = 0.01;       %--fractional bandwidth of the whole bandpass (Delta lambda / lambda0)
mp.Nsbp = 1;            %--Number of sub-bandpasses to divide the whole bandpass into for estimation and control
mp.flagParfor = false; %--whether to use parfor for Jacobian calculation


%% Generate the DM quilting maps at the current sampling and DM registrations

mp.runLabel = '_'; % temporary definition
[mp, ~] = falco_flesh_out_workspace(mp); % need to run this to fill out mp.dm1 and mp.dm2
gen_quilting_map_for_proper(mp); % no direct output--writes DM maps to files 


%% Obtain DM1 Commands to Flatten the Wavefront Prior to WFSC
% Run this section, on its own, if the OAP/DM surface data changes at
% all. Only needs to be run once to generate the correct flats_map.fits

optval = mp.full;
optval.use_pr = true;
optval.use_field_stop = false;
optval.use_lyot_stop = false;
optval.use_fpm = false;
optval.pr_pupil_diam_pix = mp.P1.full.Nbeam;

optval.use_dm1 = false;
optval.use_dm2 = false;
 
[Epup, sampling_m] = DST2(mp.lambda0, mp.full.gridsize, optval);
 
mask = zeros(size(Epup));
mask(abs(Epup) > 0.05*max(abs(Epup(:)))) = 1;
EpupPhase = mask .* angle(Epup);
surfaceToFit = -0.5*mask.*EpupPhase*(mp.lambda0/(2*pi));
 
figure(1); imagesc(abs(Epup)); axis xy equal tight; colorbar; drawnow;
figure(2); imagesc(angle(Epup)); axis xy equal tight; colorbar; drawnow;
figure(3); imagesc(surfaceToFit); axis xy equal tight; colorbar; drawnow;

% mp.runLabel = 'temp_';
% [mp, out] = falco_flesh_out_workspace(mp);
mp.dm1.inf0 = fitsread(mp.dm1.inf_fn);
mp.dm1.dm_spacing = 400e-6;
mp.dm1.dx_inf0 = mp.dm1.dm_spacing / 10;
mp.dm1.dx = mp.P2.D/mp.P1.full.Nbeam; %mp.P2.D/mp.P2.compact.Nbeam;
mp.dm1.centering = 'pixel';

flatMapMeters = falco_fit_dm_surf(mp.dm1, surfaceToFit);
figure(4); imagesc(flatMapMeters); axis xy equal tight; colorbar; drawnow;

fitswrite(flatMapMeters, [mp.full.map_dir filesep 'flat_map.fits'])


%% 
mp.full.dm1.flatmap = flatMapMeters; %fitsread([mp.full.map_dir filesep 'flat_map.fits']);
mp.full.dm2.flatmap = 0;

%% Obtain the phase retrieval phase with the flattened pupil phase.

optval = mp.full;
optval.use_pr = true;
optval.use_field_stop = false;
optval.use_lyot_stop = false;
optval.use_fpm = false;
optval.pr_pupil_diam_pix = mp.P1.compact.Nbeam;

optval.xoffset = 0; 
optval.yoffset = 0;
optval.use_dm1 = true;
optval.dm1 = flatMapMeters; %fitsread([optval.map_dir, 'flat_map.fits']); % Be sure the location to this is set properly in defaults code.
optval.dm2 = zeros(50, 50);

nout = mp.P1.full.Narr;     % nout > pupil_diam_pix
optval.output_dim = nout;   % Get the Input Pupil's E-field

if(mp.Nsbp==1)
    lambdaFacs = 1;
else
    lambdaFacs = linspace(1-mp.fracBW/2,1+mp.fracBW/2,mp.Nsbp);
end

prescription = 'DST2';

Nc = ceil_even(mp.P1.compact.Nbeam) + 20;
mp.P1.compact.E = zeros(Nc, Nc, mp.Nsbp);
for si=1:mp.Nsbp
    
    lambda_um = 1e6*mp.lambda0*lambdaFacs(si);

    [fldC, ~] = prop_run(prescription, lambda_um, mp.full.gridsize, 'quiet', 'passvalue', optval);
    fldC = pad_crop(fldC, Nc);
    
    if(mp.flagPlot)
        figure(607+si-1); imagesc(angle(fldC)); axis xy equal tight; colorbar; colormap hsv; drawnow;
        figure(608); imagesc(abs(fldC)); axis xy equal tight; colorbar; colormap parula; drawnow;
    end
    
    mp.P1.compact.E(:,:,si) = propcustom_relay(fldC, 1);
    
end


%% Generate the label associated with this trial

mp.runLabel = ['Series',num2str(mp.SeriesNum,'%04d'),'_Trial',num2str(mp.TrialNum,'%04d_'),...
    mp.coro,'_',mp.whichPupil,'_',num2str(numel(mp.dm_ind)),'DM',num2str(mp.dm1.Nact),'_z',num2str(mp.d_dm1_dm2),...
    '_IWA',num2str(mp.Fend.corr.Rin),'_OWA',num2str(mp.Fend.corr.Rout),...
    '_',num2str(mp.Nsbp),'lams',num2str(round(1e9*mp.lambda0)),'nm_BW',num2str(mp.fracBW*100),...
    '_',mp.controller];


%% Perform the Wavefront Sensing and Control

[mp, out] = falco_flesh_out_workspace(mp);

mp.P1.compact.mask = ones(size(mp.P1.compact.mask));

[mp, out] = falco_wfsc_loop(mp, out);



%% Contrast Curve Calculation
% This section takes the output of the entire simulation and
% generates a contrast curve. It takes a little while, so be patient. 

beam_ratio = mp.full.pupil_diam_pix/mp.full.gridsize;

Ndim = 1;

conversionMap = falco_compute_NI_to_contrast_map(mp, Ndim);
 
Ifinal = falco_get_summed_image(mp);
Cfinal = Ifinal ./ conversionMap;
 
figure(320); imagesc(conversionMap); axis xy equal tight; colorbar;
figure(321); imagesc(log10(Ifinal), [-10, -4]); axis xy equal tight; colorbar;

x = ((mp.full.gridsize/2):1:((mp.full.gridsize-1)/2))*beam_ratio;

figure(322); imagesc(mp.Fend.xisDL,mp.Fend.etasDL,real(log10(Cfinal)), [-10, -4]); axis xy equal tight; colorbar;
xlabel('\lambda/d'); ylabel('\lambda/d'); title('DST 2 Coronagraph Contrast')

N = 74;
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