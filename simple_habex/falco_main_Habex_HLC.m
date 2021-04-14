% Copyright 2020, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Script to perform WFSC for the Habex HLC model.
%
% Written by A.J. Riggs (JPL, CIT) in February 2020.
%
% Matlab's path must know the locations of 
% - the Habex prescription 'habex_multi_coro.m'
% - the FALCO package
%
% In falco_defaults_Habex_HLC.m, change the value of mp.full.map_dir to be
% for your computer.
% -------------------------------------------------------------------------

clear

%% Step 1: Define Necessary Paths on Your Computer System

%--Required packages are FALCO and PROPER. 
% Add FALCO to the MATLAB path with the command:  addpath(genpath(full_path_to_falco)); savepath;

%%--Output Data Directories (Comment these lines out to use defaults within falco-matlab/data/ directory.)
% mp.path.config = ; %--Location of config files and minimal output files. Default is [mp.path.falco filesep 'data' filesep 'brief' filesep]
% mp.path.ws = ; % (Mostly) complete workspace from end of trial. Default is [mp.path.falco filesep 'data' filesep 'ws' filesep];
% mp.flagSaveWS = false;  %--Whether to save out entire (large) workspace at the end of trial. Default is false


%% Step 2: Load default model parameters

falco_defaults_Habex_HLC


%% Step 3: Overwrite default values as desired

% mp.P1.full.Narr = 1024; % gridsize for PROPER

%%--Special Computational Settings
mp.flagParfor = true; %--whether to use parfor for Jacobian calculation
mp.flagPlot = true;

%--Record Keeping
mp.SeriesNum = 1;
mp.TrialNum = 1;

%%--[OPTIONAL] Start from a previous FALCO trial's DM settings
% fn_prev = 'ws_Series0002_Trial0001_HLC_WFIRST20180103_2DM48_z1_IWA2.7_OWA10_6lams575nm_BW12.5_EFC_30its.mat';
% temp = load(fn_prev,'out');
% mp.dm1.V = temp.out.DM1V;
% mp.dm2.V = temp.out.DM2V;
% clear temp

%--DEBUGGING
mp.fracBW = 0.01;       %--fractional bandwidth of the whole bandpass (Delta lambda / lambda0)
mp.Nsbp = 1;            %--Number of sub-bandpasses to divide the whole bandpass into for estimation and control
mp.flagParfor = false; %--whether to use parfor for Jacobian calculation

% mp.Nsbp = 3;

%%
if(mp.Nsbp==1)
    lambdaFacs = 1;
else
    lambdaFacs = linspace(1-mp.fracBW/2,1+mp.fracBW/2,mp.Nsbp);
end

lam_occ = lambdaFacs*mp.lambda0;

mp.F3.compact.Nxi = 250;%280;%1024;%280;%350; %--Crop down to minimum size of the spot
mp.F3.compact.Neta = mp.F3.compact.Nxi;
mp.compact.FPMcube = zeros(mp.F3.compact.Nxi,mp.F3.compact.Nxi,mp.Nsbp);
fpm_axis = 's';

for si=1:mp.Nsbp
    
    lambda_um = 1e6*mp.lambda0*lambdaFacs(si);
    fn_p_r = [mp.full.mask_dir '/run819_roman_occ_lam' num2str(lam_occ(si),12) 'theta5.5pol'   fpm_axis   '_' 'real.fits'];
    fn_p_i = [mp.full.mask_dir '/run819_roman_occ_lam' num2str(lam_occ(si),12) 'theta5.5pol'   fpm_axis   '_' 'imag.fits'];   
        
    dsfac = 4;
    occ0 = complex(fitsread(fn_p_r),fitsread(fn_p_i));
    N0 = length(occ0);
    occ0 = pad_crop(occ0, N0*4);
    focc = fft2(ifftshift(occ0));
    focc = fftshift(focc);
    focc = pad_crop(focc, N0);
    focc = fftshift(focc);
    occ = fftshift(ifft2(focc))/dsfac^2;
%     occ = pad_crop(occ, N0);
    figure(101); imagesc(abs(occ0)); axis xy equal tight; drawnow;
    figure(102); imagesc(abs(occ)); axis xy equal tight; drawnow;
    
%     occ = propcustom_relay(occ, 1);
    mp.compact.FPMcube(:,:,si) = pad_crop(occ, mp.F3.compact.Nxi);
    
end
mp.F3.compact.res = 1024/mp.P1.compact.Nbeam;

% mp.F3.compact.Nxi = 1024;%280;%1024;%280;%350; %--Crop down to minimum size of the spot
% mp.F3.compact.Neta = mp.F3.compact.Nxi;
% mp.compact.FPMcube = zeros(mp.F3.compact.Nxi,mp.F3.compact.Nxi,mp.Nsbp);
% fpm_axis = 's';
% 
% for si=1:mp.Nsbp
%     lambda_um = 1e6*mp.lambda0*lambdaFacs(si);
%     fn_p_r = [mp.full.mask_dir '/run819_roman_occ_lam' num2str(lam_occ(si),12) 'theta5.5pol'   fpm_axis   '_' 'real.fits'];
%     fn_p_i = [mp.full.mask_dir '/run819_roman_occ_lam' num2str(lam_occ(si),12) 'theta5.5pol'   fpm_axis   '_' 'imag.fits'];   
%     
%     mp.compact.FPMcube(:,:,si) = pad_crop(complex(fitsread(fn_p_r),fitsread(fn_p_i)), mp.F3.compact.Nxi);
% % %     mp.compact.FPMcube(:,:,si) = rot90(mp.compact.FPMcube(:,:,si), 2);
%     
% %     dsfac = 4;
% %     occ0 = complex(fitsread(fn_p_r),fitsread(fn_p_i));
% %     N0 = length(occ0);
% %     occ0 = pad_crop(occ0, N0*4);
% %     focc = fft2(ifftshift(occ0));
% %     focc = fftshift(focc);
% %     focc = pad_crop(focc, N0);
% %     focc = fftshift(focc);
% %     occ = fftshift(ifft2(focc))/dsfac^2;
% % %     occ = pad_crop(occ, N0);
% %     figure(101); imagesc(abs(occ0)); axis xy equal tight; drawnow;
% %     figure(102); imagesc(abs(occ)); axis xy equal tight; drawnow;
% %     
% %     mp.compact.FPMcube(:,:,si) = pad_crop(occ, mp.F3.compact.Nxi);
% 
%     
% end
% 
% mp.F3.compact.res = 4 * 1024/mp.P1.compact.Nbeam;

%% Visually check the FPM cropping
for si=1:mp.Nsbp
   figure(100); imagesc(abs(mp.compact.FPMcube(:,:,si))); axis xy equal tight; colorbar; colormap('parula'); drawnow; pause(0.1); 
%    figure(100); imagesc(angle(mp.compact.FPMcube(:,:,si))); axis xy equal tight; colorbar; colormap('hsv'); drawnow; pause(0.1); 
end

%%

%% Step 3b: Obtain the pupil phase just upstream of the FPM.

use_phase_retrieval = true;

if use_phase_retrieval
    mp.P1.compact.mask = ones(size(mp.P1.compact.mask)); % Don't double-apply the telescope aperture with phase retrieval and mp.P1.compact.mask

    %mp.full.input_field_rootname = '/Users/ajriggs/Repos/falco-matlab/data/maps/input_full';
    optval = mp.full;
    optval.use_hlc_dm_patterns = false;
    optval.xoffset = 0;
    optval.yoffset = 0;
    % optval.use_dm1 = true;
    % optval.dm1 = fitsread([optval.map_dir, 'flat_map.fits']);

    optval.use_pr = true;
    %optval.end_at_fpm_exit_pupil = 1;
    %optval.output_field_rootname = [fileparts(mp.full.input_field_rootname) filesep 'fld_at_xtPup'];
    optval.use_fpm = 0;
    nout = mp.P1.full.Narr;			% nout > pupil_diam_pix
    optval.output_dim = nout;%% Get the Input Pupil's E-field
    optval.pr_pupil_diam_pix = mp.P1.compact.Nbeam;

    if(mp.Nsbp==1)
        lambdaFacs = 1;
    else
        lambdaFacs = linspace(1-mp.fracBW/2,1+mp.fracBW/2,mp.Nsbp);
    end

    % prescription = 'habex';

    %--Get the Input Pupil's E-field
    Nf = nout; %--N full
    Nc = ceil_even( (mp.P1.compact.Nbeam/mp.P1.full.Nbeam)*Nf ); %--N compact
    mp.P1.compact.E = ones(Nc, Nc, mp.Nsbp); %--Initialize
    %mp.P1.compact.E = ones(mp.P1.compact.Nbeam+2, mp.P1.compact.Nbeam+2, mp.Nsbp); %--Initialize
    for si=1:mp.Nsbp
        lambda_um = 1e6*mp.lambda0*lambdaFacs(si);

        [fldFull, ~] = prop_run(mp.full.prescription, lambda_um, nout, 'quiet', 'passvalue', optval);
        if(mp.flagPlot)
            figure(605); imagesc(angle(fldFull)); axis xy equal tight; colorbar; colormap hsv; drawnow;
            figure(606); imagesc(abs(fldFull)); axis xy equal tight; colorbar; colormap parula; drawnow;
        end

        %lams = num2str(lambda_um, '%6.4f');
        %pols = ['polaxis'  num2str(optval.polaxis,2)];
        %fitswrite(real(fldFull), [mp.full.input_field_rootname '_' lams 'um_' pols '_real.fits' ]);
        %fitswrite(imag(fldFull), [mp.full.input_field_rootname '_' lams 'um_' pols '_imag.fits' ]);

        %%--Downsampling for the compact model
        dxF = 1;
        dxC = mp.P1.full.Nbeam/mp.P1.compact.Nbeam;

        if dxC == 1
            fldC = fldFull;
            fldC = pad_crop(fldC, Nc);
        else
            xF = (-Nf/2:(Nf/2-1))*dxF;
            xC = (-Nc/2:(Nc/2-1))*dxC;
            [Xf, Yf] = meshgrid(xF);
            [Xc, Yc] = meshgrid(xC);
            fldC = interp2(Xf, Yf, fldFull, Xc, Yc, 'cubic', 0); %--Downsample by interpolation
            %fldC = pad_crop(fldC, ceil_even(mp.P1.compact.Nbeam+1));
            if(mp.flagPlot)
                figure(607+si-1); imagesc(angle(fldC)); axis xy equal tight; colorbar; colormap hsv; drawnow;
                figure(608); imagesc(abs(fldC)); axis xy equal tight; colorbar; colormap parula; drawnow;
            end
        end

        %--Assign to initial E-field in compact model.
    % %     temp = 0*fldC;
    % %     temp(2:end,2:end) = rot90(fldC(2:end,2:end),2);
    % %     mp.P1.compact.E(:,:,si) = temp;

        mp.P1.compact.E(:,:,si) = propcustom_relay(fldC, 1);

    end

end

%%

% E =  mp.P1.compact.E(:,:,1);
% 
% pupil = max(E(:))*propcustom_relay(mp.P1.compact.mask, 0);
% 
% whos E pupil
% 
% figure(621); imagesc(abs(E)); axis xy equal tight; colorbar; colormap parula; drawnow;
% figure(622); imagesc(abs(pupil)); axis xy equal tight; colorbar; colormap parula; drawnow;
% figure(623); imagesc(abs(pupil)-abs(E)); axis xy equal tight; colorbar; colormap parula; drawnow;


% mp.dm1.V = 1e9*fitsread([mp.full.mask_dir, 'run819_roman_dm1acts.fits']);
% mp.dm2.V = 1e9*fitsread([mp.full.mask_dir, 'run819_roman_dm2acts.fits']);



%% Obtain DM1 Commands to Flatten the Wavefront Prior to WFSC

optval.use_pr = true;
optval.pr_pupil_diam_pix = mp.P1.compact.Nbeam;

optval.use_errors = 1;		%-- 1 = use optical surface errors, 0 = none
% optval.use_fpm = 0;		%-- use focal plane mask (0 = no FPM)
% optval.use_lyot_stop = 0;	%-- use Lyot stop (0 = no stop)
% optval.use_field_stop = 0;	%-- use field stop (0 = no stop)
[Epup, sampling_m]  = habex_multi_coro(mp.lambda0, mp.P1.full.Narr, optval);

mask = 0*Epup;
mask(abs(Epup) > 1e-1*max(abs(Epup(:)))) = 1;
% mask = ones(size(Epup));

surfaceToFit = -0.5*mask.*angle(Epup)*(mp.lambda0/(2*pi));

figure(1); imagesc(abs(Epup)); axis xy equal tight; colorbar; 
% title('', 'Fontsize', 16);
drawnow;
figure(2); imagesc(surfaceToFit); axis xy equal tight; colorbar; 


mp.dm1.inf0 = fitsread(mp.dm1.inf_fn);
mp.dm1.dm_spacing = 400e-6;
mp.dm1.dx_inf0 = mp.dm1.dm_spacing/10;
mp.dm1.dx = mp.P2.D/mp.P2.compact.Nbeam;
mp.dm1.centering = 'pixel';
V0 = falco_fit_dm_surf(mp.dm1,surfaceToFit);

figure(3); imagesc(V0); axis xy equal tight; colorbar; 

% fitswrite(V0, [mp.full.map_dir, 'flat_map_hlc.fits'])


%% Step 4: Generate the label associated with this trial

mp.runLabel = ['Series',num2str(mp.SeriesNum,'%04d'),'_Trial',num2str(mp.TrialNum,'%04d_'),...
    mp.coro,'_',mp.whichPupil,'_',num2str(numel(mp.dm_ind)),'DM',num2str(mp.dm1.Nact),'_z',num2str(mp.d_dm1_dm2),...
    '_IWA',num2str(mp.Fend.corr.Rin),'_OWA',num2str(mp.Fend.corr.Rout),...
    '_',num2str(mp.Nsbp),'lams',num2str(round(1e9*mp.lambda0)),'nm_BW',num2str(mp.fracBW*100),...
    '_',mp.controller];


%% Step 5: Perform the Wavefront Sensing and Control

[mp, out] = falco_flesh_out_workspace(mp);

[mp, out] = falco_wfsc_loop(mp, out);
