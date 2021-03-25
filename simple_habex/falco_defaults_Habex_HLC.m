% Copyright 2018-2020, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
% %--Initialize some structures if they don't already exist

%% Misc

%--Record Keeping
mp.SeriesNum = 867;
mp.TrialNum = 5309;

%--Special Computational Settings
mp.flagParfor = false;
mp.useGPU = false;
mp.flagPlot = false;

%--General
mp.centering = 'pixel';

%--Method of computing core throughput:
% - 'HMI' for energy within half-max isophote divided by energy at telescope pupil
% - 'EE' for encircled energy within a radius (mp.thput_radius) divided by energy at telescope pupil
mp.thput_metric = 'EE'; 
mp.thput_radius = 0.7; %--photometric aperture radius [lambda_c/D]. Used ONLY for 'EE' method.
mp.thput_eval_x = 10; % x location [lambda_c/D] in dark hole at which to evaluate throughput
mp.thput_eval_y = 0; % y location [lambda_c/D] in dark hole at which to evaluate throughput

%--Where to shift the source to compute the intensity normalization value.
mp.source_x_offset_norm = 10;  % x location [lambda_c/D] in dark hole at which to compute intensity normalization
mp.source_y_offset_norm = 0;  % y location [lambda_c/D] in dark hole at which to compute intensity normalization

%% Bandwidth and Wavelength Specs

mp.lambda0 = 575e-9;    %--Central wavelength of the whole spectral bandpass [meters]
mp.fracBW = 0.10;       %--fractional bandwidth of the whole bandpass (Delta lambda / lambda0)
mp.Nsbp = 9;            %--Number of sub-bandpasses to divide the whole bandpass into for estimation and control
mp.Nwpsbp = 1;          %--Number of wavelengths to used to approximate an image in each sub-bandpass

%% Wavefront Estimation

%--Estimator Options:
% - 'perfect' for exact numerical answer from full model
% - 'pwp-bp-square' for pairwise probing with batch process estimation in a
% square region for one star [original functionality of 'pwp-bp' prior to January 2021]
% - 'pwp-bp' for pairwise probing in the specified rectangular regions for
%    one or more stars
% - 'pwp-kf' for pairwise probing with Kalman filter [NOT TESTED YET]
mp.estimator = 'perfect';

%--For pairwise probing estimation with mp.estimator='pwp-bp-square'
mp.est.probe.Npairs = 3;     % Number of pair-wise probe PAIRS to use.
mp.est.probe.whichDM = 1;    % Which DM # to use for probing. 1 or 2. Default is 1
mp.est.probe.radius = 12;    % Max x/y extent of probed region [lambda/D].
mp.est.probe.xOffset = 0;   % offset of probe center in x [actuators]. Use to avoid central obscurations.
mp.est.probe.yOffset = 10;    % offset of probe center in y [actuators]. Use to avoid central obscurations.
mp.est.probe.axis = 'alternate';     % which axis to have the phase discontinuity along [x or y or xy/alt/alternate]
mp.est.probe.gainFudge = 1;     % empirical fudge factor to make average probe amplitude match desired value.

%% Wavefront Control: General

mp.ctrl.flagUseModel = true; %--Perform a model-based grid search using the compact model

%--Threshold for culling weak actuators from the Jacobian:
mp.logGmin = -6;  % 10^(mp.logGmin) used on the intensity of DM1 and DM2 Jacobians to weed out the weakest actuators

%--Zernikes to suppress with controller
mp.jac.zerns = 1;  %--Which Zernike modes to include in Jacobian. Given as the max Noll index. Always include the value "1" for the on-axis piston mode.
mp.jac.Zcoef = 1e-9*ones(size(mp.jac.zerns)); %--meters RMS of Zernike aberrations. (piston value is reset to 1 later)
    
%--Zernikes to compute sensitivities for
mp.eval.indsZnoll = []; %--Noll indices of Zernikes to compute values for
%--Annuli to compute 1nm RMS Zernike sensitivities over. Columns are [inner radius, outer radius]. One row per annulus.
mp.eval.Rsens = [2,3; 3,4; 4,5]; 

%--Grid- or Line-Search Settings
mp.ctrl.log10regVec = -6:1/2:-2; %--log10 of the regularization exponents (often called Beta values)
mp.ctrl.dmfacVec = 1;            %--Proportional gain term applied to the total DM delta command. Usually in range [0.5,1].
   
%--Spatial pixel weighting
mp.WspatialDef = [];% [3, 4.5, 3]; %--spatial control Jacobian weighting by annulus: [Inner radius, outer radius, intensity weight; (as many rows as desired)]

%--DM weighting
mp.dm1.weight = 1;
mp.dm2.weight = 1;

%--Voltage range restrictions
mp.dm1.maxAbsV = 1000;  %--Max absolute voltage (+/-) for each actuator [volts] %--NOT ENFORCED YET
mp.dm2.maxAbsV = 1000;  %--Max absolute voltage (+/-) for each actuator [volts] %--NOT ENFORCED YET
mp.maxAbsdV = 1000;     %--Max +/- delta voltage step for each actuator for DMs 1 and 2 [volts] %--NOT ENFORCED YET


%% Wavefront Control: Controller Specific
% Controller options: 
%  - 'gridsearchEFC' for EFC as an empirical grid search over tuning parameters
%  - 'plannedEFC' for EFC with an automated regularization schedule
mp.controller = 'gridsearchEFC';

% % % GRID SEARCH EFC DEFAULTS     
%--WFSC Iterations and Control Matrix Relinearization
mp.Nitr = 5; %--Number of estimation+control iterations to perform
mp.relinItrVec = 1:mp.Nitr;  %--Which correction iterations at which to re-compute the control Jacobian
mp.dm_ind = [1 2]; %--Which DMs to use


%% Deformable Mirrors: Influence Functions
%--Influence Function Options:
% - 'influence_dm5v2.fits' for one type of Xinetics DM
% - 'influence_BMC_2kDM_400micron_res10.fits' for BMC 2k DM
% - 'influence_BMC_kiloDM_300micron_res10_spline.fits' for BMC kiloDM
mp.dm1.inf_fn = 'influence_dm5v2.fits';
mp.dm2.inf_fn = 'influence_dm5v2.fits';

mp.dm1.dm_spacing = 400e-6; %--User defined actuator pitch [meters]
mp.dm2.dm_spacing = 400e-6; %--User defined actuator pitch [meters]
pitchRatio = 0.4/0.9906; % BMC 2k DM pitch over Xinetics DM pitch

mp.dm1.inf_sign = '+';
mp.dm2.inf_sign = '+';

%% Deformable Mirrors: Optical Layout Parameters

%--DM1 parameters
mp.dm1.Nact = 64;               % # of actuators across DM array
mp.dm1.VtoH = 1e-9*ones(mp.dm1.Nact);  % gains of all actuators [nm/V of free stroke]
mp.dm1.xtilt = 0;               % for foreshortening. angle of rotation about x-axis [degrees]
mp.dm1.ytilt = 9.65;               % for foreshortening. angle of rotation about y-axis [degrees]
mp.dm1.zrot = 0;                % clocking of DM surface [degrees]
mp.dm1.xc = (mp.dm1.Nact/2 - 1/2);       % x-center location of DM surface [actuator widths]
mp.dm1.yc = (mp.dm1.Nact/2 - 1/2);       % y-center location of DM surface [actuator widths]
mp.dm1.edgeBuffer = 1;          % max radius (in actuator spacings) outside of beam on DM surface to compute influence functions for. [actuator widths]

%--DM2 parameters
mp.dm2.Nact = 64;               % # of actuators across DM array
mp.dm2.VtoH = 1e-9*ones(mp.dm1.Nact);  % gains of all actuators [nm/V of free stroke]
mp.dm2.xtilt = 0;               % for foreshortening. angle of rotation about x-axis [degrees]
mp.dm2.ytilt = 9.65;               % for foreshortening. angle of rotation about y-axis [degrees]
mp.dm2.zrot = 0;                % clocking of DM surface [degrees]
mp.dm2.xc = (mp.dm2.Nact/2 - 1/2);       % x-center location of DM surface [actuator widths]
mp.dm2.yc = (mp.dm2.Nact/2 - 1/2);       % y-center location of DM surface [actuator widths]
mp.dm2.edgeBuffer = 1;          % max radius (in actuator spacings) outside of beam on DM surface to compute influence functions for. [actuator widths]

%--Aperture stops at DMs
mp.flagDM1stop = false; %--Whether to apply an iris or not
mp.dm1.Dstop = mp.dm1.Nact*mp.dm1.dm_spacing;  %--Diameter of iris [meters]
mp.flagDM2stop = false;  %--Whether to apply an iris or not
mp.dm2.Dstop = mp.dm2.Nact*mp.dm1.dm_spacing;   %--Diameter of iris [meters]

%--DM separations
mp.d_P2_dm1 = 0;        % distance (along +z axis) from P2 pupil to DM1 [meters]
% mp.d_dm1_dm2 = 0.32;   % distance between DM1 and DM2 [meters]
mp.d_dm1_dm2 = 1.0 * (pitchRatio*pitchRatio);  % distance between DM1 and DM2 [meters]

%% Optical Layout: All models

%--Key Optical Layout Choices
mp.flagSim = true;      %--Simulation or not
mp.layout = 'proper';  %--Which optical layout to use
mp.coro = 'hlc';
mp.flagApod = false;    %--Whether to use an apodizer or not
mp.flagDMwfe = false;  %--Whether to use BMC DM quilting maps

%--Final Focal Plane Properties
mp.Fend.res = 3; %--Sampling [ pixels per lambda0/D]
mp.Fend.FOV = 14; %--half-width of the field of view in both dimensions [lambda0/D]

%--Correction and scoring region definition
mp.Fend.corr.Rin  = 2.8;   % inner radius of dark hole correction region [lambda0/D]
mp.Fend.corr.Rout = 12;  % outer radius of dark hole correction region [lambda0/D]
mp.Fend.corr.ang  = 180;  % angular opening of dark hole correction region [degrees]

mp.Fend.score.Rin  = mp.Fend.corr.Rin;  % inner radius of dark hole scoring region [lambda0/D]
mp.Fend.score.Rout = mp.Fend.corr.Rout;  % outer radius of dark hole scoring region [lambda0/D]
mp.Fend.score.ang  = mp.Fend.corr.ang;  % angular opening of dark hole scoring region [degrees]

mp.Fend.sides = 'lr'; %--Which side(s) for correction: 'left', 'right', 'top', 'up', 'bottom', 'down', 'lr', 'rl', 'leftright', 'rightleft', 'tb', 'bt', 'ud', 'du', 'topbottom', 'bottomtop', 'updown', 'downup'

%% Optical Layout: Compact Model (and Jacobian Model)
% NOTE for HLC and LC: Lyot plane resolution must be the same as input pupil's in order to use Babinet's principle

%--Focal Lengths
mp.fl = 1; %--[meters] Focal length value used for all FTs in the compact model. Don't need different values since this is a Fourier model.

%--Pupil Plane Diameters
mp.P2.D = 1.00215*64e-3 * pitchRatio; % meters
mp.P3.D = mp.P2.D;
mp.P4.D = mp.P2.D;

%--Pupil Plane Resolutions
mp.P1.compact.Nbeam = 314.581; 
mp.P2.compact.Nbeam = mp.P1.compact.Nbeam;
mp.P3.compact.Nbeam = mp.P1.compact.Nbeam;
mp.P4.compact.Nbeam = mp.P1.compact.Nbeam;  % P4 must be the same as P1 for Vortex. 

%--Number of re-imaging relays between pupil planesin compact model. Needed
%to keep track of 180-degree rotations and (1/1j)^2 factors compared to the
%full model, which probably has extra collimated beams compared to the
%compact model.
mp.Nrelay1to2 = 1;
mp.Nrelay2to3 = 1;
mp.Nrelay3to4 = 1;
mp.NrelayFend = 0;%1; %0; %--How many times to rotate the final image by 180 degrees

%% Optical Layout: Full Model 

mp.full.flagPROPER = true; %--Whether the full model is a PROPER prescription
mp.full.prescription = 'habex_multi_coro';
mp.full.gridsize = 1024; % # of points across in PROPER model 

% Values passed into the Habex PROPER model
mp.full.cor_type = 'hlc'; mp.full.mask_dir = '/Users/ajriggs/Documents/habex/run819/';	%-- directory containing Habex HLC files

mp.full.map_dir = '/Users/ajriggs/Documents/habex/maps/';	%-- directory containing optical surface error maps
mp.full.pupil_fn = 'run819_roman_pupil_rotated.fits';
mp.full.lyot_stop_fn = 'run819_roman_lyot_rotated.fits';
% mp.full.pupil_fn = 'run819_luvoir_pupil_aj.fits';
% mp.full.lyot_stop_fn = 'run819_luvoir_lyot_aj.fits';
mp.full.nlams = mp.Nsbp;
mp.full.bw = mp.fracBW;
mp.full.fpm_aoi = '5.5'; % string, AOI in degrees.
mp.full.fpm_pol = 's';

%--Focal planes
mp.full.nout = ceil_even(1 + mp.Fend.res*(2*mp.Fend.FOV)); %  dimensions of output in pixels (overrides output_dim0)
mp.full.final_sampling_lam0 = 1/mp.Fend.res;	%   final sampling in lambda0/D
mp.full.use_field_stop = false;%true;	%-- use field stop (0 = no stop)
mp.full.field_stop_radius = mp.Fend.corr.Rout;   %-- field stop radius in lam0/D

% %--Pupil Plane Resolutions
mp.P1.full.Nbeam = 314.581; %62*7;
mp.P1.full.Narr = 2^nextpow2(mp.P1.full.Nbeam);
mp.full.lambda0_um = mp.lambda0*1e6;	%-- default reference wavelength (center of bandpass) for star offsets & field stop size
mp.full.pupil_diam_pix = mp.P1.full.Nbeam;
% mp.full.grid_size = mp.P1.full.Narr;

% mp.full.use_errors = true;
% mp.full.dm1.flatmap = fitsread([mp.full.map_dir, 'flat_map.fits']);
% mp.full.dm2.flatmap = 0;

mp.full.use_errors = false;
mp.full.dm1.flatmap = 0;%fitsread([mp.full.mask_dir, 'run819_roman_dm1acts.fits']);
mp.full.dm2.flatmap = 0;%fitsread([mp.full.mask_dir, 'run819_roman_dm2acts.fits']);

mp.full.use_hlc_dm_patterns = true;
mp.full.dm1wfe_fn = 'run819_roman_dm1wfe.fits';
mp.full.dm2wfe_fn = 'run819_roman_dm2wfe.fits';


mp.dm1.wfe = fitsread([mp.full.mask_dir, mp.full.dm1wfe_fn]);
mp.dm2.wfe = fitsread([mp.full.mask_dir, mp.full.dm2wfe_fn]);
% mp.dm1.wfe = propcustom_relay(mp.dm1.wfe, 1);
% mp.dm2.wfe = propcustom_relay(mp.dm2.wfe, 1);


%% Mask Definitions

%--Pupil definition
mp.whichPupil = 'LUVOIRAFINAL';
mp.P1.IDnorm = 0.10; %--ID of the central obscuration [diameter]. Used only for computing the RMS DM surface from the ID to the OD of the pupil. OD is assumed to be 1.
% mp.P1.ODnorm = 1.00;% Outer diameter of the telescope [diameter]
mp.P1.D = 15.0; %--telescope diameter [meters]. Used only for converting milliarcseconds to lambda0/D or vice-versa.

%--Lyot stop padding
mp.P4.IDnorm = 0.2327; %--Lyot stop ID [Dtelescope]
mp.P4.ODnorm = 0.7446312; %--Lyot stop OD [Dtelescope]
mp.P4.wStrut = 0.00933782;

mp.compact.flagGenPupil = false;
mp.compact.flagGenLS = false;
mp.P1.compact.mask = fitsread([mp.full.mask_dir, 'run819_roman_pupil.fits']);
mp.P4.compact.mask = fitsread([mp.full.mask_dir, 'run819_roman_lyot.fits']);
% mp.P1.compact.mask = fitsread([mp.full.mask_dir, 'run819_roman_pupil_rotated.fits']);
% mp.P4.compact.mask = fitsread([mp.full.mask_dir, 'run819_roman_lyot_rotated.fits']);

%--Whether to generate or load various masks: compact model
% mp.dm1.wfe = fitsread([mp.full.data_dir 'hlc_20190210/' 'run461_dm1wfe.fits']);
% mp.dm2.wfe = fitsread([mp.full.data_dir 'hlc_20190210/' 'run461_dm2wfe.fits']);
