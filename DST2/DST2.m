% Copyright 2021, by the California Institute of Technology.
% ALL RIGHTS RESERVED. United States Government Sponsorship acknowledged. 
% Any commercial use must be negotiated with the Office of Technology
% Transfer at the California Institute of Technology.
% -------------------------------------------------------------------------
% DST 2 PROPER MODEL
% By Matthew Noyes - Optical Engineer | JPL
% December 13, 2020

function [wavefront, sampling_m] = DST2(lambda_m, gridsize, optval)

% Path to the mirror surface maps
map_dir = './dst2_fits_rev_b/';

% Model Switches
xoffset = false; % Introduces wavefront tilt that shifts focus in units of lam0/d
yoffset = false; % Introduces wavefront tilt that shifts focus in units of lam0/d
use_errors = true; % Switches the use of error maps on optical surfaces
use_dm_errors = false;
use_fpm = true; % Switches the use of the Focal Plane Mask
use_pinhole = false; % Switches the use of a pinhole at the FPM
pinhole_diam = 0.5; % pinhole diameter in units of lam0/d
use_lyot_stop  = true; % Switches the use of the LYOT Stop
use_field_stop = true; % Switches the use of the Field Stop
use_dm1 = false;
use_dm2 = false;
use_pr = false; % End propagation at the pupil before the FPM
which_pupil = 'exit'; % 'pre-FPM' or 'exit'

% System Variables
dm_xc = 24.5;              % Wavefront centered at corner of DM actuator (0,0 is center of 1st actuator)
dm_yc = dm_xc;
dm_sampling = 400e-6;      % DM actuator spacing (BMC 50)
diam = 0.01544*2;          % Diameter of the beam at the first reflective pupil. Found from Zemax Model.        
lambda0_m = 550e-9; % Default reference wavelength (center of bandpass) for star offsets & field stop size
pupil_diam_pix = 512; % define sampling of pupil based on having N_pix_act pixels across each DM actuator
pr_pupil_diam_pix = 512;
field_stop_radius = 30; % Field stop radius in lam0/D
final_sampling_lam0 = 0.2;   % sampling at final image plane in lam0/D per pixel
normLyotDiam = 0.95;       % Lyot stop outer diameter normalized to the beam diameter 
vortexCharge = 6;          % charge of the vortex focal plane mask
nout = 300;

% Override default optval values if passed at the beginning
if exist('optval', 'var') == 1
    if isfield(optval, 'lambda0_m');             lambda0_m = optval.lambda0_m;  end
    if isfield(optval, 'xoffset');               xoffset = optval.xoffset;  end
    if isfield(optval, 'yoffset');               yoffset = optval.yoffset;  end
    if isfield(optval, 'use_dm1');               use_dm1 = optval.use_dm1;  end
    if isfield(optval, 'dm1');                   dm1 = optval.dm1;  end
    if isfield(optval, 'use_dm2');               use_dm2 = optval.use_dm2;  end
    if isfield(optval, 'dm2');                   dm2 = optval.dm2;  end
    if isfield(optval, 'use_errors');            use_errors = optval.use_errors;  end
    if isfield(optval, 'use_dm_errors');         use_dm_errors = optval.use_dm_errors;  end
    if isfield(optval, 'use_fpm');               use_fpm = optval.use_fpm;  end
    if isfield(optval, 'use_pinhole');           use_pinhole = optval.use_pinhole;  end
    if isfield(optval, 'pinhole_diam');          pinhole_diam = optval.pinhole_diam;  end
    if isfield(optval, 'use_lyot_stop');         use_lyot_stop = optval.use_lyot_stop;  end
    if isfield(optval, 'use_field_stop');        use_field_stop = optval.use_field_stop;  end
    if isfield(optval, 'field_stop_radius');     field_stop_radius = optval.field_stop_radius;  end
    if isfield(optval, 'final_sampling_lam0');   final_sampling_lam0 = optval.final_sampling_lam0;  end
    if isfield(optval, 'nout');                  nout = optval.nout;  end
    if isfield(optval, 'normLyotDiam');          normLyotDiam = optval.normLyotDiam;  end
    if isfield(optval, 'vortexCharge');          vortexCharge = optval.vortexCharge;  end
    if isfield(optval, 'map_dir');               map_dir = optval.map_dir;  end
    if isfield(optval, 'pupil_diam_pix');        pupil_diam_pix = optval.pupil_diam_pix;  end
    if isfield(optval, 'pr_pupil_diam_pix');     pr_pupil_diam_pix = optval.pr_pupil_diam_pix;  end
    if isfield(optval, 'use_pr');                use_pr = optval.use_pr;  end
    if isfield(optval, 'which_pupil');           which_pupil = optval.which_pupil;  end    
end

pupil_ratio = pupil_diam_pix / double(gridsize);


%% Optical prescription. Units of meters.

f_OAP1 = 01.540536;                       
d_OAP1_ReflPupil = 0.733464;       % 1st Pupil Plane
d_ReflPupil_OAP2 = 0.638094;           
f_OAP2 = 0.774577;
f_OAP3 = 0.464344;
d_OAP3_DM1 = pupilDist(f_OAP2, f_OAP3, d_ReflPupil_OAP2);         % 2nd Pupil Plane
d_DM1_DM2 = 0.3;
d_DM2_OAP4 = 0.538109;
f_OAP4 = 0.641891;             % 2nd Focal Plane
f_OAP5 = 0.774577;
d_OAP5_LYOT = pupilDist(f_OAP4, f_OAP5, d_DM1_DM2+d_DM2_OAP4);    % 3rd Pupil Plane
d_LYOT_OAP6 = 0.790014;
f_OAP6 = 0.774577;             % 3rd Pupil Plane
f_OAP7 = 0.156886;
d_OAP7_pupilBeforeDetector = pupilDist(f_OAP6, f_OAP7, d_LYOT_OAP6);      % 4th & last Pupil Plane
d_pupilBeforeDetector_OAP8 = 0.477237;             
f_OAP8 = 0.475487;             % Focal length of OAP 6


%% Propagation

wavefront = prop_begin(diam, lambda_m, gridsize, pupil_ratio); 
wavefront = prop_circular_aperture(wavefront, diam/2);
wavefront = prop_define_entrance(wavefront);

% star X,Y offset in lam0/D
if (xoffset ~= 0) || (yoffset ~= 0)
	xoffset_lam = xoffset * lambda0_m / lambda_m;
	yoffset_lam = yoffset * lambda0_m / lambda_m;
	u = ((0:(gridsize-1))-gridsize/2) / (pupil_diam_pix/2);  % IDL version: u = (dindgen(gridsize)-gridsize/2) / (pupil_diam_pix/2);
	xtilt = exp( 1i * pi * u * xoffset_lam );
	ytilt = exp( 1i * pi * u * yoffset_lam );
	wavefront = prop_multiply(wavefront, ytilt.' * xtilt); % IDL version: prop_multiply, wavefront, xtilt # ytilt
	clear xtilt ytilt
end

% Back propagate from first pupil to OAP1
wavefront = prop_propagate(wavefront, -d_OAP1_ReflPupil);
if use_errors; wavefront = prop_errormap(wavefront, [map_dir filesep 'OAP1_2-WFE.fits'], 'WAVEFRONT'); end

% Begin forward propagations
wavefront = prop_propagate(wavefront, d_OAP1_ReflPupil+d_ReflPupil_OAP2);
if use_errors; wavefront = prop_errormap(wavefront, [map_dir filesep 'OAP2_2-WFE.fits'], 'WAVEFRONT'); end
wavefront = prop_lens(wavefront, f_OAP2);           

wavefront = prop_propagate(wavefront, f_OAP2+f_OAP3);
if use_errors; wavefront = prop_errormap(wavefront, [map_dir filesep 'OAP3_1-WFE.fits'], 'WAVEFRONT'); end 
wavefront = prop_lens(wavefront, f_OAP3);

wavefront = prop_propagate(wavefront, d_OAP3_DM1); 
if use_dm1; wavefront = propcustom_dm(wavefront, dm1, dm_xc, dm_yc, dm_sampling, 'inf_file', 'influence_BMC_2kDM_400micron_res10.fits'); end
if use_dm_errors; wavefront = prop_errormap(wavefront, [map_dir filesep 'dm1_surface_quilting.fits']); end

wavefront = prop_propagate(wavefront, d_DM1_DM2);
if use_dm2;  wavefront = propcustom_dm(wavefront, dm2, dm_xc, dm_yc, dm_sampling, 'inf_file', 'influence_BMC_2kDM_400micron_res10.fits'); end
if use_dm_errors; wavefront = prop_errormap(wavefront, [map_dir filesep 'dm2_surface_quilting.fits']); end

wavefront = prop_propagate(wavefront, d_DM2_OAP4);  
if use_errors; wavefront = prop_errormap(wavefront, [map_dir filesep 'OAP4_1-WFE.fits'], 'WAVEFRONT'); end 

% Get E-field at the pre-FPM pupil and then stop
if use_pr && strcmpi(which_pupil, 'pre-fpm')
         
    wavefront = prop_propagate(wavefront, -(d_DM1_DM2+d_DM2_OAP4)); % back propagate to pupil
    Epup = ifftshift(wavefront.wf);

    mag = pr_pupil_diam_pix / pupil_diam_pix;
    wavefront = prop_magnify(Epup, mag, 'SIZE_OUT', 2*ceil(0.5*mag*gridsize), 'AMP_CONSERVE');
    sampling_m = 0; % dummy value needed as output
    
    return
end


wavefront = prop_lens(wavefront, f_OAP4);
dx = wavefront.dx;
wavefront = prop_propagate(wavefront, f_OAP4);
if use_pinhole
    % Generate pinhole
    pinhole_diam = 0.5;
    pixPerLamD = 20 * (lambda0_m / lambda_m);
    inputs.pixresFPM = pixPerLamD;
    inputs.rhoInner = pinhole_diam/2;
    inputs.rhoOuter = Inf;
    pinhole = 1 - falco_gen_annular_FPM(inputs);
    fl = 1;
    dxi = 1/pixPerLamD;
    deta = 1/pixPerLamD;
    Nxi = size(pinhole, 2);
    Neta = size(pinhole, 1);

    EpupPre = ifftshift(ifft2(wavefront.wf)); % wavefront.wf is already fftshifted
    Efoc = propcustom_mft_PtoF(EpupPre, fl, lambda_m, dx, dxi, Nxi, deta, Neta);
    EpupPre = propcustom_mft_FtoP(pinhole.*Efoc, fl, lambda_m, dxi, deta, dx, gridsize);
    wavefront.wf = fft2(fftshift(EpupPre))/gridsize^2;

elseif use_fpm
    apRad = pupil_diam_pix/2.;
    useGPU = false;
    inVal = 1;      % found empirically in lamd/D
    outVal = 5;     % found empirically in lamd/D

    % 1) IFFT to previous pupil
    % 2) Use propcustom_mft_Pup2Vortex2Pup() to go to Lyot plane
    % 3) IFFT to FPM's focal plane  
    EpupPre = ifftshift(ifft2(wavefront.wf))*gridsize; % wavefront.wf is already fftshifted
    EpupPost = propcustom_mft_Pup2Vortex2Pup(EpupPre, vortexCharge, apRad, inVal, outVal, useGPU);
    wavefront.wf = ifft2(fftshift(EpupPost))*gridsize;
end

wavefront = prop_propagate(wavefront, f_OAP5); 
if use_errors; wavefront = prop_errormap(wavefront, [map_dir filesep 'OAP5_1-WFE.fits'], 'WAVEFRONT'); end 
wavefront = prop_lens(wavefront, f_OAP5); 

wavefront = prop_propagate(wavefront, d_OAP5_LYOT);
if use_lyot_stop; wavefront = prop_circular_aperture(wavefront, normLyotDiam, 'NORM'); end 


wavefront = prop_propagate(wavefront, d_LYOT_OAP6);
if use_errors; wavefront = prop_errormap(wavefront, [map_dir filesep 'OAP6_2-WFE.fits'], 'WAVEFRONT'); end 
wavefront = prop_lens(wavefront, f_OAP6);

wavefront = prop_propagate(wavefront, f_OAP6);
if use_field_stop
    r_stop = field_stop_radius * lambda0_m / lambda_m;
    wavefront = prop_circular_aperture(wavefront, r_stop/pupil_ratio*prop_get_sampling(wavefront));
end   

wavefront = prop_propagate(wavefront, f_OAP7); 
if use_errors; wavefront = prop_errormap(wavefront, [map_dir filesep 'OAP7_3-WFE.fits'], 'WAVEFRONT'); end 
wavefront = prop_lens(wavefront, f_OAP7);

wavefront = prop_propagate(wavefront, d_OAP7_pupilBeforeDetector+d_pupilBeforeDetector_OAP8);
if use_errors; wavefront = prop_errormap(wavefront, [map_dir filesep 'OAP8_2-WFE.fits'], 'WAVEFRONT'); end 

if use_pr % Get E-field at the final pupil (including aberrations from all optics)
         
    wavefront = prop_propagate(wavefront, -d_pupilBeforeDetector_OAP8); % back propagate to pupil
    Epup = ifftshift(wavefront.wf);

    mag = pr_pupil_diam_pix / pupil_diam_pix;
    wavefront = prop_magnify(Epup, mag, 'SIZE_OUT', 2*ceil(0.5*mag*gridsize), 'AMP_CONSERVE');
    sampling_m = 0; % dummy value needed as output
    
else

    wavefront = prop_lens(wavefront, f_OAP8);
    wavefront = prop_propagate(wavefront, f_OAP8);
    [wavefront, sampling_m] = prop_end(wavefront, 'noabs');
    mag = (pupil_ratio / final_sampling_lam0) * (lambda_m / lambda0_m);
    wavefront = prop_magnify(wavefront, mag, 'SIZE_OUT', nout, 'AMP_CONSERVE');
    
end

end


function show_efield(waveStruct, label)
    % Plot phase and amplitude at the given plane
    efield = ifftshift(waveStruct.wf);
    
    figure; imagesc(abs(efield)); axis xy equal tight; colorbar;
    title([label ' amplitude'])
    drawnow;
    
    figure; imagesc(angle(efield)); axis xy equal tight; colorbar;
    title([label ' phase'])
    drawnow;
end
