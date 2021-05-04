%DST 2 PROPER MODEL - A First Pass
%By Matthew Noyes - Optical Engineer | JPL
%December 13, 2020
%NOTES:
    %This is a first pass PROPER model of DST 2. It will not contain any
    %mirror error maps, Deformable surfaces, or shaped pupil masks. Once
    %the prescription has been generated in MATLAB based off the Zemax
    %prescription, then a new model will be created to take into account
    %these omitted components. 

function [wavefront, sampling_m] = DST2(lambda_m, gridsize, optval)
% Model Switches
xoffset    = 0; %-- Introduces wavefront tilt that shifts focus in units of lam0/d
yoffset    = 0; %-- Introduces wavefront tilt that shifts focus in units of lam0/d
use_errors  = 1.*[1 1 1 1 1 1 1 1]; %-- Switches the use of error maps on optical surfaces
use_DM_errs = 0*use_errors(1);
use_fpm     = 1; %-- Swithces the use of the Focal Plane Mask
use_LYOT    = 1; %-- Switches the use of the LYOT Stop
use_FS      = 1; %-- Switches the use of the Field Stop
use_field_stop = 1;
use_pr      = false;
use_dm1 = 0;
use_dm2 = 0;

%Adding necessary functions to path, as well as the path to the mirror surface maps

addpath(genpath('C:\Users\matnoyes\Documents\GitHub\proper-models'));
addpath(genpath('C:\Users\matnoyes\Documents\Projects\HCIT_Camilo\DST2\Maps\OAP_FITS\'))

map_dir = 'C:\Users\matnoyes\Documents\Projects\HCIT_Camilo\DST2\Maps\OAP_FITS\';

% System Variables
nact = 50;                 %-- Number of actuators across DM
nact_across_pupil = 50;    %-- Number of actuators across pupil       ##### NOT THE ACTUAL NUMBER ####
dm_xc = 24.5;              %-- Wavefront centered at corner of DM actuator (0,0 is center of 1st actuator)
dm_yc = dm_xc;
dm_sampling = 400e-6;      %-- DM actuator spacing (BMC 50)
N_pix_act = 7;             %-- The number of pixels that should be on each DM atuator
diam = 0.01544*2;          %-- Diameter of the beam at the first reflective pupil. Found from Zemax Model.
lambda_um  = 0.5;          %-- Operational wavelength
lambda_m  = lambda_um/1000/1000;  
lambda0_um = 0.5;          %-- Default reference wavelength (center of bandpass) for star offsets & field stop size
lambda0_m = lambda0_um/1000/1000; 
NbeamFull = nact_across_pupil*N_pix_act;    %-- Define sampling of pupil based on having at least N_pix_act pixels across each DM actuator | Number of pixels across DM Pupil
beam_ratio = 0.25;         %-- Chosen to yeild 4 pixels per lam/D
    gridsize = 2^ceil(log2(NbeamFull/beam_ratio));
pupil_diam_pix = nact * N_pix_act;          %-- define sampling of pupil based on having N_pix_act pixels across each DM actuator
r_field_stop = 30;         %-- Field stop radius in m
final_sampling_lam0 = 1;   %-- sampling at final image plane in lam0/D
normLyotDiam = 0.95;       %-- Lyot stop outer diameter normalized to the beam diameter 
MAX_n = 0.040036173618850; %-- peak value of system WITHOUT errors, or FPM. Emperical. Used to generate Normalized Intensity Plots
MAX_y = 0.035513546229604; %-- peak value of system WITH errors, but no FPM. Emperical. Used to generate Normalized Intensity Plots
vortexCharge = 6;          %-- charge of the vortex focal plane mask
zindex = 0;                %-- vector of Zernike indices (Using Noll's sequential zernike indices)
zval = 0;      	           %-- vector of Zernike coefficients (unobscured RMS wavefront in meters)

nout = 300;


% Override default optval values if passed at the beginning
if(exist('optval','var')==1)
    if ( isfield(optval,'lam0') );                  lambda0_um = optval.lam0;  end
    if ( isfield(optval,'lambda0_um') );            lambda0_um = optval.lambda0_um;  end
%     if ( isfield(optval,'use_errors') );            use_errors = optval.use_errors;  end
    if ( isfield(optval,'zindex') );                zindex = optval.zindex;  end
    if ( isfield(optval,'zval') );                  zval = optval.zval;  end
    if ( isfield(optval,'xoffset') );               xoffset = optval.xoffset;  end
    if ( isfield(optval,'yoffset') );               yoffset = optval.yoffset;  end
    if ( isfield(optval,'use_dm1') );               use_dm1 = optval.use_dm1;  end
    if ( isfield(optval,'dm1') );                   dm1 = optval.dm1;  end
    if ( isfield(optval,'use_dm2') );               use_dm2 = optval.use_dm2;  end
    if ( isfield(optval,'dm2') );                   dm2 = optval.dm2;  end
    if ( isfield(optval,'use_fpm') );               use_fpm = optval.use_fpm;  end
    if ( isfield(optval,'use_lyot_stop') );         use_lyot_stop = optval.use_lyot_stop;  end
    if ( isfield(optval,'use_field_stop') );        use_field_stop = optval.use_field_stop;  end
    if ( isfield(optval,'field_stop_radius') );     field_stop_radius = optval.field_stop_radius;  end
    if ( isfield(optval,'final_sampling_lam0') );   final_sampling_lam0 = optval.final_sampling_lam0;  end
    if ( isfield(optval,'nout') );                  nout = optval.nout;  end
    if ( isfield(optval,'normLyotDiam') );          normLyotDiam = optval.normLyotDiam;  end
    if ( isfield(optval,'vortexCharge') );          vortexCharge = optval.vortexCharge;  end
    if ( isfield(optval,'map_dir') );               map_dir = optval.map_dir;  end
    if ( isfield(optval,'pupil_diam_pix') );        pupil_diam_pix = optval.pupil_diam_pix;  end
    if ( isfield(optval,'pr_pupil_diam_pix') );      pr_pupil_diam_pix = optval.pr_pupil_diam_pix;  end
    if ( isfield(optval,'use_pr') );                use_pr = optval.use_pr;  end
end
pr_pupil_diam_pix = pupil_diam_pix; %-- define sampling of pupil used for flattening phase with the DMs
pupil_ratio = pupil_diam_pix / double(gridsize);

%% Defining the system optical prescription | All distansces in units of [meters]
    f_OAP1 = 0.1540536;                       
d_OAP1_ReflPupil = 0.733464;       % 1st Pupil Plane
d_ReflPupil_OAP2 = 0.638094;           
    f_OAP2 = 0.774577;
    f_OAP3 = 0.464344;
d_OAP3_DM1 = pupilDist(f_OAP2,f_OAP3,d_ReflPupil_OAP2);         % 2nd Pupil Plane
d_DM1_DM2 = 0.3;
d_DM2_OAP4 = 0.538109;
    f_OAP4 = 0.641891;             % 2nd Focal Plane
    f_OAP5 = 0.774577;
d_OAP5_LYOT = pupilDist(f_OAP4,f_OAP5,d_DM1_DM2+d_DM2_OAP4);    % 3rd Pupil Plane
d_LYOT_OAP6 = 0.790014;
    f_OAP6 = 0.774577;             % 3rd Pupil Plane
    f_OAP7 = 0.156886;
d_OAP7_pupilbfimgr = pupilDist(f_OAP6,f_OAP7,d_LYOT_OAP6);      % 4th & last Pupil Plane
d_pupilbfimgr_OAP8 = 0.477237;             
    f_OAP8 = 0.475487;             % Focal length of OAP 6

%% Here begins the propogations thru the above prescription
wavefront = prop_begin(diam, lambda_m, gridsize, beam_ratio); 
wavefront = prop_circular_aperture(wavefront, diam/2);

if( (xoffset ~= 0) || (yoffset ~= 0) )
	%-- star X,Y offset in lam0/D
	xoffset_lam = xoffset * lambda0_m / lambda_m;
	yoffset_lam = yoffset * lambda0_m / lambda_m;
	u = ((0:(gridsize-1))-gridsize/2) / (pupil_diam_pix/2);  % IDL version: u = (dindgen(gridsize)-gridsize/2) / (pupil_diam_pix/2);
	xtilt = exp( 1i * pi * u * xoffset_lam );
	ytilt = exp( 1i * pi * u * yoffset_lam );
	wavefront = prop_multiply(wavefront, ytilt.' * xtilt); % IDL version: prop_multiply, wavefront, xtilt # ytilt
	clear xtilt ytilt
end

wavefront = prop_define_entrance(wavefront);

    wavefront = prop_propagate(wavefront,-d_OAP1_ReflPupil);

if(use_errors(1)); wavefront = prop_errormap(wavefront, [map_dir 'OAP1_2-WFE.fits'], 'WAVEFRONT'); end

    wavefront = prop_propagate(wavefront,d_OAP1_ReflPupil); 

    wavefront = prop_propagate(wavefront,d_ReflPupil_OAP2);

if(use_errors(2)); wavefront = prop_errormap(wavefront, [map_dir '\OAP2_2-WFE.fits'], 'WAVEFRONT'); end
wavefront = prop_lens(wavefront, f_OAP2);           

    wavefront = prop_propagate(wavefront,f_OAP2); 
                
    wavefront = prop_propagate(wavefront, f_OAP3); 
  
if(use_errors(3)); wavefront = prop_errormap(wavefront, [map_dir '\OAP3_1-WFE.fits'], 'WAVEFRONT'); end 
wavefront = prop_lens(wavefront, f_OAP3);

    wavefront = prop_propagate(wavefront, d_OAP3_DM1); 
    
if(use_dm1);  wavefront = propcustom_dm(wavefront, dm1, dm_xc, dm_yc, dm_sampling, 'inf_file', 'influence_BMC_2kDM_400micron_res10.fits'); end
if(use_DM_errs); wavefront = prop_errormap(wavefront, [map_dir '\BMC50_DM_WFE.fits'], 'wavefront'); end

    wavefront = prop_propagate(wavefront, d_DM1_DM2);
    
if(use_dm2);  wavefront = propcustom_dm(wavefront, dm2, dm_xc, dm_yc, dm_sampling, 'inf_file', 'influence_BMC_2kDM_400micron_res10.fits'); end
if(use_DM_errs); wavefront = prop_errormap(wavefront, [map_dir '\BMC50_DM_WFE.fits'], 'wavefront'); end

    wavefront = prop_propagate(wavefront, d_DM2_OAP4);  
 
if(use_errors(4)); wavefront = prop_errormap(wavefront, [map_dir '\OAP4_1-WFE.fits'], 'WAVEFRONT'); end 
wavefront = prop_lens(wavefront, f_OAP4);        

    wavefront = prop_propagate(wavefront, f_OAP4);
    
if(use_pr == false)
    if(use_fpm)
            apRad = pupil_diam_pix/2.;
            useGPU = false;
            inVal = 1;      %-- found empirically in lamd/D
            outVal = 5;     %-- found empirically in lamd/D

            % 1) IFFT to previous pupil
            % 2) Use propcustom_mft_Pup2Vortex2Pup() to go to Lyot plane
            % 3) IFFT to FPM's focal plane  
            EpupPre = ifftshift(ifft2(wavefront.wf))*gridsize; % wavefront.wf is already fftshifted
            EpupPost = propcustom_mft_Pup2Vortex2Pup(EpupPre, vortexCharge, apRad, inVal, outVal, useGPU);
            wavefront.wf = ifft2(fftshift(EpupPost))*gridsize;
    end

        wavefront = prop_propagate(wavefront, f_OAP5); 

    if(use_errors(5)); wavefront = prop_errormap(wavefront, [map_dir '\OAP5_1-WFE.fits'], 'WAVEFRONT'); end 
    wavefront = prop_lens(wavefront, f_OAP5); 

        wavefront = prop_propagate(wavefront, d_OAP5_LYOT);

    if(use_LYOT); wavefront = prop_circular_aperture(wavefront, normLyotDiam, 'NORM'); end 

        wavefront = prop_propagate(wavefront, d_LYOT_OAP6);

    if(use_errors(6)); wavefront = prop_errormap(wavefront, [map_dir '\OAP6_2-WFE.fits'], 'WAVEFRONT'); end 
    wavefront = prop_lens(wavefront, f_OAP6);

        wavefront = prop_propagate(wavefront,f_OAP6);

    if(use_field_stop)
        r_stop = r_field_stop * lambda0_m / lambda_m;
        wavefront = prop_circular_aperture(wavefront, r_stop/pupil_ratio*prop_get_sampling(wavefront));
    end
    
    Efs = ifftshift((wavefront.wf)); % wavefront.wf is already fftshifted
    fitswrite(abs(Efs),  'C:\Users\matnoyes\Documents\Projects\HCIT_Camilo\DST2\DST2_PROPER_Models\maps\Efs_abs_init_matlab.fits');
    fitswrite(angle(Efs),'C:\Users\matnoyes\Documents\Projects\HCIT_Camilo\DST2\DST2_PROPER_Models\maps\Efs_angle_init_matlab.fits');    

        wavefront = prop_propagate(wavefront,f_OAP7); 

    if(use_errors(7)); wavefront = prop_errormap(wavefront, [map_dir '\OAP7_2-WFE.fits'], 'WAVEFRONT'); end 
    wavefront = prop_lens(wavefront, f_OAP7);

        wavefront = prop_propagate(wavefront,d_OAP7_pupilbfimgr);

        wavefront = prop_propagate(wavefront,d_pupilbfimgr_OAP8); 

    if(use_errors(8)); wavefront = prop_errormap(wavefront, [map_dir '\OAP8_2-WFE.fits'], 'WAVEFRONT'); end 
    wavefront = prop_lens(wavefront, f_OAP8);

        wavefront = prop_propagate(wavefront,f_OAP8);   %Finally, Propagate to the final focal plane.
        
    [wavefront, sampling_m] = prop_end(wavefront, 'noabs');
    
    mag = (pupil_ratio / final_sampling_lam0) * (lambda_m / lambda0_m);
    wavefront = prop_magnify(wavefront, mag, 'SIZE_OUT', nout, 'AMP_CONSERVE');
    figure(1); imagesc(abs(wavefront))

 else % Get E-field at pupil preceding the FPM
     EpupBeforeFPM = ifftshift(ifft2(wavefront.wf))*gridsize; % IFFT to previous pupil
        
     mag = pr_pupil_diam_pix/pupil_diam_pix;
     wavefront = prop_magnify(EpupBeforeFPM, mag, 'SIZE_OUT', 2*ceil(0.5*mag*gridsize), 'AMP_CONSERVE');
     sampling_m = 0; %-- dummy value
end

end