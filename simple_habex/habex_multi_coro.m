% Copyright 2020, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Name:
%	habex.m
%
% Purpose:
%	Representation of the Habex telescope and coronagraph. To be called using 
%	the PROPER library procedure "prop_run".
%
% Inputs:
%	lambda_m
%	    The wavelength of propagation in meters (note that the wavelength is provided
%	    to prop_run in microns and is converted to meters in there).
%	gridsize
%	    Size of the computational grid (gridsize by gridsize elements).  Must be
%	    a power of 2. 
%
% Outputs:
%	wavefront
%	    Variable in which the computed E-field at the final image plane is returned.
%	    The field is sampled by "final_sampling_lam0" lambda_m/D over "nout" by "nout"
%	    pixels.
%	sampling_m
%	    The sampling at the final image plane in meters per pixel
%
% Optional keywords or switches:
%	optval
%	    (Optional) Structure whose fields are values 
%	    that are passed to the prescription for use as the prescription desires.
%
% Revision history:
%	Written by John Krist (Jet Propulsion Laboratory, California Inst. Technology), January 2020
%   Translated to MATLAB by A.J. Riggs (JPL, CIT), February 2020. Added an option,
%   use_pr, to retrieve the E-field at the pupil before the focal plane mask. Also 
%   added a vortex as the focal plane mask.
%%----------------------------------------------------------------------------------  

function [wavefront, sampling_m] = habex_multi_coro(lambda_m, gridsize, optval)

nact = 64;                       %-- number of actuators across DM
nact_across_pupil = 62;		%-- number of actuators across pupil       ##### NOT THE ACTUAL NUMBER ####
dm_xc = 31.5;                   %-- wavefront centered at corner of DM actuator (0,0 is center of 1st actuator)
dm_yc = 31.5;
dm_sampling = 0.4e-3;            %-- DM actuator spacing (BMC) 
% pupil_diam_pix = nact_across_pupil * 7.0; 	%-- define sampling of pupil based on having 7 pixels across each DM actuator
% pupil_ratio = pupil_diam_pix / double(gridsize);

%-- default settings (override with optval)
map_dir = '../maps/';	%-- directory containing optical surface error maps

cor_type = 'vortex';%-- default coronagraph type
lambda0_um = 0.5;	%-- default reference wavelength (center of bandpass) for star offsets & field stop size
use_errors = 1;		%-- 1 = use optical surface errors, 0 = none
zindex = 0;		%-- vector of Zernike indices (Noll ordered)
zval = 0;		%-- vector of Zernike coefficients (unobscured RMS wavefront in meters)
xoffset = 0;		%-- star X offset in lambda0/D units (must then provide lambda0_um)
yoffset = 0;		%-- star Y offset in lambda0/D units
use_dm1 = 0;		%-- use DM1 (if non-zero, must then provide pokes (meters) in "dm1" array)
use_dm2 = 0;		%-- use DM2 (if non-zero, must then provide pokes (meters) in "dm2" array)
use_fpm = 1;		%-- use focal plane mask (0 = no FPM)
use_lyot_stop = 1;	%-- use Lyot stop (0 = no stop)
use_field_stop = 1;	%-- use field stop (0 = no stop)
field_stop_radius = 25.0;   %-- field stop radius in lam0/D
final_sampling_lam0 = 0.2; %-- sampling at final image plane in lam0/D
nout = 300;		%-- output field size (nout x nout pixels)
normLyotDiam = 0.95;    %-- Lyot stop outer diameter normalized to the beam diameter
vortexCharge = 6;   %-- charge of the vortex focal plane mask
pupil_diam_pix = nact_across_pupil * 7.0; 	%-- define sampling of pupil based on having 7 pixels across each DM actuator
use_pr = false; %-- whether to return a fake phase retrieval of the pupil rather than the focal plane
use_hlc_dm_patterns = 0;	% use Dwight-generated HLC default DM wavefront patterns? 1 or 0
%-- override defaults using values passed using optval structure

if(exist('optval','var')==1)
    if ( isfield(optval,'cor_type') );  cor_type = optval.cor_type;  end
    if ( isfield(optval,'lam0') );  lambda0_um = optval.lam0;  end
    if ( isfield(optval,'lambda0_um') );  lambda0_um = optval.lambda0_um;  end
    if ( isfield(optval,'use_errors') );  use_errors = optval.use_errors;  end
    if ( isfield(optval,'zindex') );  zindex = optval.zindex;  end
    if ( isfield(optval,'zval') );  zval = optval.zval;  end
    if ( isfield(optval,'xoffset') );  xoffset = optval.xoffset;  end
    if ( isfield(optval,'yoffset') );  yoffset = optval.yoffset;  end
    if ( isfield(optval,'use_dm1') );  use_dm1 = optval.use_dm1;  end
    if ( isfield(optval,'dm1') );  dm1 = optval.dm1;  end
    if ( isfield(optval,'use_dm2') );  use_dm2 = optval.use_dm2;  end
    if ( isfield(optval,'dm2') );  dm2 = optval.dm2;  end
    if ( isfield(optval,'use_fpm') );  use_fpm = optval.use_fpm;  end
    if ( isfield(optval,'use_lyot_stop') );  use_lyot_stop = optval.use_lyot_stop;  end
    if ( isfield(optval,'use_field_stop') );  use_field_stop = optval.use_field_stop;  end
    if ( isfield(optval,'field_stop_radius') );  field_stop_radius = optval.field_stop_radius;  end
    if ( isfield(optval,'final_sampling_lam0') );  final_sampling_lam0 = optval.final_sampling_lam0;  end
    if ( isfield(optval,'nout') );  nout = optval.nout;  end
    if ( isfield(optval,'normLyotDiam') );  normLyotDiam = optval.normLyotDiam;  end
    if ( isfield(optval,'vortexCharge') );  vortexCharge = optval.vortexCharge;  end
    if ( isfield(optval,'map_dir') );  map_dir = optval.map_dir;  end
    if ( isfield(optval,'pupil_diam_pix') );  pupil_diam_pix = optval.pupil_diam_pix;  end
    if ( isfield(optval,'use_pr') );  use_pr = optval.use_pr;  end
    if ( isfield(optval,'use_hlc_dm_patterns') );   use_hlc_dm_patterns = optval.use_hlc_dm_patterns; end
end

pr_pupil_diam_pix = pupil_diam_pix; %-- define sampling of pupil used for flattening phase with the DMs
if(exist('optval','var')==1)
    if ( isfield(optval,'pr_pupil_diam_pix') );  pr_pupil_diam_pix = optval.pr_pupil_diam_pix;  end
end


lambda0_m = lambda0_um * 1.0e-6;
pupil_ratio = pupil_diam_pix / double(gridsize);

is_vc = false;
is_hlc = false;
is_blc = false;
if  strcmp(cor_type,'hlc')
    is_hlc = true;
    mask_dir = [optval.mask_dir filesep];
    hlc_factor = 1.06318242311;
    pupil_diam_pix = 314.581 / hlc_factor;
    pr_pupil_diam_pix = pupil_diam_pix;
    pupil_ratio = pupil_diam_pix / double(gridsize);
    
    dm1_ytilt_deg = 9.65;
    dm2_ytilt_deg = 9.65;
    
%     dm_sampling = 0.4e-3 * 0.94942161; %(298.67/314.581);
    
    pupil_fn = [mask_dir  optval.pupil_fn];
    lyot_stop_fn = [mask_dir  optval.lyot_stop_fn];
    %lambda0_m = 0.575e-6;
%     nlams = 9;              % number of occ trans lambda provided
%     bw = 0.10;
    nlams = optval.nlams;
    bw = optval.bw;
    fpm_aoi = optval.fpm_aoi;
    fpm_pol = optval.fpm_pol;
    lam_occ = linspace(1-bw/2, 1+bw/2, nlams)*lambda0_m; % wavelengths at which occ trans provided
    wlam = find( round(1e13*lambda_m) == round(1e13*lam_occ) ); 	% find exactly matching FPM wavelength
    occulter_file_r = [mask_dir 'run819_roman_occ_lam' num2str(lam_occ(wlam),12) sprintf('theta%spol%s_real.fits', fpm_aoi, fpm_pol)];
    occulter_file_i = [mask_dir 'run819_roman_occ_lam' num2str(lam_occ(wlam),12) sprintf('theta%spol%s_imag.fits', fpm_aoi, fpm_pol)];
    n_default = 1024;	% gridsize in non-critical areas
    if  use_fpm;    n_to_fpm = 4096; else; n_to_fpm = 1024; end
    n_from_lyotstop = 1024;
elseif strcmpi(cor_type, 'vortex')
    is_vc = true;
elseif strcmpi(cor_type, 'blc') || strcmpi(cor_type, 'band-limited')
    is_blc = true;
end

%-- define optical prescription (distances, focal lengths)

diam = 4.00;
  r_pri = 19.8;
  h_pri = 2.5;
  z_pri = h_pri^2 / (2*r_pri);
  fl_pri = sqrt(h_pri^2 + (r_pri/2-z_pri)^2);	%-- effective focal length of primary as a pure parabola
d_pri_sec = 9.172532289071727;
  d_focus_sec = fl_pri - d_pri_sec;
  d_sec_focus = 7.979857207574376844;
  fl_sec = 1 / (1/d_sec_focus - 1/d_focus_sec);
d_sec_m3 = 9.076690863872008;
  fl_m3 = d_sec_m3 - d_sec_focus;
d_m3_fold = 0.654597300210990;
d_fold_fsm = 0.577743120280288;
d_fsm_dichroic = 0.1950;
d_dichroic_m4 = 0.450;
  fl_m4 = 0.5075;
d_m4_m5 = 0.762954002022743;
  fl_m5 = d_m4_m5 - fl_m4;
d_m5_dm1 = 0.220615776458241;
d_dm1_dm2 = 0.32;
d_dm2_qwp = 0.32 + 0.157485214529470;
  fl_m6 = 1.029143136045496931;
d_qwp_m6 = fl_m6 - (d_dm1_dm2 + d_dm2_qwp);
d_m6_fpm = fl_m6;
d_fpm_m7 = 0.255580492381039;
  fl_m7 = d_fpm_m7;
d_m7_lyotstop = fl_m7; 
d_lyotstop_m8 = 0.2536;	
  fl_m8 = d_lyotstop_m8;
d_m8_fieldstop = fl_m8;
d_fieldstop_m9 = d_m8_fieldstop;
  fl_m9 = d_fieldstop_m9;
d_m9_filter = 0.296399999724129;
d_filter_m10 = 0.462615469378302;
  fl_m10 = 0.503971038519431261;
d_m10_ccd = fl_m10;	


wavefront = prop_begin(diam, lambda_m, gridsize, 'beam_diam_fraction', pupil_diam_pix/gridsize);
if is_vc || is_blc
    wavefront = prop_circular_aperture(wavefront, diam/2);
elseif is_hlc
   pupil = pad_crop(fitsread(pupil_fn), gridsize);
   wavefront = prop_multiply(wavefront, pupil);
else
    error('Coronagraph type not recognized.')
end
if( zindex(1) ~= 0);  wavefront = prop_zernikes(wavefront, zindex, zval);  end	%-- optionally add Zernikes
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

if(use_errors); wavefront = prop_errormap(wavefront, [map_dir 'habex_cycle1_PRIMARY_phase_error.fits'], 'WAVEFRONT'); end
wavefront = prop_lens(wavefront, fl_pri);
wavefront = prop_define_entrance(wavefront);

wavefront = prop_propagate(wavefront, d_pri_sec, 'SURFACE_NAME', 'secondary');
if(use_errors); wavefront = prop_errormap(wavefront, [map_dir 'habex_cycle1_SECONDARY_phase_error.fits'], 'WAVEFRONT');  end
wavefront = prop_lens(wavefront, fl_sec);

wavefront = prop_propagate(wavefront, d_sec_m3, 'SURFACE_NAME', 'M3');
if(use_errors); wavefront = prop_errormap(wavefront, [map_dir 'habex_cycle1_M3_phase_error.fits'], 'WAVEFRONT');  end
wavefront = prop_lens(wavefront, fl_m3);

wavefront = prop_propagate(wavefront, d_m3_fold, 'SURFACE_NAME', 'fold');
if(use_errors); wavefront = prop_errormap(wavefront, [map_dir 'habex_cycle1_FOLD1_phase_error.fits'], 'WAVEFRONT');  end

wavefront = prop_propagate(wavefront, d_fold_fsm, 'SURFACE_NAME', 'FSM');	%-- pupil at fast steering mirror (interface with telescope)
if(use_errors); wavefront = prop_errormap(wavefront, [map_dir 'habex_cycle1_FSM_phase_error.fits'], 'WAVEFRONT');  end

wavefront = prop_propagate(wavefront, d_fsm_dichroic, 'SURFACE_NAME', 'dichroic');
if(use_errors); wavefront = prop_errormap(wavefront, [map_dir 'habex_cycle1_DICHROIC_phase_error.fits'], 'WAVEFRONT');  end

wavefront = prop_propagate(wavefront, d_dichroic_m4, 'SURFACE_NAME', 'M4');
if(use_errors); wavefront = prop_errormap(wavefront, [map_dir 'habex_cycle1_M4_phase_error.fits'], 'WAVEFRONT');  end
wavefront = prop_lens(wavefront, fl_m4);

wavefront = prop_propagate(wavefront, d_m4_m5, 'SURFACE_NAME', 'M5');
if(use_errors); wavefront = prop_errormap(wavefront, [map_dir 'habex_cycle1_M5_phase_error.fits'], 'WAVEFRONT');  end
wavefront = prop_lens(wavefront, fl_m5);

% wavefront = prop_propagate(wavefront, d_m5_dm1, 'SURFACE_NAME', 'DM1');
% if(use_dm1);  [wavefront, map1] = prop_dm(wavefront, dm1, dm_xc, dm_yc, dm_sampling, 'ytilt', 9.65);  end
% if(use_errors); wavefront = prop_errormap(wavefront, [map_dir 'habex_cycle1_DM1_phase_error.fits'], 'WAVEFRONT');  end

% dm_sampling = 0.4e-3 * 0.94942161; %(298.67/314.581);
% if(use_dm1);  [wavefront, map1] = prop_dm(wavefront, dm1, dm_xc, dm_yc, dm_sampling, 'ytilt', 9.65);  end
% Ncrop = 400;
% map1 = pad_crop(map1, Ncrop);
% figure(601); imagesc(map1); axis xy equal tight; colorbar; drawnow;
% map1true = fitsread('/Users/ajriggs/Documents/habex/run819/run819_roman_dm1wfe.fits');
% map1true = pad_crop(map1true, Ncrop);
% figure(602); imagesc(map1true/2); axis xy equal tight; colorbar; drawnow;
% figure(603); imagesc(1e9*(map1 - map1true/2)); axis xy equal tight; colorbar; set(gca, 'Fontsize', 20); drawnow;
% diff = 1e9*(map1 - map1true/2);
% fprintf('%.3f\n', sum(abs(diff(:))));

if is_hlc
    
    % to accommodate Dwight's Run 819 design
    pitchRatio = 0.4/0.9906;
    d_dm1_dm2 = 1.0 * (pitchRatio*pitchRatio); % meters
    d_qwp_m6 = fl_m6 - (d_dm1_dm2 + d_dm2_qwp);
    
    wavefront = prop_propagate(wavefront, d_m5_dm1, 'SURFACE_NAME', 'DM1');
    if(use_dm1);  wavefront = prop_dm(wavefront, dm1, dm_xc, dm_yc, dm_sampling, 'ytilt', dm1_ytilt_deg);  end
    if (is_hlc  && use_hlc_dm_patterns)
        dm1wfe = fitsread([mask_dir  optval.dm1wfe_fn]);
        wavefront = prop_add_phase(wavefront, pad_crop(dm1wfe, gridsize));
        clear dm1wfe 
        
%         wf0 = ifftshift(wavefront.wf);
%         fitswrite(abs(wf0), '/Users/ajriggs/Downloads/absEdm1.fits');
%         fitswrite(angle(wf0), '/Users/ajriggs/Downloads/angleEdm1.fits');
        
    end
    if(use_errors); wavefront = prop_errormap(wavefront, [map_dir 'habex_cycle1_DM1_phase_error.fits'], 'WAVEFRONT');  end

    wavefront = prop_propagate(wavefront, d_dm1_dm2, 'SURFACE_NAME', 'DM2');
    if(use_dm2);  wavefront = prop_dm(wavefront, dm2, dm_xc, dm_yc, dm_sampling, 'ytilt', dm2_ytilt_deg);  end
    if (is_hlc  && use_hlc_dm_patterns)
        dm2wfe = fitsread([mask_dir  optval.dm2wfe_fn]);
        wavefront = prop_add_phase(wavefront, pad_crop(dm2wfe, gridsize));
        clear dm2wfe 
    end
    if(use_errors); wavefront = prop_errormap(wavefront, [map_dir 'habex_cycle1_DM2_phase_error.fits'], 'WAVEFRONT');  end

%     wavefront = prop_propagate(wavefront, -d_dm1_dm2, 'SURFACE_NAME', 'DM1eff');
%     wf0 = ifftshift(wavefront.wf);
%     fitswrite(abs(wf0), '/Users/ajriggs/Downloads/absEdm1eff.fits');
%     fitswrite(angle(wf0), '/Users/ajriggs/Downloads/angleEdm1eff.fits');
%     disp('stop');
    
else
    
    wavefront = prop_propagate(wavefront, d_m5_dm1, 'SURFACE_NAME', 'DM1');
    if(use_dm1);  wavefront = prop_dm(wavefront, dm1, dm_xc, dm_yc, dm_sampling, 'ytilt', 9.65);  end
    if(use_errors); wavefront = prop_errormap(wavefront, [map_dir 'habex_cycle1_DM1_phase_error.fits'], 'WAVEFRONT');  end
    
    wavefront = prop_propagate(wavefront, d_dm1_dm2, 'SURFACE_NAME', 'DM2');
    if(use_dm2);  wavefront = prop_dm(wavefront, dm2, dm_xc, dm_yc, dm_sampling);  end
    if(use_errors); wavefront = prop_errormap(wavefront, [map_dir 'habex_cycle1_DM2_phase_error.fits'], 'WAVEFRONT');  end

end

wavefront = prop_propagate(wavefront, d_dm2_qwp, 'SURFACE_NAME', 'QWP');	%-- quarter-wave plate
if(use_errors); wavefront = prop_errormap(wavefront, [map_dir 'habex_cycle1_QWP1_phase_error.fits'], 'WAVEFRONT');  end

wavefront = prop_propagate(wavefront, d_qwp_m6, 'SURFACE_NAME', 'M6');
if(use_errors); wavefront = prop_errormap(wavefront, [map_dir 'habex_cycle1_M6_phase_error.fits'], 'WAVEFRONT');  end
wavefront = prop_lens(wavefront, fl_m6);

wavefront = prop_propagate(wavefront, d_m6_fpm);

if(use_pr == false)
    if(use_fpm)

        if is_vc
            apRad = pupil_diam_pix/2.;
            useGPU = false;
            inVal = 0.3;    %-- found empirically
            outVal = 5;     %-- found empirically

            % 1) IFFT to previous pupil
            % 2) Use propcustom_mft_Pup2Vortex2Pup() to go to Lyot plane
            % 3) IFFT to FPM's focal plane  
            EpupPre = ifftshift(ifft2(wavefront.wf))*gridsize; % wavefront.wf is already fftshifted
            EpupPost = propcustom_mft_Pup2Vortex2Pup(EpupPre, vortexCharge, apRad, inVal, outVal, useGPU);
            wavefront.wf = ifft2(fftshift(EpupPost))*gridsize;
            
        elseif is_hlc
%             fprintf('  %s\n  %s\n\n', occulter_file_r, occulter_file_i)
%             occ = complex(fitsread(occulter_file_r),fitsread(occulter_file_i));
            occulterUpsampleFactor = 4;
            
            occ0 = complex(fitsread(occulter_file_r),fitsread(occulter_file_i));
            if angle(occ0(1,1)) ~= 0
                occ0 = occ0.*exp(-1j*angle(occ0(1,1))); %--Standardize the phase of the masks to be 0 for the outer glass part.
            end
            occ0 = propcustom_relay(occ0, 1);
            
            if occulterUpsampleFactor > 1
                
                wf0 = ifftshift(wavefront.wf);
                % back propagate to previous pupil (IFFT)
                EpupPre = ifftshift(ifft2(wavefront.wf))*gridsize; % wavefront.wf is already fftshifted
                % zero-pad
                EpupPrePad = pad_crop(EpupPre, gridsize*occulterUpsampleFactor);
                % forward propagate (FFT)
                EfocHighRes = ifftshift(fft2(fftshift(EpupPrePad)))/(gridsize*occulterUpsampleFactor);
                % Apply higher-res occulter
                EfocHighRes = EfocHighRes.*pad_crop(occ0, (gridsize*occulterUpsampleFactor));
                % back propagate to previous pupil (IFFT)
                EpupPrePad2 = ifftshift(ifft2(fftshift(EfocHighRes)))*(gridsize*occulterUpsampleFactor);
                % crop back down to original size
                EpupPre2 = pad_crop(EpupPrePad2, gridsize);
                % forward propagate (FFT) back to focal plane
                Efoc = fft2(fftshift(EpupPre2))/gridsize;
                % put E-field back into PROPER's wavefront structure
                wavefront.wf = Efoc;
                
%                 figure(401); imagesc(log10(abs(wf0))); axis xy equal tight; colorbar;
%                 figure(402); imagesc(log10(abs(EfocHighRes))); axis xy equal tight; colorbar;
%                 figure(410); imagesc(log10(abs(fftshift(Efoc)))); axis xy equal tight; colorbar;
% %                 figure(411); imagesc(abs(wf0)-(abs(fftshift(Efoc)))); axis xy equal tight; colorbar;
                
            else
                wavefront = prop_multiply(wavefront, pad_crop(occ0, n_default, 'extrapval', occ(1, 1)));
                clear occ
            end

%             % Fourier Downsample by Factor of 2
%             occ0 = complex(fitsread(occulter_file_r),fitsread(occulter_file_i));
%             occ0 = pad_crop(occ0, 4096);
%             nocc = size(occ0, 1);
%             focc = fft2(ifftshift(occ0));
%             focc = fftshift(focc);
%             focc = pad_crop(focc, nocc/2);
%             focc = fftshift(focc);
%             occ = fftshift(ifft2(focc))/dsfac^2;
%             occ = pad_crop(occ, nocc);
%     %         occ = pad_crop(occ, nocc, 'extrapval', occ0(1,1));
%             figure(41); imagesc(abs(occ0)); axis xy equal tight; colorbar;
%             figure(42); imagesc(abs(occ)); axis xy equal tight; colorbar;
%             
%             
%             %--DEBUGGING
%             if angle(occ(1,1)) ~= 0
%                 occ = occ.*exp(-1j*angle(occ(1,1))); %--Standardize the phase of the masks to be 0 for the outer glass part.
%             end
%             
%             wavefront = prop_multiply(wavefront, pad_crop(occ, gridsize, 'extrapval', occ(1, 1)));
%             clear occ
            
        elseif is_blc
            [wavefront, fpm] = prop_8th_order_mask(wavefront, 4.0, 'CIRCULAR');
            
        else
            error('Value of cor_type not recognized.')
        end

    end

    wavefront = prop_propagate(wavefront, d_fpm_m7, 'SURFACE_NAME', 'M7');
    if(use_errors); wavefront = prop_errormap(wavefront, [map_dir 'habex_cycle1_M7_phase_error.fits'], 'WAVEFRONT');  end
    wavefront = prop_lens(wavefront, fl_m7);

    wavefront = prop_propagate(wavefront, d_m7_lyotstop, 'SURFACE_NAME', 'Lyot stop');
    
    if(use_errors); wavefront = prop_errormap(wavefront, [map_dir 'habex_cycle1_QWP2_phase_error.fits'], 'WAVEFRONT');  end
    if(use_lyot_stop)
        if is_vc || is_blc
            wavefront = prop_circular_aperture(wavefront, normLyotDiam, 'NORM');
        elseif is_hlc
           lyot = pad_crop(fitsread(lyot_stop_fn), gridsize);
           wavefront = prop_multiply(wavefront, lyot);
        else
            error('Coronagraph type not recognized.')
        end
    end

    wavefront = prop_propagate(wavefront, d_lyotstop_m8, 'SURFACE_NAME', 'M8');
    if(use_errors); wavefront = prop_errormap(wavefront, [map_dir 'habex_cycle1_M8_phase_error.fits'], 'WAVEFRONT');  end
    wavefront = prop_lens(wavefront, fl_m8);

    wavefront = prop_propagate(wavefront, prop_get_distancetofocus(wavefront), 'SURFACE_NAME', 'field stop');

    if(use_field_stop)
        r_stop = field_stop_radius * lambda0_m / lambda_m;
        wavefront = prop_circular_aperture(wavefront, r_stop/pupil_ratio*prop_get_sampling(wavefront));
    end


%     Efs = ifftshift((wavefront.wf)); % wavefront.wf is already fftshifted
%     fitswrite(abs(Efs),'/Users/ajriggs/Downloads/Efs_abs_init_matlab.fits');
%     fitswrite(angle(Efs),'/Users/ajriggs/Downloads/Efs_angle_init_matlab.fits');      
    
    wavefront = prop_propagate(wavefront, d_fieldstop_m9, 'SURFACE_NAME', 'M9');
    if(use_errors); wavefront = prop_errormap(wavefront, [map_dir 'habex_cycle1_M9_phase_error.fits'], 'WAVEFRONT');  end
    wavefront = prop_lens(wavefront, fl_m9);

    wavefront = prop_propagate(wavefront, d_m9_filter, 'SURFACE_NAME', 'filter');
    if(use_errors); wavefront = prop_errormap(wavefront, [map_dir 'habex_cycle1_FILTER_phase_error.fits'], 'WAVEFRONT');  end

    wavefront = prop_propagate(wavefront, d_filter_m10, 'SURFACE_NAME', 'M10');
    if(use_errors); wavefront = prop_errormap(wavefront, [map_dir 'habex_cycle1_M10_phase_error.fits'], 'WAVEFRONT');  end
    wavefront = prop_lens(wavefront, fl_m10);
    
    wavefront = prop_propagate(wavefront, prop_get_distancetofocus(wavefront), 'SURFACE_NAME', 'CCD');

    [wavefront, sampling_m] = prop_end(wavefront, 'noabs');

    %-- rescale to "final_sampling_lam0" lam0/D per pixel
    mag = (pupil_ratio / final_sampling_lam0) * (lambda_m / lambda0_m);
    if is_hlc
        mag = mag * hlc_factor;
    end
    wavefront = prop_magnify(wavefront, mag, 'SIZE_OUT', nout, 'AMP_CONSERVE');
    
else % Get E-field at pupil preceding the FPM
        EpupBeforeFPM = ifftshift(ifft2(wavefront.wf))*gridsize; % IFFT to previous pupil
        
        mag = pr_pupil_diam_pix / pupil_diam_pix;
        wavefront = prop_magnify(EpupBeforeFPM, mag, 'SIZE_OUT', 2*ceil(0.5*mag*gridsize), 'AMP_CONSERVE');
        sampling_m = 0; %-- dummy value
end

end
