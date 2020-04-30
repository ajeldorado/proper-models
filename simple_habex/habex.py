# Copyright 2020, by the California Institute of Technology. ALL RIGHTS
# RESERVED. United States Government Sponsorship acknowledged. Any
# commercial use must be negotiated with the Office of Technology Transfer
# at the California Institute of Technology.
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Name:
#	habex.m
#
# Purpose:
#	Representation of the Habex telescope and coronagraph. To be called using 
#	the PROPER library procedure "proper.prop_run".
#
# Inputs:
#	lambda_m
#	    The wavelength of propagation in meters (note that the wavelength is provided
#	    to proper.prop_run in microns and is converted to meters in there).
#	gridsize
#	    Size of the computational grid (gridsize by gridsize elements).  Must be
#	    a power of 2. 
#
# Outputs:
#	wavefront
#	    Variable in which the computed E-field at the final image plane is returned.
#	    The field is sampled by "final_sampling_lam0" lambda_m/D over "nout" by "nout"
#	    pixels.
#	sampling_m
#	    The sampling at the final image plane in meters per pixel
#
# Optional keywords or switches:
#	optval
#	    (Optional) Structure whose fields are values 
#	    that are passed to the prescription for use as the prescription desires.
#
# Revision history:
#	Written by John Krist (Jet Propulsion Laboratory, California Inst. Technology), January 2020
#   Translated to Python by A.J. Riggs (JPL, CIT), February 2020. Added an option,
#   use_pr, to retrieve the E-field at the pupil before the focal plane mask. Also
#   added the vortex as the focal plane mask.
##----------------------------------------------------------------------------------  

import numpy as np
#import matplotlib.pyplot as plt # For Debugging
#from astropy.io import fits # For Debugging
import proper # Use v3.2 or higher
import falco # FALCO needed for propagation to/from vortex

def habex(lambda_m, gridsize, PASSVALUE={'dummy':0}):

    nact = 64;                       #-- number of actuators across DM
    nact_across_pupil = 62;		#-- number of actuators across pupil
    dm_xc = 31.5;                   #-- wavefront centered at corner of DM actuator (0,0 is center of 1st actuator)
    dm_yc = 31.5;
    dm_sampling = 0.4e-3;            #-- DM actuator spacing (BMC) 

    #-- default settings (override with optval)
    map_dir = '../maps/';	#-- directory containing optical surface error maps

    lambda0_um = 0.5;	#-- default reference wavelength (center of bandpass) for star offsets & field stop size
    use_errors = 1;		#-- 1 = use optical surface errors, 0 = none
    zindex = np.array([0,])		#-- vector of Zernike indices (Noll ordered)
    zval = np.array([0,])		#-- vector of Zernike coefficients (unobscured RMS wavefront in meters)
    xoffset = 0;		#-- star X offset in lambda0/D units (must then provide lambda0_um)
    yoffset = 0;		#-- star Y offset in lambda0/D units
    use_dm1 = 0;		#-- use DM1 (if non-zero, must then provide pokes (meters) in "dm1" array)
    use_dm2 = 0;		#-- use DM2 (if non-zero, must then provide pokes (meters) in "dm2" array)
    use_fpm = 1;		#-- use focal plane mask (0 = no FPM)
    use_lyot_stop = 1;	#-- use Lyot stop (0 = no stop)
    use_field_stop = 1;	#-- use field stop (0 = no stop)
    field_stop_radius = 25.0;   #-- field stop radius in lam0/D
    final_sampling_lam0 = 0.2; #-- sampling at final image plane in lam0/D
    nout = 300;		#-- output field size (nout x nout pixels)
    normLyotDiam = 0.95;    #-- Lyot stop outer diameter normalized to the beam diameter
    vortexCharge = 6;   #-- charge of the vortex focal plane mask
    pupil_diam_pix = nact_across_pupil * 7 	#-- define sampling of pupil based on having 7 pixels across each DM actuator
    pr_pupil_diam_pix = pupil_diam_pix; #-- define sampling of pupil used for flattening phase with the DMs
    use_pr = False #-- whether to return a fake phase retrieval of the pupil rather than the focal plane

    #-- override defaults using values passed using optval structure
    if 'PASSVALUE' in locals():
        if 'lam0' in PASSVALUE: lamba0_um = PASSVALUE['lam0']
        if 'lambda0_um' in PASSVALUE: lambda0_um = PASSVALUE['lambda0_um']
        if 'use_errors' in PASSVALUE: use_errors = PASSVALUE['use_errors']
        if 'zindex' in PASSVALUE: zindex = PASSVALUE['zindex']
        if 'zval' in PASSVALUE: zval = PASSVALUE['zval']
        if 'xoffset' in PASSVALUE: xoffset = PASSVALUE['xoffset']
        if 'yoffset' in PASSVALUE: yoffset = PASSVALUE['yoffset']
        if 'use_dm1' in PASSVALUE: use_dm1 = PASSVALUE['use_dm1']
        if 'dm1' in PASSVALUE: dm1 = PASSVALUE['dm1']
        if 'use_dm2' in PASSVALUE: use_dm2 = PASSVALUE['use_dm2']
        if 'dm2' in PASSVALUE: dm2 = PASSVALUE['dm2']
        if 'use_fpm' in PASSVALUE: use_fpm = PASSVALUE['use_fpm']
        if 'use_lyot_stop' in PASSVALUE: use_lyot_stop = PASSVALUE['use_lyot_stop']
        if 'use_field_stop' in PASSVALUE: use_field_stop = PASSVALUE['use_field_stop']
        if 'field_stop_radius' in PASSVALUE: field_stop_radius = PASSVALUE['field_stop_radius']
        if 'final_sampling_lam0' in PASSVALUE: final_sampling_lam0 = PASSVALUE['final_sampling_lam0']
        if 'nout' in PASSVALUE: nout = PASSVALUE['nout']
        if 'normLyotDiam' in PASSVALUE: normLyotDiam = PASSVALUE['normLyotDiam']
        if 'vortexCharge' in PASSVALUE: vortexCharge = PASSVALUE['vortexCharge']
        if 'map_dir' in PASSVALUE: map_dir = PASSVALUE['map_dir']
        if 'pupil_diam_pix' in PASSVALUE: pupil_diam_pix = PASSVALUE['pupil_diam_pix']
        if 'pr_pupil_diam_pix' in PASSVALUE: pr_pupil_diam_pix = PASSVALUE['pr_pupil_diam_pix']
        if 'use_pr' in PASSVALUE: use_pr = PASSVALUE['use_pr']

    # Convert 0 and 1 to False and True
    use_errors = bool(use_errors)
    use_dm1 = bool(use_dm1)
    use_dm2 = bool(use_dm2)
    use_fpm = bool(use_fpm)
    use_lyot_stop = bool(use_lyot_stop)
    use_field_stop = bool(use_field_stop)
    use_pr = bool(use_pr)
    
    if(np.isscalar(zindex)):
        zindex = np.asarray((zindex,))
    else: # Check if iterable. If not, then make an array containing 0
        try:
            temp = zindex[0]
        except:
            zindex = np.array([0])          
    
    lambda0_m = lambda0_um * 1.0e-6;
    pupil_ratio = pupil_diam_pix / float(gridsize)

    #-- define optical prescription (distances, focal lengths)
    diam = 4.00;
    r_pri = 19.8;
    h_pri = 2.5;
    z_pri = h_pri**2 / (2*r_pri)
    fl_pri = np.sqrt(h_pri**2 + (r_pri/2-z_pri)**2)	#-- effective focal length of primary as a pure parabola
    d_pri_sec = 9.172532289071727;
    d_focus_sec = fl_pri - d_pri_sec;
    d_sec_focus = 7.979857207574376844;
    fl_sec = 1 / (1/d_sec_focus - 1/d_focus_sec)
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
    d_qwp_m6 = fl_m6 - (d_dm1_dm2 + d_dm2_qwp)
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


    wavefront = proper.prop_begin(diam, lambda_m, gridsize, pupil_diam_pix/gridsize)
    proper.prop_circular_aperture(wavefront, diam/2)
    if not zindex[0] == 0:   
        proper.prop_zernikes(wavefront, zindex, zval) 	#-- optionally add Zernikes

    if( (xoffset != 0) or (yoffset != 0) ):
        #-- star X,Y offset in lam0/D
        xoffset_lam = xoffset * lambda0_m / lambda_m;
        yoffset_lam = yoffset * lambda0_m / lambda_m;
        u = (np.arange(gridsize)-gridsize/2.) / (pupil_diam_pix/2.)  # IDL version: u = (dindgen(gridsize)-gridsize/2) / (pupil_diam_pix/2)
        xtilt = np.exp( 1j * np.pi * u * xoffset_lam )
        ytilt = np.exp( 1j * np.pi * u * yoffset_lam )
        proper.prop_multiply(wavefront, ytilt.reshape((gridsize, 1)) @ xtilt.reshape((1, gridsize))) # IDL version: proper.prop_multiply, wavefront, xtilt # ytilt

    if(use_errors):  proper.prop_errormap(wavefront, map_dir+'habex_cycle1_PRIMARY_phase_error.fits', WAVEFRONT=True)
    proper.prop_lens(wavefront, fl_pri)
    proper.prop_define_entrance(wavefront)

    proper.prop_propagate(wavefront, d_pri_sec, 'secondary')
    if(use_errors):  proper.prop_errormap(wavefront, map_dir+'habex_cycle1_SECONDARY_phase_error.fits', WAVEFRONT=True)
    proper.prop_lens(wavefront, fl_sec)

    proper.prop_propagate(wavefront, d_sec_m3, 'M3')
    if(use_errors):  proper.prop_errormap(wavefront, map_dir+'habex_cycle1_M3_phase_error.fits', WAVEFRONT=True)
    proper.prop_lens(wavefront, fl_m3)

    proper.prop_propagate(wavefront, d_m3_fold, 'fold')
    if(use_errors):  proper.prop_errormap(wavefront, map_dir+'habex_cycle1_FOLD1_phase_error.fits', WAVEFRONT=True)
    
    proper.prop_propagate(wavefront, d_fold_fsm, 'FSM')	#-- pupil at fast steering mirror (interface with telescope)
    if(use_errors):  proper.prop_errormap(wavefront, map_dir+'habex_cycle1_FSM_phase_error.fits', WAVEFRONT=True)

    proper.prop_propagate(wavefront, d_fsm_dichroic, 'dichroic')
    if(use_errors):  proper.prop_errormap(wavefront, map_dir+'habex_cycle1_DICHROIC_phase_error.fits', WAVEFRONT=True)

    proper.prop_propagate(wavefront, d_dichroic_m4, 'M4')
    if(use_errors):  proper.prop_errormap(wavefront, map_dir+'habex_cycle1_M4_phase_error.fits', WAVEFRONT=True)
    proper.prop_lens(wavefront, fl_m4)

    proper.prop_propagate(wavefront, d_m4_m5, 'M5')
    if(use_errors):  proper.prop_errormap(wavefront, map_dir+'habex_cycle1_M5_phase_error.fits', WAVEFRONT=True)
    proper.prop_lens(wavefront, fl_m5)

    proper.prop_propagate(wavefront, d_m5_dm1, 'DM1')
    if(use_dm1):  proper.prop_dm(wavefront, dm1, dm_xc, dm_yc, dm_sampling)
    if(use_errors):  proper.prop_errormap(wavefront, map_dir+'habex_cycle1_DM1_phase_error.fits', WAVEFRONT=True)
    
    proper.prop_propagate(wavefront, d_dm1_dm2, 'DM2')
    if(use_dm2):  proper.prop_dm(wavefront, dm2, dm_xc, dm_yc, dm_sampling)
    if(use_errors):  proper.prop_errormap(wavefront, map_dir+'habex_cycle1_DM2_phase_error.fits', WAVEFRONT=True)

    proper.prop_propagate(wavefront, d_dm2_qwp, 'QWP')	#-- quarter-wave plate
    if(use_errors):  proper.prop_errormap(wavefront, map_dir+'habex_cycle1_QWP1_phase_error.fits', WAVEFRONT=True)

    proper.prop_propagate(wavefront, d_qwp_m6, 'M6')
    if(use_errors):  proper.prop_errormap(wavefront, map_dir+'habex_cycle1_M6_phase_error.fits', WAVEFRONT=True)
    proper.prop_lens(wavefront, fl_m6)

    proper.prop_propagate(wavefront, d_m6_fpm)

    if not use_pr:

#        if use_fpm:  
#            fpm = proper.prop_8th_order_mask(wavefront, 4.0, CIRCULAR=True)   #--Band-limited mask
        
        if use_fpm:

            apRad = pupil_diam_pix/2.
            inVal = 0.3    #-- found empirically
            outVal = 5     #-- found empirically

            # 1) IFFT to previous pupil from FPM's plane
            # 2) Use propcustom_mft_Pup2Vortex2Pup() to go to Lyot plane
            # 3) IFFT to FPM's focal plane  
            EpupPre = np.fft.ifftshift(np.fft.ifft2(wavefront.wfarr))*gridsize # wavefront.wf is already fftshifted
            EpupPost = falco.prop.mft_p2v2p(EpupPre, vortexCharge, apRad, inVal, outVal)
            wavefront.wfarr = np.fft.ifft2(np.fft.fftshift(EpupPost))*gridsize

        proper.prop_propagate(wavefront, d_fpm_m7, 'M7')
        if(use_errors):  proper.prop_errormap(wavefront, map_dir+'habex_cycle1_M7_phase_error.fits', WAVEFRONT=True)
        proper.prop_lens(wavefront, fl_m7)

        proper.prop_propagate(wavefront, d_m7_lyotstop, 'Lyot stop')

        if(use_errors):  proper.prop_errormap(wavefront, map_dir+'habex_cycle1_QWP2_phase_error.fits', WAVEFRONT=True)
        if(use_lyot_stop):  proper.prop_circular_aperture(wavefront, normLyotDiam, NORM=True)

        proper.prop_propagate(wavefront, d_lyotstop_m8, 'M8')
        if(use_errors):  proper.prop_errormap(wavefront, map_dir+'habex_cycle1_M8_phase_error.fits', WAVEFRONT=True)
        proper.prop_lens(wavefront, fl_m8)

        proper.prop_propagate(wavefront, proper.prop_get_distancetofocus(wavefront), 'field stop')

        if(use_field_stop):
            r_stop = field_stop_radius * lambda0_m / lambda_m;
            proper.prop_circular_aperture(wavefront, r_stop/pupil_ratio*proper.prop_get_sampling(wavefront))

        proper.prop_propagate(wavefront, d_fieldstop_m9, 'M9')
        if(use_errors):  proper.prop_errormap(wavefront, map_dir+'habex_cycle1_M9_phase_error.fits', WAVEFRONT=True)
        proper.prop_lens(wavefront, fl_m9)

        proper.prop_propagate(wavefront, d_m9_filter, 'filter')
        if(use_errors):  proper.prop_errormap(wavefront, map_dir+'habex_cycle1_FILTER_phase_error.fits', WAVEFRONT=True)

        proper.prop_propagate(wavefront, d_filter_m10, 'M10')
        if(use_errors):  proper.prop_errormap(wavefront, map_dir+'habex_cycle1_M10_phase_error.fits', WAVEFRONT=True)

        proper.prop_lens(wavefront, fl_m10)

        proper.prop_propagate(wavefront, proper.prop_get_distancetofocus(wavefront), 'CCD')

        [wavefront, sampling_m] = proper.prop_end(wavefront, NOABS=True)

        #-- rescale to "final_sampling_lam0" lam0/D per pixel
        mag = (pupil_ratio / final_sampling_lam0) * (lambda_m / lambda0_m)
        wavefront = proper.prop_magnify(wavefront, mag, nout, AMP_CONSERVE=True)

    else: #-- Perfect-knowledge phase retrieval

            # 1) IFFT to previous pupil
            EpupBeforeFPM = np.fft.ifftshift(np.fft.ifft2(wavefront.wfarr))*gridsize

            mag = pr_pupil_diam_pix/pupil_diam_pix;
            wavefront = proper.prop_magnify(EpupBeforeFPM, mag, 2*np.ceil(0.5*mag*gridsize), AMP_CONSERVE=True)
            sampling_m = 0 #-- dummy value

    return wavefront, sampling_m
