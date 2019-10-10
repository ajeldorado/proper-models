#   Copyright 2019 California Institute of Technology
# ------------------------------------------------------------------

# Version 1.0, 3 Jan 2019, JEK
# Version 1.0a, 9 Jan 2019, JEK
#   Changes: If reading in input field using input_field_rootname, the
#            user can specify polaxis, and the appropriate file will
#            be used.  By default, polaxis = 0.  This is the only
#            time polaxis can be specified for the compact model.
#   Version 1.0c, 13 February 2019, JEK 
#      Changes: Changed HLC design to HLC-20190210; added 'fpm_axis' parameter to specify
#                 axis of HLC FPM ('s' or 'p') due to FPM tilt; changed 'spc' coronagraph
#                 option to 'cor_type' to 'spc-ifs' and changed default SPC to SPC-20190130 and
#                 updated associated parameters; added 'spc-wide' option to 'cor_type' to use
#                 SPC-20191220. Changed DM tilt to rotation around Y axis.
#   Version 1.2, 28 May 2019, JEK
#       Changed lib_dir to data_dir 
#   Vesion 1.4, 7 Oct 2019, JEK
#       Aliased spc-spec_* to spc-ifs_*

import proper
import numpy as np
from scipy.interpolate import interp1d
import math
import wfirst_phaseb_proper
from wfirst_phaseb_proper import trim, ffts, mft2

# Written by John Krist


##########################################################################
def angle( n ):
    x = np.arange( n ) - int(n)//2
    return np.arctan2( x[:,np.newaxis], x ) * (180 / np.pi)

##########################################################################
def wfirst_phaseb_compact( lambda_m, output_dim0, PASSVALUE={'dummy':0} ):

    # "output_dim" is used to specify the output dimension in pixels at the final image plane.
    # Computational grid sizes are hardcoded for each coronagraph.
    # Based on Zemax prescription "WFIRST_CGI_DI_LOWFS_Sep24_2018.zmx" by Hong Tang.

    data_dir = wfirst_phaseb_proper.data_dir
    if 'PASSVALUE' in locals():
        if 'data_dir' in PASSVALUE: data_dir = PASSVALUE['data_dir']

    cor_type = 'hlc'            # coronagraph type ('hlc', 'spc', 'none')
    input_field_rootname = ''   # rootname of files containing aberrated pupil
    polaxis = 0                 # polarization condition (only used with input_field_rootname)
    source_x_offset = 0         # source offset in lambda0_m/D radians (tilt applied at primary)
    source_y_offset = 0                 
    use_hlc_dm_patterns = 0     # use Dwight's HLC default DM wavefront patterns? 1 or 0
    use_dm1 = 0                 # use DM1? 1 or 0
    use_dm2 = 0                 # use DM2? 1 or 0
    dm_sampling_m = 0.9906e-3   # actuator spacing in meters
    dm1_xc_act = 23.5           # for 48x48 DM, wavefront centered at actuator intersections: (0,0) = 1st actuator center
    dm1_yc_act = 23.5              
    dm1_xtilt_deg = 0           # tilt around X axis (deg)
    dm1_ytilt_deg = 5.7         # effective DM tilt in deg including 9.65 deg actual tilt and pupil ellipticity
    dm1_ztilt_deg = 0           # rotation of DM about optical axis (deg)
    dm2_xc_act = 23.5           # for 48x48 DM, wavefront centered at actuator intersections: (0,0) = 1st actuator center
    dm2_yc_act = 23.5              
    dm2_xtilt_deg = 0           # tilt around X axis (deg)
    dm2_ytilt_deg = 5.7         # effective DM tilt in deg including 9.65 deg actual tilt and pupil ellipticity
    dm2_ztilt_deg = 0           # rotation of DM about optical axis (deg)
    fpm_axis = 'p'              # HLC FPM axis: '', 's', 'p'
    final_sampling_lam0 = 0     # final sampling in lambda0/D
    output_dim = output_dim0    # dimension of output in pixels (overrides output_dim0)

    if 'PASSVALUE' in locals():
        if 'cor_type' in PASSVALUE: cor_type = PASSVALUE['cor_type']
        if 'fpm_axis' in PASSVALUE: fpm_axis = PASSVALUE['fpm_axis']

    is_hlc = False
    is_spc = False

    if cor_type == 'hlc':
        is_hlc = True
        file_directory = data_dir + '/hlc_20190210/'         # must have trailing "/"
        prefix = file_directory + 'run461_' 
        pupil_diam_pix = 309.0
        pupil_file = prefix + 'pupil_rotated.fits'
        lyot_stop_file = prefix + 'lyot.fits'
        lambda0_m = 0.575e-6
        lam_occ = [     5.4625e-07, 5.49444444444e-07, 5.52638888889e-07, 5.534375e-07, 5.55833333333e-07, 5.59027777778e-07, 5.60625e-07, 
                        5.62222222222e-07, 5.65416666667e-07, 5.678125e-07, 5.68611111111e-07, 5.71805555556e-07, 5.75e-07, 5.78194444444e-07, 
                        5.81388888889e-07, 5.821875e-07, 5.84583333333e-07, 5.87777777778e-07, 5.89375e-07, 5.90972222222e-07, 5.94166666667e-07, 
                        5.965625e-07, 5.97361111111e-07, 6.00555555556e-07, 6.0375e-07 ]
        lam_occs = [    '5.4625e-07', '5.49444444444e-07', '5.52638888889e-07', '5.534375e-07', '5.55833333333e-07', '5.59027777778e-07', 
                        '5.60625e-07', '5.62222222222e-07', '5.65416666667e-07', '5.678125e-07', '5.68611111111e-07', '5.71805555556e-07', 
                        '5.75e-07', '5.78194444444e-07', '5.81388888889e-07', '5.821875e-07', '5.84583333333e-07', '5.87777777778e-07', 
                        '5.89375e-07', '5.90972222222e-07', '5.94166666667e-07', '5.965625e-07', '5.97361111111e-07', '6.00555555556e-07', '6.0375e-07' ]
        lam_occs = [ prefix + 'occ_lam' + s + 'theta6.69pol' + fpm_axis + '_' for s in lam_occs ]
        # find nearest matching FPM wavelength
        wlam = (np.abs(lambda_m-np.array(lam_occ))).argmin()
        occulter_file_r = lam_occs[wlam] + 'real_rotated.fits'
        occulter_file_i = lam_occs[wlam] + 'imag_rotated.fits'
        n_small = 1024                  # gridsize in non-critical areas
        n_big = 2048                    # gridsize to/from FPM
    elif cor_type == 'hlc_erkin':
        is_hlc = True
        file_directory = data_dir + '/hlc_20190206_v3/'         # must have trailing "/"
        prefix = file_directory + 'dsn17d_run2_pup310_fpm2048_'
        pupil_diam_pix = 310.0
        pupil_file = prefix + 'pupil.fits'
        lyot_stop_file = prefix + 'lyot.fits'
        lambda0_m = 0.575e-6
        lam_occ = [     5.4625e-07, 5.4944e-07, 5.5264e-07, 5.5583e-07, 5.5903e-07, 5.6222e-07, 5.6542e-07,
                        5.6861e-07, 5.7181e-07, 5.75e-07, 5.7819e-07, 5.8139e-07, 5.8458e-07, 5.8778e-07,
                        5.9097e-07, 5.9417e-07, 5.9736e-07, 6.0056e-07, 6.0375e-07 ]
        lam_occs =  [   '5.4625e-07', '5.4944e-07', '5.5264e-07', '5.5583e-07', '5.5903e-07', '5.6222e-07', '5.6542e-07', 
                        '5.6861e-07', '5.7181e-07', '5.75e-07', '5.7819e-07', '5.8139e-07', '5.8458e-07', '5.8778e-07',
                        '5.9097e-07', '5.9417e-07', '5.9736e-07', '6.0056e-07', '6.0375e-07' ]
        fpm_axis = 's'
        lam_occs = [ prefix + 'occ_lam' + s + 'theta6.69pol' + fpm_axis + '_' for s in lam_occs ]
        # find nearest matching FPM wavelength
        wlam = (np.abs(lambda_m-np.array(lam_occ))).argmin()
        occulter_file_r = lam_occs[wlam] + 'real.fits'
        occulter_file_i = lam_occs[wlam] + 'imag.fits'
        n_small = 1024                  # gridsize in non-critical areas
        n_big = 2048                    # gridsize to/from FPM
    elif cor_type == 'spc-ifs_short' or cor_type == 'spc-ifs_long' or cor_type == 'spc-spec_short' or cor_type == 'spc-spec_long':
        is_spc = True
        file_dir = data_dir + '/spc_20190130/' # must have trailing "/"
        pupil_diam_pix = 1000.0
        pupil_file = file_dir + 'pupil_SPC-20190130_rotated.fits'
        pupil_mask_file = file_dir + 'SPM_SPC-20190130_rotated.fits'
        fpm_file = file_dir + 'fpm_0.05lamdivD.fits' 
        fpm_sampling = 0.05    # sampling in lambda0/D of FPM mask 
        if cor_type == 'spc-ifs_short' or cor_type == 'spc-spec_short':
            fpm_sampling_lambda_m = 0.66e-6
            lambda0_m = 0.66e-6
        else:
            fpm_sampling_lambda_m = 0.73e-6
            lambda0_m = 0.73e-6     # FPM scaled for this central wavelength
        lyot_stop_file = file_dir + 'lyotstop_0.5mag.fits'
        n_small = 2048              # gridsize in non-critical areas
        n_big = 1400                # gridsize to FPM (propagation to/from FPM handled by MFT)
    elif cor_type == 'spc-wide':
        is_spc = True
        file_dir = data_dir + '/spc_20181220/' # must have trailing "/"
        pupil_diam_pix = 1000.0
        pupil_file = file_dir + 'pupil_SPC-20181220_1k_rotated.fits'
        pupil_mask_file = file_dir + 'SPM_SPC-20181220_1000_rounded9_gray_rotated.fits'
        fpm_file = file_dir + 'fpm_0.05lamdivD.fits'
        fpm_sampling = 0.05    # sampling in lambda0/D of FPM mask 
        lyot_stop_file = file_dir + 'LS_half_symm_CGI180718_Str3.20pct_38D91_N500_pixel.fits'
        fpm_sampling_lambda_m = 0.825e-6
        lambda0_m = 0.825e-6        # FPM scaled for this central wavelength
        n_small = 2048              # gridsize in non-critical areas
        n_big = 1400                # gridsize to FPM (propagation to/from FPM handled by MFT)
    elif cor_type == 'none':
        file_directory = data_dir + '/hlc_20190210/'         # must have trailing "/"
        prefix = file_directory + 'run461_' 
        pupil_diam_pix = 309.0
        pupil_file = prefix + 'pupil_rotated.fits'
        use_fpm = 0
        use_lyot_stop = 0
        n_small = 1024
        n_big = 1024
    else:
        raise Exception( 'wfirst_phaseb_compact: Unsuported cor_type: '+cor_type )

    if 'PASSVALUE' in locals():
        if 'lam0' in PASSVALUE: lamba0_m = PASSVALUE['lam0'] * 1.0e-6
        if 'lambda0_m' in PASSVALUE: lambda0_m = PASSVALUE['lambda0_m']
        if 'input_field_rootname' in PASSVALUE: input_field_rootname = PASSVALUE['input_field_rootname'] 
        if 'polaxis' in PASSVALUE: polaxis = PASSVALUE['polaxis']
        if 'source_x_offset' in PASSVALUE: source_x_offset = PASSVALUE['source_x_offset']
        if 'source_y_offset' in PASSVALUE: source_y_offset = PASSVALUE['source_y_offset']
        if 'use_hlc_dm_patterns' in PASSVALUE: use_hlc_dm_patterns = PASSVALUE['use_hlc_dm_patterns']
        if 'use_dm1' in PASSVALUE: use_dm1 = PASSVALUE['use_dm1'] 
        if 'dm1_m' in PASSVALUE: dm1_m = PASSVALUE['dm1_m']
        if 'dm1_xc_act' in PASSVALUE: dm1_xc_act = PASSVALUE['dm1_xc_act']
        if 'dm1_yc_act' in PASSVALUE: dm1_yc_act = PASSVALUE['dm1_yc_act']
        if 'dm1_xtilt_deg' in PASSVALUE: dm1_xtilt_deg = PASSVALUE['dm1_xtilt_deg']
        if 'dm1_ytilt_deg' in PASSVALUE: dm1_ytilt_deg = PASSVALUE['dm1_ytilt_deg']
        if 'dm1_ztilt_deg' in PASSVALUE: dm1_ztilt_deg = PASSVALUE['dm1_ztilt_deg']
        if 'use_dm2' in PASSVALUE: use_dm2 = PASSVALUE['use_dm2']
        if 'dm2_m' in PASSVALUE: dm2_m = PASSVALUE['dm2_m']
        if 'dm2_xc_act' in PASSVALUE: dm2_xc_act = PASSVALUE['dm2_xc_act']
        if 'dm2_yc_act' in PASSVALUE: dm2_yc_act = PASSVALUE['dm2_yc_act']
        if 'dm2_xtilt_deg' in PASSVALUE: dm2_xtilt_deg = PASSVALUE['dm2_xtilt_deg']
        if 'dm2_ytilt_deg' in PASSVALUE: dm2_ytilt_deg = PASSVALUE['dm2_ytilt_deg']
        if 'dm2_ztilt_deg' in PASSVALUE: dm2_ztilt_deg = PASSVALUE['dm2_ztilt_deg']
        if 'final_sampling_lam0' in PASSVALUE: final_sampling_lam0 = PASSVALUE['final_sampling_lam0']
        if 'output_dim' in PASSVALUE: output_dim = PASSVALUE['output_dim']

    if polaxis != 0 and input_field_rootname == '':
        raise Exception( 'wfirst_phaseb_compact: polaxis can only be defined when input_field_rootname is given' )

    diam_at_dm1 = 0.0463
    d_dm1_dm2 = 1.0

    n = n_small        # start off with less padding
 
    wavefront = proper.prop_begin( diam_at_dm1, lambda_m, n, float(pupil_diam_pix)/n )
    if input_field_rootname == '':
        pupil = proper.prop_fits_read( pupil_file )
        proper.prop_multiply( wavefront, trim(pupil,n) )
        pupil = 0
    else:
        lams = format( lambda_m*1e6, "6.4f" )
        pols = format( int(round(polaxis)) )
        rval = proper.prop_fits_read( input_field_rootname+'_'+lams+'um_'+pols+'_real.fits' )
        ival = proper.prop_fits_read( input_field_rootname+'_'+lams+'um_'+pols+'_imag.fits' )
        proper.prop_multiply( wavefront, trim(rval + 1j * ival, n) )
        rval = 0
        ival = 0
    proper.prop_define_entrance( wavefront )
    if source_x_offset != 0 or source_y_offset != 0:
        # compute tilted wavefront to offset source by xoffset,yoffset lambda0_m/D
        xtilt_lam = -source_x_offset * lambda0_m / lambda_m
        ytilt_lam = -source_y_offset * lambda0_m / lambda_m
        x = np.tile( (np.arange(n)-n//2)/(pupil_diam_pix/2.0), (n,1) )
        y = np.transpose(x)
        proper.prop_multiply( wavefront, np.exp(complex(0,1) * np.pi * (xtilt_lam * x + ytilt_lam * y)) )
        x = 0
        y = 0
    
    if use_dm1 != 0: prop_dm( wavefront, dm1_m, dm1_xc_act, dm1_yc_act, dm_sampling_m, XTILT=dm1_xtilt_deg, YTILT=dm1_ytilt_deg, ZTILT=dm1_ztilt_deg )
    if is_hlc == True and use_hlc_dm_patterns == 1:
        dm1wfe = proper.prop_fits_read( prefix+'dm1wfe.fits' )
        proper.prop_add_phase( wavefront, trim(dm1wfe, n) )
        dm1wfe = 0

    proper.prop_propagate( wavefront, d_dm1_dm2, 'DM2' )
    if use_dm2 == 1: prop_dm( wavefront, dm2_m, dm2_xc_act, dm2_yc_act, dm_sampling_m, XTILT=dm2_xtilt_deg, YTILT=dm2_ytilt_deg, ZTILT=dm2_ztilt_deg )
    if is_hlc == True:
        if use_hlc_dm_patterns == 1:
            dm2wfe = proper.prop_fits_read( prefix+'dm2wfe.fits' )
            proper.prop_add_phase( wavefront, trim(dm2wfe, n) )
            dm2wfe = 0
        dm2mask = proper.prop_fits_read( prefix+'dm2mask.fits' )
        proper.prop_multiply( wavefront, trim(dm2mask, n) )
        dm2mask = 0

    proper.prop_propagate( wavefront, -d_dm1_dm2, 'back to DM1' )

    (wavefront, sampling_m) = proper.prop_end( wavefront, NOABS=True )

    # apply shape pupil mask

    if is_spc == True:
        pupil_mask = proper.prop_fits_read( pupil_mask_file )
        wavefront *= trim(pupil_mask,n)
        pupil_mask = 0

    # propagate to FPM and apply FPM

    if is_hlc == True:
        n = n_big
        wavefront = trim(wavefront,n)
        wavefront = ffts(wavefront,-1)  # to focus
        occ_r = proper.prop_fits_read( occulter_file_r )
        occ_i = proper.prop_fits_read( occulter_file_i )
        occ = np.array( occ_r + 1j * occ_i, dtype=np.complex128 )
        wavefront *= trim(occ,n)
        occ_r = 0
        occ_i = 0
        occ = 0
        wavefront = ffts(wavefront,+1)  # to lyot stop
    elif is_spc == True:
        n = n_big
        wavefront = trim(wavefront,n)
        fpm = proper.prop_fits_read( fpm_file )
        nfpm = fpm.shape[1]
        fpm_sampling_lam = fpm_sampling * fpm_sampling_lambda_m / lambda_m
        wavefront = mft2(wavefront, fpm_sampling_lam, pupil_diam_pix, nfpm, -1)   # MFT to highly-sampled focal plane
        wavefront *= fpm
        fpm = 0
        pupil_diam_pix = pupil_diam_pix / 2.0   # Shrink pupil by 1/2
        wavefront = mft2(wavefront, fpm_sampling_lam, pupil_diam_pix, int(pupil_diam_pix), +1)  # MFT to Lyot stop with 1/2 magnification

    n = n_small
    wavefront = trim(wavefront,n)
    lyot = proper.prop_fits_read( lyot_stop_file )
    wavefront *= trim(lyot,n)
    lyot = 0

    wavefront *= n
    wavefront = ffts(wavefront,-1)    # to focus

    # rotate to convention used by full prescription

    wavefront[:,:] = np.rot90( wavefront, 2 )
    wavefront[:,:] = np.roll( wavefront, 1, axis=0 ) 
    wavefront[:,:] = np.roll( wavefront, 1, axis=1 )
 
    if final_sampling_lam0 != 0:
        mag = (float(pupil_diam_pix)/n) / final_sampling_lam0 * (lambda_m/lambda0_m)
        wavefront = proper.prop_magnify( wavefront, mag, output_dim, AMP_CONSERVE=True )
    else:
        wavefront = trim(wavefront, output_dim)

    sampling_m = 0.0
    return wavefront, sampling_m

