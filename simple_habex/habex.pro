; Copyright 2020, by the California Institute of Technology. ALL RIGHTS
; RESERVED. United States Government Sponsorship acknowledged. Any
; commercial use must be negotiated with the Office of Technology Transfer
; at the California Institute of Technology.
;++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
; Name:
;	habex.pro
;
; Purpose:
;	Representation of the Habex telescope and coronagraph. To be called using 
;	the PROPER library procedure "prop_run".
;
; Inputs:
;	lambda_m
;	    The wavelength of propagation in meters (note that the wavelength is provided
;	    to prop_run in microns and is converted to meters in there).
;	gridsize
;	    Size of the computational grid (gridsize by gridsize elements).  Must be
;	    a power of 2. 
;
; Outputs:
;	wavefront
;	    Variable in which the computed E-field at the final image plane is returned.
;	    The field is sampled by "final_sampling_lam0" lambda_m/D over "nout" by "nout"
;	    pixels.
;	sampling_m
;	    The sampling at the final image plane in meters per pixel
;
; Optional keywords or switches:
;	PASSVALUE=optval
;	    (Optional) Points to a value (which could be a constant or a variable) 
;	    that is passed to the prescription for use as the prescription desires.
;
; Revision history:
;	Written by John Krist (Jet Propulsion Laboratory, California Inst. Technology), January 2020 
;----------------------------------------------------------------------------------  

pro habex, wavefront, lambda_m, gridsize, sampling_m, PASSVALUE=optval

map_dir = '/home/krist/habex/maps/'	;-- directory containing optical surface error maps

nact = 64                       ;-- number of actuators across DM
nact_across_pupil = 62		;-- number of actuators across pupil
dm_xc = 31.5d                   ;-- wavefront centered at corner of DM actuator (0,0 is center of 1st actuator)
dm_yc = 31.5d
dm_sampling = 0.4d-3            ;-- DM actuator spacing (BMC) 
pupil_diam_pix = nact_across_pupil * 7.0 	;-- define sampling of pupil based on having 7 pixels across each DM actuator
pupil_ratio = pupil_diam_pix / double(gridsize)

;-- default settings (override with optval)

lambda0_um = 0.5d	;-- default reference wavelength (center of bandpass) for star offsets & field stop size
use_errors = 1		;-- 1 = use optical surface errors, 0 = none
zindex = 0		;-- vector of Zernike indices (Noll ordered)
zval = 0		;-- vector of Zernike coefficients (unobscured RMS wavefront in meters)
xoffset = 0		;-- star X offset in lambda0/D units (must then provide lambda0_um)
yoffset = 0		;-- star Y offset in lambda0/D units
use_dm1 = 0		;-- use DM1 (if non-zero, must then provide pokes (meters) in "dm1" array)
use_dm2 = 0		;-- use DM2 (if non-zero, must then provide pokes (meters) in "dm2" array)
use_fpm = 1		;-- use focal plane mask (0 = no FPM)
use_lyot_stop = 1	;-- use Lyot stop (0 = no stop)
use_field_stop = 1	;-- use field stop (0 = no stop)
field_stop_radius = 25.0   ;-- field stop radius in lam0/D
final_sampling_lam0 = 0.2d ;-- sampling at final image plane in lam0/D
nout = 300		;-- output field size (nout x nout pixels)

;-- override defaults using values passed using optval structure

if ( n_elements(optval) ne 0 ) then begin
	if ( tag_exists('lam0',optval) ) then lambda0_um = optval.lam0 
	if ( tag_exists('lambda0_um',optval) ) then lambda0_um = optval.lambda0_um 
	if ( tag_exists('use_errors',optval) ) then use_errors = optval.use_errors
	if ( tag_exists('zindex',optval) ) then zindex = optval.zindex
	if ( tag_exists('zval',optval) ) then zval = optval.zval
	if ( tag_exists('xoffset',optval) ) then xoffset = optval.xoffset
	if ( tag_exists('yoffset',optval) ) then yoffset = optval.yoffset
	if ( tag_exists('use_dm1',optval) ) then use_dm1 = optval.use_dm1
	if ( tag_exists('dm1',optval) ) then dm1 = optval.dm1
	if ( tag_exists('use_dm2',optval) ) then use_dm2 = optval.use_dm2
	if ( tag_exists('dm2',optval) ) then dm2 = optval.dm2
	if ( tag_exists('use_fpm',optval) ) then use_fpm = optval.use_fpm
	if ( tag_exists('use_lyot_stop',optval) ) then use_lyot_stop = optval.use_lyot_stop
	if ( tag_exists('use_field_stop',optval) ) then use_field_stop = optval.use_field_stop
	if ( tag_exists('field_stop_radius',optval) ) then field_stop_radius = optval.field_stop_radius
	if ( tag_exists('final_sampling_lam0',optval) ) then final_sampling_lam0 = optval.final_sampling_lam0
	if ( tag_exists('nout',optval) ) then nout = optval.nout
endif

lambda0_m = lambda0_um * 1.0d-6

;-- define optical prescription (distances, focal lengths)

diam = 4.00d
  r_pri = 19.8d
  h_pri = 2.5d
  z_pri = h_pri^2 / (2*r_pri)
  fl_pri = sqrt(h_pri^2 + (r_pri/2-z_pri)^2)	;-- effective focal length of primary as a pure parabola
d_pri_sec = 9.172532289071727d
  d_focus_sec = fl_pri - d_pri_sec
  d_sec_focus = 7.979857207574376844d
  fl_sec = 1 / (1/d_sec_focus - 1/d_focus_sec)  
d_sec_m3 = 9.076690863872008d
  fl_m3 = d_sec_m3 - d_sec_focus
d_m3_fold = 0.654597300210990d
d_fold_fsm = 0.577743120280288d
d_fsm_dichroic = 0.1950d
d_dichroic_m4 = 0.450d
  fl_m4 = 0.5075d
d_m4_m5 = 0.762954002022743d
  fl_m5 = d_m4_m5 - fl_m4	
d_m5_dm1 = 0.220615776458241d
d_dm1_dm2 = 0.32d
d_dm2_qwp = 0.32d + 0.157485214529470d
  fl_m6 = 1.029143136045496931d
d_qwp_m6 = fl_m6 - (d_dm1_dm2 + d_dm2_qwp)
d_m6_fpm = fl_m6
d_fpm_m7 = 0.255580492381039d
  fl_m7 = d_fpm_m7
d_m7_lyotstop = fl_m7 
d_lyotstop_m8 = 0.2536d	
  fl_m8 = d_lyotstop_m8
d_m8_fieldstop = fl_m8
d_fieldstop_m9 = d_m8_fieldstop
  fl_m9 = d_fieldstop_m9
d_m9_filter = 0.296399999724129d
d_filter_m10 = 0.462615469378302d
  fl_m10 = 0.503971038519431261d
d_m10_ccd = fl_m10	


prop_begin, wavefront, diam, lambda_m, gridsize, pupil_diam_pix/gridsize
   prop_circular_aperture, wavefront, diam/2
   if ( zindex[0] ne 0 ) then prop_zernikes, wavefront, zindex, zval	;-- optionally add Zernikes
   if ( xoffset ne 0 or yoffset ne 0 ) then begin
	;-- star X,Y offset in lam0/D
	i = dcomplex(0,1)
	xoffset_lam = xoffset * lambda0_m / lambda_m
	yoffset_lam = yoffset * lambda0_m / lambda_m
	u = (dindgen(gridsize)-gridsize/2) / (pupil_diam_pix/2)
	xtilt = exp( i * !dpi * u * xoffset_lam )
	ytilt = exp( i * !dpi * u * yoffset_lam )
	prop_multiply, wavefront, xtilt # ytilt
	xtilt = 0
	ytilt = 0
   endif
   if ( use_errors ) then prop_errormap, wavefront, map_dir+'habex_cycle1_PRIMARY_phase_error.fits', /WAVEFRONT
   prop_lens, wavefront, fl_pri
   prop_define_entrance, wavefront

prop_propagate, wavefront, d_pri_sec, 'secondary'
   if ( use_errors ) then prop_errormap, wavefront, map_dir+'habex_cycle1_SECONDARY_phase_error.fits', /WAVEFRONT
   prop_lens, wavefront, fl_sec

prop_propagate, wavefront, d_sec_m3, 'M3'
   if ( use_errors ) then prop_errormap, wavefront, map_dir+'habex_cycle1_M3_phase_error.fits', /WAVEFRONT
   prop_lens, wavefront, fl_m3

prop_propagate, wavefront, d_m3_fold, 'fold'
   if ( use_errors ) then prop_errormap, wavefront, map_dir+'habex_cycle1_FOLD1_phase_error.fits', /WAVEFRONT

prop_propagate, wavefront, d_fold_fsm, 'FSM'	;-- pupil at fast steering mirror (interface with telescope)
   if ( use_errors ) then prop_errormap, wavefront, map_dir+'habex_cycle1_FSM_phase_error.fits', /WAVEFRONT

prop_propagate, wavefront, d_fsm_dichroic, 'dichroic'
   if ( use_errors ) then prop_errormap, wavefront, map_dir+'habex_cycle1_DICHROIC_phase_error.fits', /WAVEFRONT

prop_propagate, wavefront, d_dichroic_m4, 'M4'
   if ( use_errors ) then prop_errormap, wavefront, map_dir+'habex_cycle1_M4_phase_error.fits', /WAVEFRONT
   prop_lens, wavefront, fl_m4

prop_propagate, wavefront, d_m4_m5, 'M5'
   if ( use_errors ) then prop_errormap, wavefront, map_dir+'habex_cycle1_M5_phase_error.fits', /WAVEFRONT
   prop_lens, wavefront, fl_m5

prop_propagate, wavefront, d_m5_dm1, 'DM1'
   if ( use_dm1 ) then prop_dm, wavefront, dm1, dm_xc, dm_yc, dm_sampling
   if ( use_errors ) then prop_errormap, wavefront, map_dir+'habex_cycle1_DM1_phase_error.fits', /WAVEFRONT

prop_propagate, wavefront, d_dm1_dm2, 'DM2'
   if ( use_dm2 ) then prop_dm, wavefront, dm2, dm_xc, dm_yc, dm_sampling
   if ( use_errors ) then prop_errormap, wavefront, map_dir+'habex_cycle1_DM2_phase_error.fits', /WAVEFRONT

prop_propagate, wavefront, d_dm2_qwp, 'QWP'	;-- quarter-wave plate
   if ( use_errors ) then prop_errormap, wavefront, map_dir+'habex_cycle1_QWP1_phase_error.fits', /WAVEFRONT

prop_propagate, wavefront, d_qwp_m6, 'M6'
   if ( use_errors ) then prop_errormap, wavefront, map_dir+'habex_cycle1_M6_phase_error.fits', /WAVEFRONT
   prop_lens, wavefront, fl_m6

prop_propagate, wavefront, d_m6_fpm
   if ( use_fpm ) then prop_8th_order_mask, wavefront, 4.0, /CIRCULAR	

prop_propagate, wavefront, d_fpm_m7, 'M7'
   if ( use_errors ) then prop_errormap, wavefront, map_dir+'habex_cycle1_M7_phase_error.fits', /WAVEFRONT
   prop_lens, wavefront, fl_m7

prop_propagate, wavefront, d_m7_lyotstop, 'Lyot stop'
   if ( use_errors ) then prop_errormap, wavefront, map_dir+'habex_cycle1_QWP2_phase_error.fits', /WAVEFRONT
   if ( use_lyot_stop ) then prop_circular_aperture, wavefront, 0.54, /NORM

prop_propagate, wavefront, d_lyotstop_m8, 'M8'
   if ( use_errors ) then prop_errormap, wavefront, map_dir+'habex_cycle1_M8_phase_error.fits', /WAVEFRONT
   prop_lens, wavefront, fl_m8

prop_propagate, wavefront, prop_get_distancetofocus(wavefront), 'field stop'
   if ( use_field_stop ) then begin
   	r_stop = field_stop_radius * lambda0_m / lambda_m
   	prop_circular_aperture, wavefront, r_stop/pupil_ratio*prop_get_sampling(wavefront)
   endif

prop_propagate, wavefront, d_fieldstop_m9, 'M9'
   if ( use_errors ) then prop_errormap, wavefront, map_dir+'habex_cycle1_M9_phase_error.fits', /WAVEFRONT
   prop_lens, wavefront, fl_m9

prop_propagate, wavefront, d_m9_filter, 'filter'
   if ( use_errors ) then prop_errormap, wavefront, map_dir+'habex_cycle1_FILTER_phase_error.fits', /WAVEFRONT

prop_propagate, wavefront, d_filter_m10, 'M10'
   if ( use_errors ) then prop_errormap, wavefront, map_dir+'habex_cycle1_M10_phase_error.fits', /WAVEFRONT
   prop_lens, wavefront, fl_m10

prop_propagate, wavefront, prop_get_distancetofocus(wavefront), 'CCD'

prop_end, wavefront, sampling_m, /NOABS

;-- rescale to "final_sampling_lam0" lam0/D per pixel

mag = (pupil_ratio / final_sampling_lam0) * (lambda_m / lambda0_m)
wavefront = prop_magnify( wavefront, mag, nout, /AMP_CONSERVE )

return
end
