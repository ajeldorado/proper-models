;-- Version 1.0, 3 January 2019, JEK
;-- Version 1.0a, 9 January 2019, JEK 
;--    Changes: 1) now append polaxis to output_field_rootname 
;--             2) can output field at FPM exit pupil (via output_field_rootname) without 
;--                having to stop there
;-- Version 1.0b, 16 January 2019, JEK 
;--    Changes: Multiplied wfirst_phaseb_GROUND_TO_ORBIT_phase_error_V1.0.fits by 4.2x to get the
;--		  the RMS WFE at the CGI entrance to be 76.4 nm RMS.  The new maps is in 
;--		  result to wfirst_phaseb_GROUND_TO_ORBIT_4.2X_phase_error_V1.0.fits 
;-- Version 1.0c, 13 February 2019, JEK 
;--    Changes: Changed HLC design to HLC-20190210; added 'fpm_axis' parameter to specify
;--		  axis of HLC FPM ('s' or 'p') due to FPM tilt; changed 'spc' coronagraph
;--		  option to 'cor_type' to 'spc-ifs' and changed default SPC to SPC-20190130 and
;--		  updated associated parameters; added 'spc-wide' option to 'cor_type' to use
;--		  SPC-20181220. Changed DM tilt to rotation around Y axis.

pro wfirst_phaseb, wavefront, lambda_m, output_dim, sampling_m, PASSVALUE=optval


;-- "output_dim" is used to specify the output dimension in pixels at the final image plane.  
;-- The computational grid sizes are hardcoded for each coronagraph.
;-- Based on Zemax prescription "WFIRST_CGI_DI_LOWFS_Sep24_2018.zmx" by Hong Tang.
;--
;-- by John Krist

map_dir = './maps/'		;-- directory for surface maps
polfile = './new_toma'		;-- polarization aberration table rootname

cor_type = 'hlc'		;-- 'hlc', 'spc-ifs', 'spc-wide', or 'none' (none = clear aperture, no coronagraph)
source_x_offset_mas = 0		;-- source offset in milliarcsec (tilt applied at primary)
source_y_offset_mas = 0
source_x_offset = 0		;-- source offset in lambda0_m/D radians (tilt applied at primary)
source_y_offset = 0
polaxis = 0			;-- polarization axis aberrations: 
				;--     -2 = -45d in, Y out 
				;--     -1 = -45d in, X out 
				;--      1 = +45d in, X out 
				;--      2 = +45d in, Y out 
				;--      5 = mean of modes -1 & +1 (X channel polarizer)
				;--      6 = mean of modes -2 & +2 (Y channel polarizer)
				;--     10 = mean of all modes (no polarization filtering)
use_errors = 1			;-- use optical surface phase errors? 1 or 0
zindex = 0			;-- array of Zernike polynomial indices
zval_m = 0			;-- array of Zernike coefficients (meters RMS WFE)
use_aperture = 0		;-- use apertures on all optics? 1 or 0
end_at_fsm = 0			;-- end propagation after propagating to FSM (no FSM errors)
cgi_x_shift_pupdiam = 0		;-- X,Y shift of wavefront at FSM (bulk displacement of CGI); normalized relative to pupil diameter
cgi_y_shift_pupdiam = 0
fsm_x_offset_mas = 0		;-- offset source using FSM (milliarcsec)
fsm_y_offset_mas = 0
fsm_x_offset = 0		;-- offset source using FSM (lambda0/D radians)
fsm_y_offset = 0
focm_z_shift_m = 0		;-- offset (meters) of focus correction mirror (+ increases path length)
use_hlc_dm_patterns = 1		;-- use Dwight-generated HLC default DM wavefront patterns? 1 or 0
use_dm1 = 0			;-- use DM1? 1 or 0
use_dm2 = 0			;-- use DM2? 1 or 0
dm_sampling_m = 1.0d-3         	;-- actuator spacing in meters; default is 1 mm
dm1_xc_act = 23.5d		;-- for 48x48 DM, wavefront centered at actuator intersections: (0,0) = 1st actuator center
dm1_yc_act = 23.5d
dm1_xtilt_deg = 0   		;-- tilt around X axis
dm1_ytilt_deg = 5.7d		;-- effective DM tilt in deg including 9.65 deg actual tilt and pupil ellipticity
dm1_ztilt_deg = 0
dm2_xc_act = 23.5d		
dm2_yc_act = 23.5d
dm2_xtilt_deg = 0   
dm2_ytilt_deg = 5.7d
dm2_ztilt_deg = 0
mask_x_shift_pupdiam = 0	;-- X,Y shift of shaped pupil mask; normalized relative to pupil diameter
mask_y_shift_pupdiam = 0
use_fpm = 1			;-- use occulter? 1 or 0
fpm_axis = 'p'                  ;-- HLC FPM axis: '', 's', 'p'
fpm_x_offset = 0		;-- FPM x,y offset in lambda0/D
fpm_y_offset = 0
fpm_z_shift_m = 0		;-- FPM offset in meters along optical axis (+ = away from prior optics)
pinhole_diam_m = 0		;-- pinhole diameter in meters at FPM
end_at_fpm_exit_pupil = 0	;-- return field at FPM exit pupil?
output_field_rootname = ''	;-- rootname of FPM exit pupil field file (must set end_at_fpm_exit_pupil=1)
use_lyot_stop = 1		;-- use Lyot stop? 1 or 0
lyot_x_shift_pupdiam = 0	;-- X,Y shift of Lyot stop mask; normalized relative to pupil diameter
lyot_y_shift_pupdiam = 0
use_field_stop = 1		;-- use field stop (HLC)? 1 or 0
field_stop_radius_lam0 = 0	;-- field stop radius in lambda0/D (HLC or SPC-wide mask only)
use_pupil_lens = 0		;-- use pupil imaging lens?
use_defocus_lens = 0		;-- use defocusing lens? Options are 1, 2, 3, 4, corresponding to +18.0, +9.0, -4.0, -8.0 waves P-V @ 550 nm 
defocus = 0			;-- instead of specific lens, defocus in waves P-V @ 550 nm (-8.7 to 42.0 waves)
final_sampling_lam0 = 0		;-- final sampling in lambda0/D

if ( n_elements(optval) ne 0 ) then begin
	if ( tag_exists('cor_type',optval) ) then cor_type = optval.cor_type
	if ( tag_exists('use_fpm',optval) ) then use_fpm = optval.use_fpm
	if ( tag_exists('fpm_axis',optval) ) then fpm_axis = optval.fpm_axis
endif

if ( cor_type eq 'hlc' ) then begin
        file_directory = './hlc_20190210/'         ;-- must have trailing "/"
        prefix = file_directory + 'run461_nro_' 
        pupil_diam_pix = 309.0
        pupil_file = prefix + 'pupil_rotated.fits'
        lyot_stop_file = prefix + 'lyot.fits'
        lambda0_m = 0.575d-6
	lam_occ =                       [  5.4625d-07,  5.534375d-07,  5.60625d-07,  5.678125d-07,  5.75d-07,  5.821875e-07,  5.89375d-07,  5.965625d-07,  6.0375d-07 ]
	lam_occs = prefix + 'occ_lam' + [ '5.4625e-07','5.534375e-07','5.60625e-07','5.678125e-07','5.75e-07','5.821875e-07','5.89375e-07','5.965625e-07','6.0375e-07' ] 
	lam_occs = lam_occs + 'theta6.69pol' + fpm_axis + '_'
        ;-- find nearest matching FPM wavelength
        diff = abs(lambda_m-lam_occ)
        wlam = (where(diff eq min(diff)))[0]
        occulter_file_r = lam_occs[wlam] + 'real.fits'
        occulter_file_i = lam_occs[wlam] + 'imag.fits'
	n_default = 1024	;-- gridsize in non-critical areas
	if ( use_fpm ne 0 ) then n_to_fpm = 2048 else n_to_fpm = 1024
	n_from_lyotstop = 1024
	;-- these commented-out DM settings match the testbed rather than flight
        ;dm_sampling_m = 0.9906d-3
        ;dm1_xtilt_deg = 0.0
        ;dm1_ytilt_deg = 10.0
        ;dm1_ztilt_deg = -0.88
        ;dm1_xtilt_deg = 0.0
        ;dm1_ytilt_deg = 8.0
        ;dm1_ztilt_deg =  0.45
	field_stop_radius_lam0 = 9.0
endif else if ( cor_type eq 'spc-ifs' ) then begin
        file_dir = './spc_20190130/'       ;-- must have trailing "/"
        pupil_diam_pix = 1000.0
        pupil_file = file_dir + 'pupil_SPC-20190130.fits'
        pupil_mask_file = file_dir + 'SPM_SPC-20190130_rotated.fits'
        fpm_file = file_dir + 'fpm_0.05lamdivD.fits'
        fpm_sampling_lam0 = 0.05d       ;-- sampling in lambda0/D of FPM mask 
        lyot_stop_file = file_dir + 'LS_SPC-20190130.fits'
        lambda0_m = 0.73d-6     ;-- FPM scaled for this central wavelength
	n_default = 2048	;-- gridsize in non-critical areas
	n_to_fpm = 2048	;-- gridsize to/from FPM
	n_mft = 1400
	n_from_lyotstop = 4096 
endif else if ( cor_type eq 'spc-wide' ) then begin
        file_dir = './spc_20181220/'       ;-- must have trailing "/"
        pupil_diam_pix = 1000.0
        pupil_file = file_dir + 'pupil_SPC-20181220_1k.fits'
        pupil_mask_file = file_dir + 'SPM_SPC-20181220_1000_rounded9_gray_rotated.fits'
        fpm_file = file_dir + 'fpm_0.05lamdivD.fits'
        fpm_sampling_lam0 = 0.05d       ;-- sampling in lambda0/D of FPM mask 
        lyot_stop_file = file_dir + 'LS_SPC-20181220_1k_rotated.fits'
        lambda0_m = 0.825d-6     ;-- FPM scaled for this central wavelength
	n_default = 2048	;-- gridsize in non-critical areas
	n_to_fpm = 2048	;-- gridsize to/from FPM
	n_mft = 1400
	n_from_lyotstop = 4096 
	field_stop_radius_lam0 = 20.0
endif else begin
	use_fpm = 0
	use_lyot_stop = 0
	pupil_diam_pix = 300.0
	n_default = 1024
	n_to_fpm = 1024
	n_from_lyotstop = 1024
endelse


if ( n_elements(optval) ne 0 ) then begin
	if ( tag_exists('lam0',optval) ) then lambda0_m = optval.lam0 * 1.0e-6
	if ( tag_exists('lambda0_m',optval) ) then lambda0_m = optval.lambda0_m
	mas_per_lamD = lambda0_m * 360.0 * 3600.0 / (2 * !dpi * 2.363e6) * 1000	;-- mas per lambda0/D
	if ( tag_exists('source_x_offset',optval) ) then source_x_offset = optval.source_x_offset 
	if ( tag_exists('source_y_offset',optval) ) then source_y_offset = optval.source_y_offset 
	if ( tag_exists('source_x_offset_mas',optval) ) then source_x_offset = optval.source_x_offset_mas / mas_per_lamD 
	if ( tag_exists('source_y_offset_mas',optval) ) then source_y_offset = optval.source_y_offset_mas / mas_per_lamD 
	if ( tag_exists('use_errors',optval) ) then use_errors = optval.use_errors 
	if ( tag_exists('polaxis',optval) ) then polaxis = optval.polaxis 
	if ( tag_exists('zindex',optval) ) then zindex = optval.zindex 
	if ( tag_exists('zval_m',optval) ) then zval_m = optval.zval_m 
        if ( tag_exists('end_at_fsm',optval) ) then end_at_fsm = optval.end_at_fsm
        if ( tag_exists('cgi_x_shift_pupdiam',optval) ) then cgi_x_shift_pupdiam = optval.cgi_x_shift_pupdiam
        if ( tag_exists('cgi_y_shift_pupdiam',optval) ) then cgi_y_shift_pupdiam = optval.cgi_y_shift_pupdiam
	if ( tag_exists('fsm_x_offset',optval) ) then fsm_x_offset = optval.fsm_x_offset
	if ( tag_exists('fsm_y_offset',optval) ) then fsm_y_offset = optval.fsm_y_offset
	if ( tag_exists('fsm_x_offset_mas',optval) ) then fsm_x_offset = optval.fsm_x_offset_mas / mas_per_lamD
	if ( tag_exists('fsm_y_offset_mas',optval) ) then fsm_y_offset = optval.fsm_y_offset_mas / mas_per_lamD
	if ( tag_exists('focm_z_shift_m',optval) ) then focm_z_shift_m = optval.focm_z_shift_m
	if ( tag_exists('use_hlc_dm_patterns',optval) ) then use_hlc_dm_patterns = optval.use_hlc_dm_patterns
	if ( tag_exists('use_dm1',optval) ) then use_dm1 = optval.use_dm1
	if ( tag_exists('dm1_m',optval) ) then dm1_m = optval.dm1_m
	if ( tag_exists('dm1_xc_act',optval) ) then dm1_xc_act = optval.dm1_xc_act
	if ( tag_exists('dm1_yc_act',optval) ) then dm1_yc_act = optval.dm1_yc_act
	if ( tag_exists('dm1_xtilt_deg',optval) ) then dm1_xtilt_deg = optval.dm1_xtilt_deg
	if ( tag_exists('dm1_ytilt_deg',optval) ) then dm1_ytilt_deg = optval.dm1_ytilt_deg
	if ( tag_exists('dm1_ztilt_deg',optval) ) then dm1_ztilt_deg = optval.dm1_ztilt_deg
	if ( tag_exists('use_dm2',optval) ) then use_dm2 = optval.use_dm2
	if ( tag_exists('dm2_m',optval) ) then dm2_m = optval.dm2_m
	if ( tag_exists('dm2_xc_act',optval) ) then dm2_xc_act = optval.dm2_xc_act
	if ( tag_exists('dm2_yc_act',optval) ) then dm2_yc_act = optval.dm2_yc_act
	if ( tag_exists('dm2_xtilt_deg',optval) ) then dm2_xtilt_deg = optval.dm2_xtilt_deg
	if ( tag_exists('dm2_ytilt_deg',optval) ) then dm2_ytilt_deg = optval.dm2_ytilt_deg
	if ( tag_exists('dm2_ztilt_deg',optval) ) then dm2_ztilt_deg = optval.dm2_ztilt_deg
        if ( tag_exists('mask_x_shift_pupdiam',optval) ) then mask_x_shift_pupdiam = optval.mask_x_shift_pupdiam
        if ( tag_exists('mask_y_shift_pupdiam',optval) ) then mask_y_shift_pupdiam = optval.mask_y_shift_pupdiam
	if ( tag_exists('fpm_x_offset',optval) ) then fpm_x_offset = optval.fpm_x_offset
	if ( tag_exists('fpm_y_offset',optval) ) then fpm_y_offset = optval.fpm_y_offset
	if ( tag_exists('fpm_z_shift_m',optval) ) then fpm_z_shift_m = optval.fpm_z_shift_m
	if ( tag_exists('pinhole_diam_m',optval) ) then pinhole_diam_m = optval.pinhole_diam_m 
	if ( tag_exists('end_at_fpm_exit_pupil',optval) ) then end_at_fpm_exit_pupil = optval.end_at_fpm_exit_pupil
	if ( tag_exists('output_field_rootname',optval) ) then output_field_rootname = optval.output_field_rootname
	if ( tag_exists('use_lyot_stop',optval) ) then use_lyot_stop = optval.use_lyot_stop
        if ( tag_exists('lyot_x_shift_pupdiam',optval) ) then lyot_x_shift_pupdiam = optval.lyot_x_shift_pupdiam
        if ( tag_exists('lyot_y_shift_pupdiam',optval) ) then lyot_y_shift_pupdiam = optval.lyot_y_shift_pupdiam
	if ( tag_exists('use_field_stop',optval) ) then use_field_stop = optval.use_field_stop
	if ( tag_exists('field_stop_radius_lam0',optval) ) then field_stop_radius_lam0 = optval.field_stop_radius_lam0
	if ( tag_exists('use_pupil_lens',optval) ) then use_pupil_lens = optval.use_pupil_lens
	if ( tag_exists('use_defocus_lens',optval) ) then use_defocus_lens = optval.use_defocus_lens
	if ( tag_exists('defocus',optval) ) then defocus = optval.defocus
	if ( tag_exists('final_sampling_lam0',optval) ) then final_sampling_lam0 = optval.final_sampling_lam0
endif


diam = 2.3633372d
  fl_pri = 2.83459423440d * 1.0013d
d_pri_sec = 2.285150515460035d
  d_focus_sec = d_pri_sec - fl_pri 
  fl_sec = -0.653933011d * 1.0004095d 
  d_sec_focus = 3.580188916677103d 	; 1 / (1/fl_sec - 1/(d_pri_sec-fl_pri)) 
  diam_sec = 0.58166d
d_sec_fold1 = 2.993753476654728d
  d_fold1_focus = 0.586435440022375d	; d_sec_focus - d_sec_fold1
  diam_fold1 = 0.09
d_fold1_m3 = 1.680935841598811d
  fl_m3 = 0.430216463069001d
  d_focus_m3 = 1.094500401576436d	; d_fold1_m3 - d_fold1_focus
  d_m3_pupil = 0.469156807701977d	; double_lens(d_pri_sec,fl_sec,d_sec_fold1+d_fold1_m3,fl_m3) 
  d_m3_focus = 0.708841602661368d	; 1 / (1/fl_m3 - 1/d_focus_m3)
  diam_m3 = 0.2
d_m3_m4 = 0.943514749358944d
  fl_m4 = 0.116239114833590d
  d_focus_m4 = 0.234673014520402d	; d_m3_m4 - d_m3_focus
  d_pupil_m4 = 0.474357941656967d	; d_m3_m4 - d_m3_pupil
  d_m4_focus = 0.230324117970585d	; 1 / (1/fl_m4 - 1/d_focus_m4)
  diam_m4 = 0.07
d_m4_m5 = 0.429145636743193d
  d_m5_focus = 0.198821518772608d	; d_m4_m5 - d_m4_focus
  fl_m5 = 0.198821518772608d		; d_m5_focus
  d_m5_pupil = 0.716529242882632d	; double_lens(d_pupil_m4,fl_m4,d_m4_m5,fl_m5)
  diam_m5 = 0.07
d_m5_fold2 = 0.351125431220770d
  diam_fold2 = 0.06
d_fold2_fsm = 0.365403811661862d	; d_m5_pupil - d_m5_fold2 ;0.390956443540524d
d_fsm_oap1 = 0.354826767220001d
  fl_oap1 = 0.503331895563883d
  diam_oap1 = 0.06
d_oap1_focm = 0.768005607094041d
d_focm_oap2 = 0.314483210543378d
  fl_oap2 = 0.579156922073536d
  diam_oap2 = 0.06
d_oap2_dm1 = 0.775775726154228d		; double_lens(d_fsm_oap1,fl_oap1,d_oap1_focm+d_focm_oap2,fl_oap2)
d_dm1_dm2 = 1.0d
d_dm2_oap3 = 0.394833855161549d
  fl_oap3 = 1.217276467668519d
  diam_oap3 = 0.06
d_oap3_fold3 = 0.505329955078121d
  diam_fold3 = 0.06
d_fold3_oap4 = 1.158897671642761d
  fl_oap4 = 0.446951159052363d
  diam_oap4 = 0.06
d_oap4_pupilmask = 0.423013568764728d	; double_lens(d_dm1_dm2+d_dm2_oap3,fl_oap3,d_oap3_fold3+d_fold3_oap4,fl_oap4)
d_pupilmask_oap5 = 0.408810648253099d 	; 0.548189351937178d
  fl_oap5 =  0.548189351937178d
  diam_oap5 = 0.06
d_oap5_fpm = 0.548189083164429d   
d_fpm_oap6 = 0.548189083164429d		; d_oap5_fpm
  fl_oap6 = 0.548189083164429d		; d_fpm_oap6
  diam_oap6 = 0.06
d_oap6_lyotstop = 0.687567667550736d	; double_lens(d_pupilmask_oap5,fl_oap5,d_oap5_fpm+d_fpm_oap6,fl_oap6)
d_lyotstop_oap7 = 0.401748843470518d
  fl_oap7 = 0.708251083480054d
  diam_oap7 = 0.06
d_oap7_fieldstop = 0.708251083480054d	; fl_oap7
d_fieldstop_oap8 = 0.210985967281651d
  fl_oap8 = 0.210985967281651d		; d_fieldstop_oap8
  diam_oap8 = 0.06
  d_oap8_pupil = 0.238185804200797d	; double_lens( d_lyotstop_oap7, fl_oap7, d_oap7_fieldstop+d_fieldstop_oap8, fl_oap8 )
d_oap8_filter = 0.368452268225530d
  diam_filter = 0.01
d_filter_lens = 0.170799548215162d
  fl_lens = 0.246017378417573d + 0.050001306014153d
  diam_lens = 0.01
d_lens_fold4 = 0.246017378417573d
  diam_fold4 = 0.02
d_fold4_image = 0.050001578514650d

;s1 = d_oap8_filter + d_filter_lens - d_oap8_pupil
;s2 = d_lens_fold4 + d_fold4_image
fl_pupillens = 0.149260576823040d	; 1 / (1/s1 + 1/s2)

n = n_default	;-- start off with less padding
 
prop_begin, wavefront, diam, lambda_m, n, double(pupil_diam_pix)/n
   if ( cor_type ne 'none' ) then begin 
	fits_read, pupil_file, pupil
   	prop_multiply, wavefront, trim(pupil,n)
   	pupil = 0
   endif else begin
	prop_circular_aperture, wavefront, diam/2
   endelse
   if ( polaxis ne 0 ) then polmap, wavefront, polfile, pupil_diam_pix, polaxis
   prop_define_entrance, wavefront
   prop_lens, wavefront, fl_pri, 'primary'
   if ( source_x_offset ne 0.0 or source_y_offset ne 0.0) then begin
        ;-- compute tilted wavefront to offset source by source_x_offset,source_y_offset lambda0_m/D
        xtilt_lam = -source_x_offset * lambda0_m / lambda_m
        ytilt_lam = -source_y_offset * lambda0_m / lambda_m
        x = ((dindgen(n)-n/2)/(pupil_diam_pix/2)) # replicate(1.0d,n)
        y = transpose(x)
        prop_multiply, wavefront, exp(dcomplex(0,1) * !dpi * (xtilt_lam * x + ytilt_lam * y))
        x = 0
        y = 0
   endif
   if ( zindex[0] ne 0 ) then prop_zernikes, wavefront, zindex, zval_m
   if ( use_errors ) then begin
	prop_errormap, wavefront, map_dir+'wfirst_phaseb_PRIMARY_phase_error_V1.0.fits', /WAVEFRONT
 	prop_errormap, wavefront, map_dir+'wfirst_phaseb_GROUND_TO_ORBIT_4.2X_phase_error_V1.0.fits', /WAVEFRONT
   endif

prop_propagate, wavefront, d_pri_sec, 'secondary'
   prop_lens, wavefront, fl_sec, 'secondary'
   if ( use_errors ) then prop_errormap, wavefront, map_dir+'wfirst_phaseb_SECONDARY_phase_error_V1.0.fits', /WAVEFRONT
   if ( use_aperture ) then prop_circular_aperture, wavefront, diam_sec/2  

prop_propagate, wavefront, d_sec_fold1, 'FOLD_1'
   if ( use_errors ) then prop_errormap, wavefront, map_dir+'wfirst_phaseb_FOLD1_phase_error_V1.0.fits', /WAVEFRONT
   if ( use_aperture ) then prop_circular_aperture, wavefront, diam_fold1/2 

prop_propagate, wavefront, d_fold1_m3, 'M3'
   prop_lens, wavefront, fl_m3
   if ( use_errors ) then prop_errormap, wavefront, map_dir+'wfirst_phaseb_M3_phase_error_V1.0.fits', /WAVEFRONT
   if ( use_aperture ) then prop_circular_aperture, wavefront, diam_m3/2 

prop_propagate, wavefront, d_m3_m4, 'M4' 
   prop_lens, wavefront, fl_m4
   if ( use_errors ) then prop_errormap, wavefront, map_dir+'wfirst_phaseb_M4_phase_error_V1.0.fits', /WAVEFRONT
   if ( use_aperture ) then prop_circular_aperture, wavefront, diam_m4/2 

prop_propagate, wavefront, d_m4_m5, 'M5' 
   prop_lens, wavefront, fl_m5
   if ( use_errors ) then prop_errormap, wavefront, map_dir+'wfirst_phaseb_M5_phase_error_V1.0.fits', /WAVEFRONT
   if ( use_aperture ) then prop_circular_aperture, wavefront, diam_m5/2 
 
prop_propagate, wavefront, d_m5_fold2, 'FOLD_2'
   if ( use_errors ) then prop_errormap, wavefront, map_dir+'wfirst_phaseb_FOLD2_phase_error_V1.0.fits', /WAVEFRONT
   if ( use_aperture ) then prop_circular_aperture, wavefront, diam_fold2/2

prop_propagate, wavefront, d_fold2_fsm, 'FSM'
   if ( end_at_fsm ) then begin
	prop_end, wavefront, sampling_m, /NOABS
	return
   endif
   if ( cgi_x_shift_pupdiam ne 0 or cgi_y_shift_pupdiam ne 0 ) then begin  ;-- bulk coronagraph pupil shear
        ;-- shears are normalized to pupil diameter
        ;-- FFT the field, apply a tilt, FFT back
        xt = -cgi_x_shift_pupdiam * pupil_diam_pix * double(pupil_diam_pix)/n 
        yt = -cgi_y_shift_pupdiam * pupil_diam_pix * double(pupil_diam_pix)/n
        x = (dindgen(n)-n/2) / (pupil_diam_pix/2) # replicate(1.0d, n)
        y = transpose(x)
        tilt = dcomplex(0,1)*!dpi*(x*xt + y*yt)
        x = 0
        y = 0
        wavefront0 = prop_get_wavefront(wavefront)
        wavefront0 = fftsi(wavefront0,-1)
        wavefront0 = wavefront0 * exp(tilt)
        wavefront0 = fftsi(wavefront0,1)
        wavefront.wavefront = prop_shift_center(wavefront0)
   endif
   if ( use_errors ) then prop_errormap, wavefront, map_dir+'wfirst_phaseb_FSM_phase_error_V1.0.fits', /WAVEFRONT
   if ( use_aperture ) then prop_circular_aperture, wavefront, diam_fsm/2
   if ( fsm_x_offset ne 0.0 or fsm_y_offset ne 0.0) then begin
        ;-- compute tilted wavefront to offset source by fsm_x_offset,fsm_y_offset lambda0_m/D
        xtilt_lam = fsm_x_offset * lambda0_m / lambda_m
        ytilt_lam = fsm_y_offset * lambda0_m / lambda_m
        x = ((dindgen(n)-n/2)/(pupil_diam_pix/2)) # replicate(1.0d,n)
        y = transpose(x)
        prop_multiply, wavefront, exp(dcomplex(0,1) * !dpi * (xtilt_lam * x + ytilt_lam * y))
        x = 0
        y = 0
   endif

prop_propagate, wavefront, d_fsm_oap1, 'OAP1'
   prop_lens, wavefront, fl_oap1
   if ( use_errors ) then prop_errormap, wavefront, map_dir+'wfirst_phaseb_OAP1_phase_error_V1.0.fits', /WAVEFRONT
   if ( use_aperture ) then prop_circular_aperture, wavefront, diam_oap1/2  

prop_propagate, wavefront, d_oap1_focm+focm_z_shift_m, 'FOCM'
   if ( use_errors ) then prop_errormap, wavefront, map_dir+'wfirst_phaseb_FOCM_phase_error_V1.0.fits', /WAVEFRONT
   if ( use_aperture ) then prop_circular_aperture, wavefront, diam_focm/2

prop_propagate, wavefront, d_focm_oap2+focm_z_shift_m, 'OAP2'
   prop_lens, wavefront, fl_oap2
   if ( use_errors ) then prop_errormap, wavefront, map_dir+'wfirst_phaseb_OAP2_phase_error_V1.0.fits', /WAVEFRONT
   if ( use_aperture ) then prop_circular_aperture, wavefront, diam_oap2/2  

prop_propagate, wavefront, d_oap2_dm1, 'DM1'
   if ( use_dm1 ) then prop_dm, wavefront, dm1_m, dm1_xc_act, dm1_yc_act, dm_sampling_m, XTILT=dm1_xtilt_deg, YTILT=dm1_ytilt_deg, ZTILT=dm1_ztilt_deg
   if ( use_errors ) then prop_errormap, wavefront, map_dir+'wfirst_phaseb_DM1_phase_error_V1.0.fits', /WAVEFRONT
   if ( cor_type eq 'hlc' and use_hlc_dm_patterns ) then begin
   	fits_read, prefix+'dm1wfe.fits', dm1wfe
    	prop_add_phase, wavefront, trim(dm1wfe, n)
   	dm1wfe = 0
   endif

prop_propagate, wavefront, d_dm1_dm2, 'DM2'
   if ( use_dm2 ) then prop_dm, wavefront, dm2_m, dm2_xc_act, dm2_yc_act, dm_sampling_m, XTILT=dm2_xtilt_deg, YTILT=dm2_ytilt_deg, ZTILT=dm2_ztilt_deg
   if ( use_errors ) then prop_errormap, wavefront, map_dir+'wfirst_phaseb_DM2_phase_error_V1.0.fits', /WAVEFRONT
   if ( cor_type eq 'hlc' ) then begin
	if ( use_hlc_dm_patterns ) then begin
   		fits_read, prefix+'dm2wfe.fits', dm2wfe
    		prop_add_phase, wavefront, trim(dm2wfe, n)
   		dm2wfe = 0
	endif
   	fits_read, prefix+'dm2mask.fits', dm2mask
    	prop_multiply, wavefront, trim(dm2mask, n)
   	dm2mask = 0
   endif 

prop_propagate, wavefront, d_dm2_oap3, 'OAP3'
   prop_lens, wavefront, fl_oap3
   if ( use_errors ) then prop_errormap, wavefront, map_dir+'wfirst_phaseb_OAP3_phase_error_V1.0.fits', /WAVEFRONT
   if ( use_aperture ) then prop_circular_aperture, wavefront, diam_oap3/2 

prop_propagate, wavefront, d_oap3_fold3, 'FOLD_3'
   if ( use_errors ) then prop_errormap, wavefront, map_dir+'wfirst_phaseb_FOLD3_phase_error_V1.0.fits', /WAVEFRONT
   if ( use_aperture ) then prop_circular_aperture, wavefront, diam_fold3/2

prop_propagate, wavefront, d_fold3_oap4, 'OAP4'
   prop_lens, wavefront, fl_oap4
   if ( use_errors ) then prop_errormap, wavefront, map_dir+'wfirst_phaseb_OAP4_phase_error_V1.0.fits', /WAVEFRONT
   if ( use_aperture ) then prop_circular_aperture, wavefront, diam_oap4/2

prop_propagate, wavefront, d_oap4_pupilmask, 'PUPIL_MASK'	;-- flat/reflective shaped pupil 
   if ( cor_type eq 'spc-ifs' or cor_type eq 'spc-wide' ) then begin
	if ( mask_x_shift_pupdiam ne 0 or mask_y_shift_pupdiam ne 0 ) then begin
        	;-- offsets are normalized to pupil diameter
		;-- shift wavefront, apply SP mask, unshift wavefront by
        	;-- FFTing the field, applying a tilt, FFTing back to mask, applying centered mask,
		;-- then FFTing, taking out tilt, and FFTing back to pupil
        	xt = mask_x_shift_pupdiam * pupil_diam_pix * double(pupil_diam_pix)/n
        	yt = mask_y_shift_pupdiam * pupil_diam_pix * double(pupil_diam_pix)/n
        	x = (dindgen(n)-n/2) / (pupil_diam_pix/2) # replicate(1.0d, n)
        	y = transpose(x)
        	tilt = dcomplex(0,1)*!dpi*(x*xt + y*yt)
        	x = 0
        	y = 0
		wavefront0 = prop_get_wavefront(wavefront)
        	wavefront0 = fftsi(wavefront0,-1)
        	wavefront0 = wavefront0 * exp(tilt)
        	wavefront0 = fftsi(wavefront0,1)
		wavefront.wavefront = prop_shift_center(wavefront0)
		wavefront0 = 0
	endif
	fits_read, pupil_mask_file, pupil_mask
	prop_multiply, wavefront, trim(pupil_mask,n)
	pupil_mask = 0
	if ( mask_x_shift_pupdiam ne 0 or mask_y_shift_pupdiam ne 0 ) then begin
		wavefront0 = prop_get_wavefront(wavefront)
        	wavefront0 = fftsi(wavefront0,-1)
        	wavefront0 = wavefront0 * exp(-tilt)
        	wavefront0 = fftsi(wavefront0,1)
		wavefront.wavefront = prop_shift_center(wavefront0)
		wavefront0 = 0
        	tilt = 0
	endif
   endif
   if ( use_errors ) then prop_errormap, wavefront, map_dir+'wfirst_phaseb_PUPILMASK_phase_error_V1.0.fits', /WAVEFRONT
   ;-- while at a pupil, use more padding to provide better sampling at FPM
   diam = 2 * prop_get_beamradius(wavefront)
   prop_end, wavefront, /NOABS
   n = n_to_fpm
   wavefront0 = trim(wavefront,n)
   prop_begin, wavefront, diam, lambda_m, n, double(pupil_diam_pix)/n
   wavefront.wavefront = prop_shift_center(wavefront0)
   wavefront0 = 0

prop_propagate, wavefront, d_pupilmask_oap5, 'OAP5'
   prop_lens, wavefront, fl_oap5
   if ( use_errors ) then prop_errormap, wavefront, map_dir+'wfirst_phaseb_OAP5_phase_error_V1.0.fits', /WAVEFRONT
   if ( use_aperture ) then prop_circular_aperture, wavefront, diam_oap5/2  

prop_propagate, wavefront, d_oap5_fpm+fpm_z_shift_m, 'FPM', /TO_PLANE
   if ( use_fpm ) then begin
	if ( fpm_x_offset ne 0 or fpm_y_offset ne 0 ) then begin	
		;-- to reproduce X,Y offset of FPM, offset focus, apply FPM, remove offset by
		;-- FFT of wavefront to pupil, apply tilt, FFT back to focus, apply FPM, FFT 
		;-- back to pupil, remove tilt, FFT back to focus
		;-- shifts are specified in lambda0/D
		x_offset_lamD = fpm_x_offset * lambda0_m / lambda_m
		y_offset_lamD = fpm_y_offset * lambda0_m / lambda_m
        	x = (dindgen(n)-n/2) / (pupil_diam_pix/2) # replicate(1.0d, n)
        	y = transpose(x)
        	tilt = dcomplex(0,1)*!dpi*(x*x_offset_lamD + y*y_offset_lamD)
        	x = 0
        	y = 0
		wavefront0 = prop_get_wavefront(wavefront)
        	wavefront0 = fftsi(wavefront0,-1)
        	wavefront0 = wavefront0 * exp(tilt)
        	wavefront0 = fftsi(wavefront0,1)
		wavefront.wavefront = prop_shift_center(wavefront0)
		wavefront0 = 0
	endif
	if ( cor_type eq 'hlc' ) then begin
        	fits_read, occulter_file_r, occ_r
        	fits_read, occulter_file_i, occ_i
		occ = dcomplex(occ_r,occ_i)
          	prop_multiply, wavefront, trim(occ,n) 
		occ_r = 0
		occ_i = 0
		occ = 0
	endif else if ( cor_type eq 'spc-ifs' or cor_type eq 'spc-wide' ) then begin
		;-- super-sample FPM
		wavefront0 = prop_get_wavefront( wavefront )
                wavefront0 = fftsi(wavefront0, +1)                              ;-- to virtual pupil
		wavefront0 = trim(wavefront0, n_mft)
        	fits_read, fpm_file, fpm
        	nfpm = (size(fpm))[1]
        	fpm_sampling_lam = fpm_sampling_lam0 * lambda0_m / lambda_m
        	wavefront0 = mft2(wavefront0, fpm_sampling_lam, pupil_diam_pix, nfpm, -1)         ;-- MFT to highly-sampled focal plane
        	wavefront0 = temporary(wavefront0) * fpm
        	fpm = 0
        	wavefront0 = mft2(wavefront0, fpm_sampling_lam, pupil_diam_pix, n, +1)  ;-- MFT to virtual pupil 
		wavefront0 = fftsi(wavefront0, -1)				;-- back to normally-sampled focal plane
        	wavefront.wavefront = prop_shift_center(wavefront0)
		wavefront0 = 0
	endif
	if ( fpm_x_offset ne 0 or fpm_y_offset ne 0 ) then begin	
		wavefront0 = prop_get_wavefront(wavefront)
        	wavefront0 = fftsi(wavefront0,-1)
        	wavefront0 = wavefront0 * exp(-tilt)
        	wavefront0 = fftsi(wavefront0,1)
		wavefront.wavefront = prop_shift_center(wavefront0)
		wavefront0 = 0
		tilt = 0
	endif
   endif
   if ( pinhole_diam_m ne 0 ) then begin
	;-- "pinhole_diam_m" is pinhole diameter in meters
	dx_m = prop_get_sampling(wavefront)
	dx_pinhole_diam_m = pinhole_diam_m / 101.0		;-- 101 samples across pinhole
	n_out = 105
	m_per_lamD = dx_m * n / double(pupil_diam_pix)		;-- current focal plane sampling in lambda_m/D
	dx_pinhole_lamD = dx_pinhole_diam_m / m_per_lamD	;-- pinhole sampling in lambda_m/D
       	n_in = round(pupil_diam_pix * 1.2)
	wavefront0 = prop_get_wavefront(wavefront)
	wavefront0 = fftsi(wavefront0, +1)			;-- to virtual pupil
	wavefront0 = trim(wavefront0, n_in)
       	wavefront0 = mft2( wavefront0, dx_pinhole_lamD, pupil_diam_pix, n_out, -1 )		;-- MFT to highly-sampled focal plane
       	p = radius(n_out,n_out) * dx_pinhole_diam_m le pinhole_diam_m/2.0
       	wavefront0 = wavefront0 * p
       	wavefront0 = mft2( wavefront0, dx_pinhole_lamD, pupil_diam_pix, n, +1 )	;-- MFT back to virtual pupil
	wavefront0 = fftsi(wavefront0, -1)			;-- back to normally-sampled focal plane
       	wavefront.wavefront = prop_shift_center(wavefront0)
	wavefront0 = 0
   endif

prop_propagate, wavefront, d_fpm_oap6-fpm_z_shift_m, 'OAP6'
   prop_lens, wavefront, fl_oap6
   if ( use_errors and not end_at_fpm_exit_pupil ) then prop_errormap, wavefront, map_dir+'wfirst_phaseb_OAP6_phase_error_V1.0.fits', /WAVEFRONT
   if ( use_aperture ) then prop_circular_aperture, wavefront, diam_oap6/2  

prop_propagate, wavefront, d_oap6_lyotstop, 'LYOT_STOP'
   ;-- while at a pupil, switch grid size (smaller for HLC for better speed, larger for SPC for better sampling)
   diam = 2 * prop_get_beamradius(wavefront)
   prop_end, wavefront, sampling_m, /NOABS
   n = n_from_lyotstop
   wavefront = trim(wavefront,n)
   wavefront0 = wavefront
   if ( output_field_rootname ne '' ) then begin
	lams = string(lambda_m*1d6, form='(f6.4)')
	pols = 'polaxis' + strtrim(string(round(polaxis)),2)
	writefits, output_field_rootname+'_'+lams+'um_'+pols+'_real.fits', double(wavefront)
	writefits, output_field_rootname+'_'+lams+'um_'+pols+'_imag.fits', imaginary(wavefront)
   endif
   if ( end_at_fpm_exit_pupil ) then return
   prop_begin, wavefront, diam, lambda_m, n, double(pupil_diam_pix)/n
   wavefront.wavefront = prop_shift_center(wavefront0)
   wavefront0 = 0

   if ( use_lyot_stop ) then begin
   	if ( lyot_x_shift_pupdiam ne 0 or lyot_y_shift_pupdiam ne 0 ) then begin
        	;-- offsets are normalized to pupil diameter
		;-- apply shift to lyot stop by shifting wavefront, applying Lyot stop, then unshifting wavefront
        	;-- FFT the field, apply a tilt, FFT back
        	xt = lyot_x_shift_pupdiam * pupil_diam_pix * double(pupil_diam_pix)/n
        	yt = lyot_y_shift_pupdiam * pupil_diam_pix * double(pupil_diam_pix)/n
        	x = (dindgen(n)-n/2) / (pupil_diam_pix/2) # replicate(1.0d, n)
        	y = transpose(x)
        	tilt = dcomplex(0,1)*!dpi*(x*xt + y*yt)
        	x = 0
        	y = 0
		wavefront0 = prop_get_wavefront(wavefront)
        	wavefront0 = fftsi(wavefront0,-1)
        	wavefront0 = wavefront0 * exp(tilt)
        	wavefront0 = fftsi(wavefront0,1)
       		wavefront.wavefront = prop_shift_center(wavefront0)
		wavefront0 = 0
   	endif
	fits_read, lyot_stop_file, lyot
	prop_multiply, wavefront, trim(lyot,n)
	lyot = 0
   	if ( lyot_x_shift_pupdiam ne 0 or lyot_y_shift_pupdiam ne 0 ) then begin
		wavefront0 = prop_get_wavefront(wavefront)
        	wavefront0 = fftsi(wavefront0,-1)
        	wavefront0 = wavefront0 * exp(-tilt)
        	wavefront0 = fftsi(wavefront0,1)
       		wavefront.wavefront = prop_shift_center(wavefront0)
		wavefront0 = 0
        	tilt = 0
   	endif
   endif
   if ( use_pupil_lens eq 1 or pinhole_diam_m ne 0 ) then prop_circular_aperture, wavefront, 1.1, /NORM

prop_propagate, wavefront, d_lyotstop_oap7, 'OAP7'
   prop_lens, wavefront, fl_oap7
   if ( use_errors ) then prop_errormap, wavefront, map_dir+'wfirst_phaseb_OAP7_phase_error_V1.0.fits', /WAVEFRONT
   if ( use_aperture ) then prop_circular_aperture, wavefront, diam_oap7/2  

prop_propagate, wavefront, d_oap7_fieldstop, 'FIELD_STOP'
   if ( use_field_stop and (cor_type eq 'hlc' or cor_type eq 'spc-wide') ) then begin
        sampling_lamD = double(pupil_diam_pix) / n      ;-- sampling at focus in lambda_m/D
        stop_radius = field_stop_radius_lam0 / sampling_lamD * (lambda0_m/lambda_m) * prop_get_sampling(wavefront)
        prop_circular_aperture, wavefront, stop_radius
   endif

prop_propagate, wavefront, d_fieldstop_oap8, 'OAP8'
   prop_lens, wavefront, fl_oap8
   if ( use_errors ) then prop_errormap, wavefront, map_dir+'wfirst_phaseb_OAP8_phase_error_V1.0.fits', /WAVEFRONT
   if ( use_aperture ) then prop_circular_aperture, wavefront, diam_oap8/2  

prop_propagate, wavefront, d_oap8_filter, 'filter'
   if ( use_errors ) then prop_errormap, wavefront, map_dir+'wfirst_phaseb_FILTER_phase_error_V1.0.fits', /WAVEFRONT
   if ( use_aperture ) then prop_circular_aperture, wavefront, diam_filter/2 

prop_propagate, wavefront, d_filter_lens, 'LENS'
   if ( use_pupil_lens eq 0 and use_defocus_lens eq 0 and defocus eq 0 ) then begin   
	   ;-- use imaging lens to create normal focus
   	   prop_lens, wavefront, fl_lens
   	   if ( use_errors ) then prop_errormap, wavefront, map_dir+'wfirst_phaseb_LENS_phase_error_V1.0.fits', /WAVEFRONT
   endif else begin			
	   if ( use_pupil_lens ) then begin 
		;-- use pupil imaging lens
   		prop_lens, wavefront, fl_pupillens
   		if ( use_errors ) then prop_errormap, wavefront, map_dir+'wfirst_phaseb_PUPILLENS_phase_error_V1.0.fits', /WAVEFRONT
	   endif else begin
		;-- table is waves P-V @ 575 nm
		z4_pv_waves = [    -8.993, -8.868, -8.539, -8.336, -7.979, -7.461, -6.802, -5.877, -5.030, -3.845, -2.493,  0.000,  3.011,  9.943, 20.414, 28.687, 43.354, 55.956 ]
		fl_defocus_lens = [ 5.000,  3.500,  2.000,  1.600,  1.200,  0.900,  0.700,  0.550,  0.470,  0.400,  0.350,  0.296,  0.260,  0.220,  0.195,  0.185,  0.175, 0.170 ]
		z6_rms_waves = [    0.000,  0.030,  0.030,  0.029,  0.027,  0.026,  0.023,  0.020,  0.017,  0.013,  0.008, -0.002, -0.013, -0.038, -0.076, -0.107, -0.160, -0.206 ]
		if ( use_defocus_lens ne 0 ) then begin	;-- defocus lens (1-4)
			;-- use one of 4 defocusing lenses
			defocus = [ 18.0, 9.0, -4.0, -8.0 ]	;-- waves P-V @ 550
			lens_fl = interpol( fl_defocus_lens, z4_pv_waves, defocus, /SPLINE )	
   			prop_lens, wavefront, lens_fl[use_defocus_lens-1] 
   			if ( use_errors ) then prop_errormap, wavefront, map_dir+'wfirst_phaseb_DEFOCUSLENS' + $
				strtrim(string(use_defocus_lens),2) + '_phase_error_V1.0.fits', /WAVEFRONT
	   	endif else begin
			;-- specify amount of defocus (P-V waves @ 575 nm)
			lens_fl = interpol( fl_defocus_lens, z4_pv_waves, defocus, /SPLINE )	
   			prop_lens, wavefront, lens_fl 
   			if ( use_errors ) then prop_errormap, wavefront, map_dir+'wfirst_phaseb_DEFOCUSLENS1_phase_error_V1.0.fits', /WAVEFRONT
		endelse
	   endelse 
   endelse
   if ( use_aperture ) then prop_circular_aperture, wavefront, diam_lens/2 

prop_propagate, wavefront, d_lens_fold4, 'FOLD_4'
   if ( use_errors ) then prop_errormap, wavefront, map_dir+'wfirst_phaseb_FOLD4_phase_error_V1.0.fits', /WAVEFRONT
   if ( use_aperture ) then prop_circular_aperture, wavefront, diam_fold4/2 

prop_propagate, wavefront, d_fold4_image, 'IMAGE'

prop_end, wavefront, sampling_m, /NOABS

if ( final_sampling_lam0 ne 0 ) then begin
	mag = (double(pupil_diam_pix)/n) / final_sampling_lam0 * (lambda_m/lambda0_m)
	wavefront = prop_magnify( wavefront, mag, output_dim, /AMP_CONSERVE )
	sampling_m = sampling_m / mag
endif else begin
	wavefront = trim(wavefront, output_dim)
endelse

return
end
