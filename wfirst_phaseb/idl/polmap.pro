pro polmap, wavefront, polfile, pupil_diam_pix, condition, MUF=muf

;-- wavefront: current wavefront structure
;-- polfile: rootname of file containing polarization coefficients
;-- pupil_diam_pix: diameter of pupil in pixels
;-- condition: polarization circumstance:
;--		-2: -45 deg in, Y out
;--		-1: -45 deg in, X out
;--		 1: +45 deg in, X out
;--		 2: +45 deg in, Y out
;--		 5: X polarization (mean of +/-45 deg in, X out)
;--		 6: Y polarization (mean of +/-45 deg in, X out)
;--		10: All polarizations (mean of +/-45 deg in, X&Y out)
;--	NOTE: the mean conditions (5,6,10) should only be used for sensing;
;--	contrast evaluation must be done by computing each in/out condition separately

n = prop_get_gridsize()
lambda_m = prop_get_wavelength(wavefront)

if ( n_elements(muf) eq 0 ) then muf = 1.0

if ( condition le 2 ) then begin
	polab, polfile, lambda_m, pupil_diam_pix, condition, amp, pha
endif else if ( condition eq 5 ) then begin
	polab, polfile, lambda_m, pupil_diam_pix, -1, amp_m45_x, pha_m45_x
	polab, polfile, lambda_m, pupil_diam_pix, +1, amp_p45_x, pha_p45_x
	amp = (amp_m45_x + amp_p45_x) / 2
	pha = (pha_m45_x + pha_p45_x) / 2
endif else if ( condition eq 6 ) then begin
	polab, polfile, lambda_m, pupil_diam_pix, -2, amp_m45_y, pha_m45_y
	polab, polfile, lambda_m, pupil_diam_pix, +2, amp_p45_y, pha_p45_y
	amp = (amp_m45_y + amp_p45_y) / 2
	pha = (pha_m45_y + pha_p45_y) / 2
endif else if ( condition eq 10 ) then begin
	polab, polfile, lambda_m, pupil_diam_pix, -1, amp_m45_x, pha_m45_x
	polab, polfile, lambda_m, pupil_diam_pix, +1, amp_p45_x, pha_p45_x
	polab, polfile, lambda_m, pupil_diam_pix, -2, amp_m45_y, pha_m45_y
	polab, polfile, lambda_m, pupil_diam_pix, +2, amp_p45_y, pha_p45_y
	amp = (amp_m45_x + amp_p45_x + amp_m45_y + amp_p45_y) / 4
	pha = (pha_m45_x + pha_p45_x + pha_m45_y + pha_p45_y) / 4
endif else begin
	print, 'POLMAP: unmatched condition = ', condition
	stop
endelse

prop_multiply, wavefront, trim(amp,n) 
prop_add_phase, wavefront, trim(muf*pha,n) 

return
end
