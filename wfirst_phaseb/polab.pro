pro polab, polfile, lambda_m, pupil_diam_pix, condition, amp, pha

;-- polfile: rootname of file containing polarization coefficients
;-- lambda_m: wavelength in meters
;-- pupil_diam_pix: diameter of pupil in pixels
;-- condition: polarization circumstance:
;--		-2: -45 deg in, Y out
;--		-1: -45 deg in, X out
;--		 1: +45 deg in, X out
;--		 2: +45 deg in, Y out
;-- amp, pha: returned aberration maps (pha is WFE in meters)

if ( abs(condition) eq 1 ) then dir_out = 0 else dir_out = 1  ;-- dir_out: output polarization (1=X, 2=Y)
if ( condition lt 0 ) then dir_in = 0 else dir_in = 1	      ;-- dir_in: input polarization (negative=-45d, positive=+45d)

;-- zernike coefficient files are [nzer, nlam, ndir_in, ndir_out]
;--	nzer = 22 (number of zernikes)
;--	nlam = 6 or 11 (450 - 950 nm in 100 or 50 nm steps)
;--	ndir_in = 2 (input polarization direction, 0=-45 deg, 1=+45 deg)
;--	ndir_out = 2 (output polarization direction, 0=X, 1=Y)

fits_read, polfile+'_amp.fits', zamp_array
fits_read, polfile+'_pha.fits', zpha_array
s = size(zamp_array)
if ( s[2] eq 6 ) then lam_array_m = (indgen(6) * 100 + 450) * 1.0e-9 else lam_array_m = (indgen(11) * 50 + 450) * 1.0e-9

;-- interpolate to get zernikes at specified wavelength

zamp = dblarr(22)
zpha = dblarr(22)
for iz = 0, 21 do begin 
	zamp[iz] = spline( lam_array_m, reform(zamp_array[iz,*,dir_in,dir_out]), lambda_m )
	zpha[iz] = spline( lam_array_m, reform(zpha_array[iz,*,dir_in,dir_out]), lambda_m )
endfor

n = round(pupil_diam_pix * 1.1)
n = fix(n / 2) * 2 	;-- force even 
x = (dindgen(n) - n/2) / (pupil_diam_pix / 2.)

amp = dblarr(n,n)
pha = dblarr(n,n)

for j = 0, n-1 do begin
	y = x[j]
	r2 = x^2 + y^2
	r = sqrt(r2)
	r3 = r^3
	r4 = r^4
	r5 = r^5
	r6 = r^6
	t = atan(y,x)

	for itype = 0, 1 do begin		;-- 0 = amp, 1 = phase
		map = dblarr(n)

		if ( itype eq 0 ) then begin
			z = zamp
			map = map + z[0]	;-- include piston if amplitude map
		endif else begin
			z = zpha
		endelse

		map = map + z[1] * 2 * x				;-- x tilt
		map = map + z[2] * 2 * y				;-- y tilt
		map = map + z[3] * sqrt(3) * (2*r2 - 1)			;-- focus
		map = map + z[4] * sqrt(6) * r2 * sin(2*t)		;-- 45 deg astig
		map = map + z[5] * sqrt(6) * r2 * cos(2*t)		;-- 0 deg astig
		map = map + z[6] * sqrt(8) * (3*r3 - 2*r) * sin(t)  	;-- y coma
		map = map + z[7] * sqrt(8) * (3*r3 - 2*r) * cos(t)	;-- x coma
		map = map + z[8] * sqrt(8) * r3 * sin(3*t)		;-- y trefoil 
		map = map + z[9] * sqrt(8) * r3 * cos(3*t)		;-- x trefoil 
		map = map + z[10] * sqrt(5) * (6*r4 - 6*r2 + 1)		;-- spherical
      		map = map + z[11] * sqrt(10) * (4*r4 - 3*r2) * cos(2*t)
      		map = map + z[12] * sqrt(10) * (4*r4 - 3*r2) * sin(2*t)
      		map = map + z[13] * sqrt(10) * r4 * cos(4*t)
      		map = map + z[14] * sqrt(10) * r4 * sin(4*t)
      		map = map + z[15] * sqrt(12) * (10*r5 - 12*r3 + 3*r) * cos(t)
      		map = map + z[16] * sqrt(12) * (10*r5 - 12*r3 + 3*r) * sin(t)
      		map = map + z[17] * sqrt(12) * (5*r5 - 4*r3) * cos(3*t)
      		map = map + z[18] * sqrt(12) * (5*r5 - 4*r3) * sin(3*t)
      		map = map + z[19] * sqrt(12) * r5 * cos(5*t)
      		map = map + z[20] * sqrt(12) * r5 * sin(5*t)
      		map = map + z[21] * sqrt(7) * (20*r6 - 30*r4 + 12*r2 - 1)

		if ( itype eq 0 ) then amp[0,j] = map else pha[0,j] = map
	endfor
endfor

return
end
