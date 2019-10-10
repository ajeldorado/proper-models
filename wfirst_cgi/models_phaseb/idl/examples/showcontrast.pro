;   Copyright 2019 California Institute of Technology
; ------------------------------------------------------------------

pro showcontrast, contrast, pupil_ratio, inner_radius, outer_radius, MAG=MAG, xoffset=xoffset, yoffset=yoffset, COLOR_TABLE=COLOR_TABLE, $
	MIN_CONTRAST=min_contrast, MAX_CONTRAST=max_contrast, MATLAB=MATLAB, NOBAR=NOBAR, CIRCULAR=CIRCULAR, C_SHAPED=C_SHAPED, D_SHAPED=D_SHAPED, $
	XBAROFFSET=XBAROFFSET, YBAROFFSET=YBAROFFSET, MASK=mask0, GRIDSIZE=gridsize, TEXTCOLOR=textcolor, OPENING_ANGLE=opening_angle


if ( n_elements(xoffset) eq 0 ) then xoffset=0
if ( n_elements(yoffset) eq 0 ) then yoffset=0

if ( n_elements(min_contrast) eq 0 ) then min_contrast = 1e-11
if ( n_elements(max_contrast) eq 0 ) then max_contrast = 1e-5

if ( n_elements(COLOR_TABLE) eq 0 ) then COLOR_TABLE = 5
s = size(contrast)
nx = s(1)
ny = s(2)
if ( n_elements(GRIDSIZE) eq 0 ) then n = 512 else n = gridsize


if ( not keyword_set(MATLAB) ) then begin
	loadct, COLOR_TABLE
	tvlct,rr,gg,bb,/get
endif else begin
	matlab_colortable, rr,gg,bb
	tvlct,rr,gg,bb
endelse

c = contrast
c = alog10(c>1e-20)

if ( n_elements(mask0) ne 0 ) then begin
	s = size(mask0)
	if ( s(1) ne nx ) then begin
		mask = trim(mask0,nx)
	endif else begin
		mask = mask0
	endelse
endif

if ( keyword_set(MAG) ) then begin
	x = (findgen(n)-n/2) / MAG + nx/2
	y = (findgen(n)-n/2) / MAG + ny/2
	c = interpolate(c,x,y,/GRID,cub=-0.5,miss=0)
	if ( n_elements(mask0) ne 0 ) then mask = interpolate(1.*mask0,x,y,/GRID,miss=0)
endif else begin
	MAG = 1
	c = c(nx/2-n/2:nx/2+n/2-1,ny/2-n/2:ny/2+n/2-1)
	if ( n_elements(mask0) ne 0 ) then mask = mask(nx/2-n/2:nx/2+n/2-1,ny/2-n/2:ny/2+n/2-1)
endelse

lmin = round(alog10(min_contrast))
lmax = round(alog10(max_contrast))
c(0) = lmin
c(1) = lmax
c = c > lmin
c = c < lmax

c = bytscl( c )
if ( n_elements(mask0) ne 0 ) then c = c * mask

tv3, rr(c), gg(c), bb(c), xoffset, yoffset

if ( n_elements(inner_radius) ne 0 ) then begin
	if ( not keyword_set(CIRCULAR) and not keyword_set(D_SHAPED) and not keyword_set(C_SHAPED) ) then begin	;-- rectangular
		plots, [inner_radius,outer_radius,outer_radius,inner_radius,inner_radius]*MAG/pupil_ratio+n/2+xoffset, $
		       [-outer_radius,-outer_radius,outer_radius,outer_radius,-outer_radius]*MAG/pupil_ratio+n/2+yoffset, $
			/dev, thick=2
		plots, -[inner_radius,outer_radius,outer_radius,inner_radius,inner_radius]*MAG/pupil_ratio+n/2+xoffset, $
		       [-outer_radius,-outer_radius,outer_radius,outer_radius,-outer_radius]*MAG/pupil_ratio+n/2+yoffset, $
			/dev, thick=2
	endif else begin
		t = findgen(500)/499 * 2 * !pi
		plots, MAG*outer_radius/pupil_ratio*cos(t)+n/2+xoffset,MAG*outer_radius/pupil_ratio*sin(t)+n/2+yoffset,/dev, thick=3
		if ( keyword_set(D_SHAPED) ) then begin
			;-- linear inner edge
			plots, MAG*[inner_radius,inner_radius]/pupil_ratio+n/2+xoffset,[0,n]+yoffset, /dev, thick=2
			plots, -MAG*[inner_radius,inner_radius]/pupil_ratio+n/2+xoffset,[0,n]+yoffset, /dev, thick=2
		endif else begin
			;-- circular inner edge
   			if ( inner_radius ne 0 ) then begin
			   plots, MAG*inner_radius/pupil_ratio*cos(t)+n/2+xoffset,MAG*inner_radius/pupil_ratio*sin(t)+n/2+yoffset,/dev, thick=3
			endif
		endelse
	endelse
endif

if ( n_elements(opening_angle) ne 0 ) then begin
	y = n/2 * atan(opening_angle/2.0*!pi/180.0)
	plots, [0,n-1]+xoffset, [-y,y]+n/2+yoffset, /dev, thick=3
	plots, [0,n-1]+xoffset, [y,-y]+n/2+yoffset, /dev, thick=3
endif

ticks = indgen(lmax-lmin+1) + lmin
labels = strtrim(string(ticks),2)
if ( not keyword_set(NOBAR) ) then begin
	if ( n_elements(XBAROFFSET) eq 0 ) then xbaroffset = 0
	if ( n_elements(YBAROFFSET) eq 0 ) then ybaroffset = 0
	if ( n_elements(TEXTCOLOR) eq 0 ) then begin
		level_bar, xoffset+xbaroffset, n/2-128+yoffset+ybaroffset, 20, 256, min(ticks), max(ticks), ticks, labels
	endif else begin
		level_bar, xoffset+xbaroffset, n/2-128+yoffset+ybaroffset, 20, 256, min(ticks), max(ticks), ticks, labels, TEXTCOLOR=textcolor
	endelse
endif

empty

return
end
