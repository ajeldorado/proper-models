;   Copyright 2019 California Institute of Technology
; ------------------------------------------------------------------


function angle, nx, ny, xc, yc, RADIANS=RADIANS

;-- returns angle in degrees (or radians if /RADIANS)

x = findgen(nx) # replicate(1.0,ny)
y = replicate(1.0,nx) # findgen(ny)

if ( n_elements(xc) eq 0 ) then xc = nx/2
if ( n_elements(yc) eq 0 ) then yc = ny/2

x = temporary(x) - xc
y = temporary(y) - yc

if ( keyword_set(RADIANS) ) then begin
	return, atan(y,x) 
endif else begin
	return, atan(y,x) * (180/!pi)
endelse

end
