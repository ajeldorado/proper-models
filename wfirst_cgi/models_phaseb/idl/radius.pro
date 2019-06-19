;   Copyright 2019 California Institute of Technology
; ------------------------------------------------------------------


function radius, nx, ny, xc, yc

if ( n_params() eq 2 ) then begin
	xc = fix(nx) / 2
	yc = fix(ny) / 2
endif

x = dindgen(nx) # replicate(1.0,ny)
y = replicate(1.0,nx) # dindgen(ny)
x = temporary(x) - xc
y = temporary(y) - yc
r = sqrt(x*x + y*y)

return, r
end
