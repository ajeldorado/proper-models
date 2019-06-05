;-----------------------------------------------------------------------
; mtf:
;   Compute a matrix fourier transform.  Based on Soummer et al. 2007.
;
; Input parameters:
;    in : 2-D wavefront to transform
;     dout: sampling in lambda/D of output (if pupil-to-focus) or input (if focus-to-pupil)  
;        D: pupil size in pixels
;    nout : dimensions of output array (nout by nout)
;  direction : direction of transform (-1 or +1) 
;
; Optional input parameters:
;   xoffset, yoffset : offsets in output field in cycles/D
;  Returns:
;    2-D Fourier transform of input array.
;
;  Written by Dimitri Mawet (JPL)
;  March 2010
;  Copyright 2012 California Institute of Technology
;------------------------------------------------------------------------
function mft2, in, dout, D, nout, direction, xoffset, yoffset, xc, yc

if ( n_elements(xoffset) eq 0 ) then xoffset = 0.0
if ( n_elements(yoffset) eq 0 ) then yoffset = 0.0
if ( n_elements(xc) eq 0 ) then xc = 0.0
if ( n_elements(yc) eq 0 ) then yc = 0.0

nout = fix(nout)
nin = (size(in))[1]
nin = fix(nin)

x = reform(findgen(nin)-nin/2-xc, 1, nin)
y = findgen(nin)-nin/2-yc 

u = (findgen(nout)-nout/2 - xoffset/dout)*float(dout/D)
v = reform((findgen(nout)-nout/2 - yoffset/dout)*float(dout/D), 1, nout)

if ( direction eq -1 ) then begin
	out=dout/D * exp((-2*!pi*complex(0,1)*u)#x) # in # exp((-2*!pi*complex(0,1)*y)#v)
endif else begin
	out=dout/D * exp((+2*!pi*complex(0,1)*u)#x) # in # exp((+2*!pi*complex(0,1)*y)#v)
endelse

return, out
end
