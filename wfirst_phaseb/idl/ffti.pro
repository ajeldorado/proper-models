;   Copyright 2019 California Institute of Technology
; ------------------------------------------------------------------


;++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
; Name:
;	ffti
;
; Purpose:
;	Compute FFT of wavefront array using Intel FFT library routines
;
; Calling sequence:
;	ffti, array, direction
;
; Inputs:
;	array
;	   2-D data array; this array will be made double-complex prior
;	   to calling the FFT routine, if it is not already
;	direction
;	   Transform direction (-1 for forward, +1 for inverse)
;
; Revision history:
;	Written by John Krist (JPL), February 2006
;----------------------------------------------------------------------------------  
pro ffti, array, direction, NTHREADS=nthreads

if ( n_elements(nthreads) eq 0 ) then nthreads = 0
if ( n_elements(direction) eq 0 ) then direction = -1

s = size(array)

ndim = s(0)
if ( ndim ne 2 ) then begin
	print, 'Error in FFTI: input array is not 2D'
        print,'Stopping.'
	stop
endif

nx = long(s(1))
ny = long(s(2))
datatype = long(s(3))
if ( datatype eq 4 ) then begin
	array = complex(array)
	datatype = 6L
endif else if ( datatype eq 5 ) then begin
	array = dcomplex(array)
	datatype = 9L
endif

stat = call_external( './ffti.so', 'ffti', datatype, array, long(direction), nx, ny, long(nthreads) )

return
end
