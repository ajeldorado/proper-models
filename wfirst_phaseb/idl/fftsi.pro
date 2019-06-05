;   Copyright 2019 California Institute of Technology
; ------------------------------------------------------------------


function fftsi, data, direction, NTHREADS=NTHREADS

if ( n_elements(direction) eq 0 ) then direction=-1

s = size(data)
ndim = s(0)
nx = s(1)
if ( ndim gt 1 ) then ny = s(2)

a = shift(data, -nx/2, -ny/2)

if ( n_elements(NTHREADS) eq 0 ) then ffti, a, direction else ffti, a, direction, NTHREADS=NTHREADS

a = shift(a, nx/2, ny/2)

return, a
end
