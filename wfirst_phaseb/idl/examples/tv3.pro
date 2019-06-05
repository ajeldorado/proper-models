;   Copyright 2019 California Institute of Technology
; ------------------------------------------------------------------

pro tv3, red, green, blue, xoff, yoff, SCALE=SCALE

if ( n_elements(red) eq 1 ) then begin
	s = size(green)
	red = bytarr(s(1),s(2))
endif else if ( n_elements(green) eq 1 ) then begin
	s = size(red)
	green = bytarr(s(1),s(2))
endif else if ( n_elements(blue) eq 1 ) then begin
	s = size(red)
	blue = bytarr(s(1),s(2))
endif

if ( n_elements(xoff) eq 0 ) then begin
	xoff = 0
	yoff = 0
endif

if ( keyword_set(SCALE) ) then begin
	tv,[[[bytscl(red)]],[[bytscl(green)]],[[bytscl(blue)]]],xoff,yoff,true=3
endif else begin
	tv,[[[red]],[[green]],[[blue]]],xoff,yoff,true=3
endelse

return
end
