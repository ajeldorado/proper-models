pro level_bar, xbar, ybar, dxbar, dybar, minbar, maxbar, ticks, labels, BIT8=BIT8, TEXTCOLOR=textcolor

if ( dxbar gt dybar ) then begin
	xmarg1 = 5
	xmarg2 = 5
	ymarg1 = 15 
	ymarg2 = 15 
	bar = bytscl(indgen(dxbar)) # replicate(1,dybar)
endif else begin
	xmarg1 = 2
	xmarg2 = 32
	ymarg1 = 7
	ymarg2 = 7
	bar = replicate(1,dxbar) # bytscl(indgen(dybar))
endelse

tv, bytarr(xmarg1+dxbar+xmarg2,ymarg1+dybar+ymarg2), xbar, ybar

tvlct, rr, gg, bb, /get
if ( not keyword_set(BIT8) ) then begin
	tv3, rr(bar), gg(bar), bb(bar), xbar+xmarg1, ybar+ymarg1
endif else begin
	tv, bar, xbar+xmarg1, ybar+ymarg1
endelse

device, get_current_font=old_font
;device, set_font='-adobe-new century schoolbook-medium-r-normal--12-120-75-75-p-70-iso8859-1'

for i = 0, n_elements(ticks)-1 do begin
	c = (ticks(i) - float(minbar)) / (maxbar - minbar) 
	if ( dxbar gt dybar ) then begin
		xtick = c * (dxbar-1) + xbar + xmarg1
		if ( n_elements(textcolor) eq 0 ) then begin
			plots, [xtick,xtick], [0,-3]+ybar+ymarg1, /dev
			xyouts, xtick, ybar+ymarg1-15, labels(i), align=0.5, /dev, font=0
		endif else begin
			plots, [xtick,xtick], [0,-3]+ybar+ymarg1, /dev, color=textcolor
			xyouts, xtick, ybar+ymarg1-15, labels(i), align=0.5, /dev, font=0, color=textcolor
		endelse
	endif else begin 
		ytick = c * (dybar-1) + ybar + ymarg1
		if ( n_elements(textcolor) eq 0 ) then begin
			plots, [0,2]+xmarg1+xbar+dxbar, [ytick,ytick], /dev
			xyouts, xmarg1+xbar+dxbar+8, ytick-4, labels(i), /dev, font=0
		endif else begin
			plots, [0,2]+xmarg1+xbar+dxbar, [ytick,ytick], /dev, color=textcolor
			xyouts, xmarg1+xbar+dxbar+8, ytick-4, labels(i), /dev, font=0, color=textcolor
		endelse
	endelse
endfor

if ( n_elements(ticks) eq 2 ) then begin
	y = alog10(findgen(10)+1) * (dybar-1) + ybar + ymarg1
	for i = 1, 8 do plots, [0,2]+xmarg1+xbar+dxbar, [1,1]*y(i), /dev
endif
 
device, set_font=old_font

return
end
