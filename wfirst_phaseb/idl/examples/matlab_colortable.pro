;   Copyright 2019 California Institute of Technology
; ------------------------------------------------------------------

pro matlab_colortable, red, green, blue

x = findgen(256)/255.0 * 2 - 1

red = fltarr(256)
blue = red
green = red

y0 = 0.0 
x0 = -0.75 
y1 = 1 
x1 = -0.25
t = x + 0.5
blue = blue + ((t-x0)*(y1-y0)/(x1-x0) + y0) * (t gt -0.75 and t le -0.25)
t = x - 0.5
red  = red  + ((t-x0)*(y1-y0)/(x1-x0) + y0) * (t gt -0.75 and t le -0.25)
t = x
green= green+ ((t-x0)*(y1-y0)/(x1-x0) + y0) * (t gt -0.75 and t le -0.25)

t = x + 0.5
blue = blue + (t gt -0.25 and t le 0.25)
t = x - 0.5
red = red + (t gt -0.25 and t le 0.25)
t = x
green= green+ (t gt -0.25 and t le 0.25)

y0 = 1.0 
x0 = 0.25 
y1 = 0 
x1 = 0.75
t = x + 0.5
blue = blue + ((t-x0)*(y1-y0)/(x1-x0) + y0) * (t gt 0.25 and t le 0.75)
t = x - 0.5
red  = red  + ((t-x0)*(y1-y0)/(x1-x0) + y0) * (t gt 0.25 and t le 0.75)
t = x
green= green+ ((t-x0)*(y1-y0)/(x1-x0) + y0) * (t gt 0.25 and t le 0.75)

red = fix(red*255)
blue = fix(blue*255)
green = fix(green*255)

return
end 
