pro ctcustom

;          black lblue blue yellow orange purple red white
x_sample = [  0,   30,   60,   90,  120,  160,  200,  255]
r_sample = [  0,    0,    0,  255,  255,  153,  255,  255]
g_sample = [  0,  255,    0,  215,  150,    0,    0,  255]
b_sample = [  0,  255,  255,    0,    0,  153,    0,  255]
 
;          black  blue lblue green orange yellow red purple white
x_sample = [  0,   20,   70,   85,  110,  140,  160,  220,  255]
r_sample = [  0,    0,    0,    0,  200,  255,  255,  223,  255]
g_sample = [  0,    0,  255,  215,  120,  215,    0,    0,  255]
b_sample = [  0,  255,  255,    0,    0,    0,    0,  223,  255]

; make orange less brown -- looks better, but is less discriminating
if 1 then begin
    r_sample[4]=255
    g_sample[4]=150
endif
 
x = findgen(256)
r = interpol( r_sample, x_sample, x )
g = interpol( g_sample, x_sample, x )
b = interpol( b_sample, x_sample, x )
tvlct, r, g, b

end
