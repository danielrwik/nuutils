pro ctplusminus

;          white lblue blue black   red yellow white
x_sample = [  0,   42,   85,  127,  170,  212,  255]
r_sample = [255,    0,    0,    0,  225,  255,  255]
g_sample = [255,  255,    0,    0,    0,  255,  255]
b_sample = [255,  255,  255,    0,    0,    0,  255]
 
x = findgen(256)
r = interpol( r_sample, x_sample, x )
g = interpol( g_sample, x_sample, x )
b = interpol( b_sample, x_sample, x )
tvlct, r, g, b

end
