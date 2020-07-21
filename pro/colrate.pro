function colrate,e1,e2,f,h1,h2

if h2 eq 1.0 then h2=0.99

c = (1.0-h1)/(1.0+h1)
d = (1.0+h2)/(1.0-h2)

if fix(e1) eq 4 and fix(e2) eq 6 then rate=f*c/(1.0+c+d) else $
      if fix(e1) eq 6 and fix(e2) eq 12 then rate=f/(1.0+c+d) else $
      if fix(e1) eq 12 and fix(e2) eq 25 then rate=f*d/(1.0+c+d) else $
      stop,'colrate function: invalid energy range: '+str(e1)+'-'+str(e2)+' keV'

return,rate*1e-3

end
