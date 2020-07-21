function getnuarf,cldir,obsid,ab,xpos,ypos,mjd, $
      optax,xoptmin,yoptmin,xoptax,yoptax,xapstop,yapstop,pa, $
      offang=offang,weight=weight,apstopcorr=apstopcorr

plscl = 60.*6.828076e-4
mmscl = 0.12096
if size(optax,/type) eq 0 then begin
    undefine,xoptmin,yoptmin,xoptax,yoptax,xapstop,yapstop
    gti=mrdfits(cldir+'nu'+obsid+ab+'01_gti.fits',1,gh,/silent)
    hpa=headfits(cldir+'nu'+obsid+ab+'01_cl.evt',exten=1,/silent)
    pa=sxpar(hpa,'PA_PNT')
    aprot=pa+120.+180.
    if aprot ge 360. then aprot=(aprot-360.)*!pi/180. else aprot=aprot*!pi/180.
    aprot*=-1.
    oa=mrdfits(cldir+'nu'+obsid+ab+'_oa.fits',1,oh,/silent)
    for i=0L,n_elements(gti.start)-1 do begin
        ii=where(oa.time ge gti[i].start and oa.time lt gti[i].stop)
        push,xoptax,oa[ii].x_oa
        push,yoptax,oa[ii].y_oa
        push,xapstop,((oa[ii].x_apstop*cos(aprot)-oa[ii].y_apstop*sin(aprot))-$
              (oa[ii].x_oa*cos(aprot)-oa[ii].y_oa*sin(aprot)))*mmscl
        push,yapstop,-((oa[ii].x_apstop*sin(aprot)+oa[ii].y_apstop*cos(aprot))-$
              (oa[ii].x_oa*sin(aprot)+oa[ii].y_oa*cos(aprot)))*mmscl
    endfor

    xoptax-=1    ; assume IDL convention for x,y
    yoptax-=1

    xoptmin=floor(min(xoptax))
    yoptmin=floor(min(yoptax))
    nxopt=ceil(max(xoptax))-xoptmin
    nyopt=ceil(max(yoptax))-yoptmin
    optax=lonarr(nxopt,nyopt)
    for i=0,nxopt-1 do for j=0,nyopt-1 do begin
        ii=where(floor(xoptax) eq i+xoptmin and floor(yoptax) eq j+yoptmin)
        if ii[0] ne -1 then optax[i,j]=n_elements(ii)
    endfor
endif
iiopt=where(optax gt 0)
;fopt=float(optax)/total(optax*1.)
fopt=float(optax[iiopt])/total(optax[iiopt]*1.)
ii2d=array_indices(optax,iiopt)



vignstr=mrdfits(getcaldbfile('vign',ab),$
;refmjd=nuimylze_obs[iobs].mjd),$
      1,hvign,/silent)
vazarr=vignstr.azimuth
vtheta=vignstr[0].theta
arfstr=mrdfits(getcaldbfile('arf',ab,refmjd=mjd),1,harf,/silent)
apststr=mrdfits(getcaldbfile('apstop',ab), $
;,refmjd=nuimylze_obs[iobs].mjd),$
      1,hapstop,/silent)
azarr=apststr.azimuth
phi=atan(ypos-yoptax,xpos-xoptax)*180./!pi
ii=where(phi lt 0.)
if ii[0] ne -1 then phi[ii]+=360.
azval=mean(300.-phi+pa)
if azval gt 360. then azval-=360.
blah=min(abs(azval-vazarr),vaz)
v=vignstr[vaz].vignet
blah=size(v)
maxoff=blah[2]
blah=min(abs(azval-azarr),az)
;print,'AZ: ',azval,vaz,az
;print,'dX: ',mean(xapstop)
;print,'dY: ',mean(yapstop)
apstarr=apststr[az].aperture
apstx=apststr[az].deltax_cen
apsty=apststr[az].deltay_cen
apstoffang=apststr[az].theta
apste1=apststr[az].energ_lo
apste2=apststr[az].energ_hi
apste=(apste1+apste2)/2.
napste=n_elements(apste)
noffang=n_elements(apstoffang)
nxy=n_elements(apstx)

undefine,avgoffaxis
earf=1.62+0.04*indgen(4096)
ee=where(earf lt apste2[napste-1])
vign=fltarr(4096)

if not keyword_set(apstopcorr) then begin

offaxis=sqrt((xpos-(ii2d[0,*]+xoptmin))^2+ $
          (ypos-(ii2d[1,*]+yoptmin))^2)*plscl
offang=total(offaxis*fopt)
for i=0,n_elements(vtheta)-2 do begin
    jj=where(offaxis ge vtheta[i] and offaxis lt vtheta[i+1])
    if jj[0] ne -1 then begin
        slope=(v[*,i+1]-v[*,i])/(vtheta[i+1]-vtheta[i])
        for j=0,n_elements(jj)-1 do $
              vign+=(v[*,i]+slope*(offaxis[jj[j]]-vtheta[i]))*fopt[jj[j]]
    endif
endfor

;offaxis=fix(sqrt((xpos-(ii2d[0,*]+xoptmin))^2+ $
;          (ypos-(ii2d[1,*]+yoptmin))^2)*plscl*6.)
;minoff=min(offaxis)
;noff=max(offaxis)-minoff+1
;if noff+minoff gt maxoff then noff=maxoff-minoff
;for i=0,noff-1 do begin
;    jj=where(offaxis eq minoff+i)
;    if jj[0] ne -1 then begin
;        vign+=v[*,minoff+i]*total(fopt[jj])
;    endif
;endfor
;offang=total(offaxis*fopt)/6.

endif else begin

for i=0,n_elements(iiopt)-1 do begin
    ii2d=array_indices(optax,iiopt[i])
    offaxis=sqrt((xpos-(ii2d[0]+xoptmin))^2+(ypos-(ii2d[1]+yoptmin))^2)*plscl
    blah=min(abs(apstoffang-offaxis),oo)
    push,avgoffaxis,offaxis
    offaxis*=6.
    joffax=floor(offaxis)
    ii=where(floor(xoptax) eq ii2d[0]+xoptmin and floor(yoptax) eq ii2d[1]+yoptmin)
    xas=xapstop[ii]
    yas=yapstop[ii]
    w=0L
    apstspec=fltarr(napste)
    for j=min(round(xas)),max(round(xas)) do $
          for k=min(round(yas)),max(round(yas)) do begin
              aa=where(round(xas) eq j and round(yas) eq k)
              if aa[0] ne -1 then begin
                    blah=min(abs(j-apstx),xx)
                    blah=min(abs(k-apsty),yy)
                    index=long(long((oo*nxy+xx)*nxy)+yy)*long(napste)
                    apstspec+=float(apstarr[index:index+long(napste)-1])*$
                          n_elements(aa)/1000.
                    w+=n_elements(aa)
              endif
          endfor
    apstspec/=float(w)
    apst=fltarr(n_elements(vign))
    apst[*]=1.0
    apst[ee]=interpol(apstspec,apste,earf[ee])
    vign+=(v[*,joffax]*(offaxis-joffax)+v[*,joffax+1]*(joffax+1-offaxis))*$
          fopt[i]*apst
endfor
offang=total(avgoffaxis*fopt)
;print,'OffAng: ',offang

endelse

;210 is 10 keV
;1210 is 50 keV
weight=vign[1210]

return,vign*arfstr.specresp

end
