function getnupsf,cldir,obsid,ab,xpos,ypos,energy,$
      optax,xoptmin,yoptmin,pa,psfsmooth=psfsmooth,dirty=dirty,nodev=nodev,$
      meanoffaxis=meanoffaxis

;psfsmooth=1.1

if not keyword_set(nodev) then nodev=0
;nodev=1

plscl = 60.*6.828076e-4
xpsf=325
if size(optax,/type) eq 0 then begin
    gti=mrdfits(cldir+'nu'+obsid+ab+'01_gti.fits',1,gh,/silent)
    oa=mrdfits(cldir+'nu'+obsid+ab+'_oa.fits',1,oh,/silent)
    for i=0L,n_elements(gti.start)-1 do begin
        ii=where(oa.time ge gti[i].start and oa.time lt gti[i].stop)
        push,xoptax,oa[ii].x_oa
        push,yoptax,oa[ii].y_oa
    endfor
    xoptax-=1    ; assume IDL convention for x,y
    yoptax-=1
    xoptmin=floor(min(xoptax))
    yoptmin=floor(min(yoptax))
    nxopt=ceil(max(xoptax))-xoptmin
    nyopt=ceil(max(yoptax))-yoptmin
    optax=intarr(nxopt,nyopt)
    for i=0,nxopt-1 do for j=0,nyopt-1 do begin
        ii=where(floor(xoptax) eq i+xoptmin and floor(yoptax) eq j+yoptmin)
        if ii[0] ne -1 then optax[i,j]=n_elements(ii)
    endfor
endif
iiopt=where(optax gt 0)
;fopt=float(optax)/total(optax*1.)
fopt=float(optax[iiopt])/total(optax[iiopt]*1.)
ii2d=array_indices(optax,iiopt)

offaxis=fix(sqrt((xpos-(ii2d[0,*]+xoptmin))^2+ $
          (ypos-(ii2d[1,*]+yoptmin))^2)*plscl*2.)
meanoffaxis=mean(offaxis)/2.
jj=where(offaxis gt 20)
if jj[0] ne -1 then print,"GETNUPSF: WARNING: "+$
      str(n_elements(jj)/n_elements(iiopt)*100.,format="(F6.2)")+ $
      "% of positions are beyond PSF library (r>9')"
jj=where(offaxis ge 18)
if jj[0] ne -1 then offaxis[jj]=17
phi=round(-atan(ypos-(ii2d[1,*]+yoptmin),xpos-(ii2d[0,*]+xoptmin))*180./!pi)

minoff=min(offaxis)
noff=max(offaxis)-minoff+1
hist=intarr(noff)
for i=0,noff-1 do begin
    ii=where(offaxis ge minoff+i and offaxis lt minoff+i+1)
    hist[i]=n_elements(ii)
endfor
blah=max(hist,offmost)
minphi=min(phi)
nphi=max(phi)-minphi+1
hist=intarr(nphi)
for i=0,nphi-1 do begin
    ii=where(phi ge minphi+i and phi lt minphi+i+1)
    hist[i]=n_elements(ii)
endfor
blah=max(hist,phimost)

psfall=fltarr(xpsf,xpsf,18)
for i=0,17 do psfall[*,*,i]=mrdfits(getcaldbfile('psf',ab,energy),i+1,ph,/silent)
;      refmjd=.mjd),i+1,ph,/silent)

;;;  PSF energies 1:3-4.5, 2:4.5-6, 3:6-8, 4:8-12, 5:12-20, 6:20-80  ;;;
if energy lt 4.5 then estr='1' else if energy lt 6. then estr='2' else $
      if energy lt 8. then estr='3' else if energy lt 12. then estr='4' else $
      if energy lt 20. then estr='5' else estr='6'
if energy lt 4.5 then estr2='3to4.5' else if energy lt 6. then estr2='4.5to6' else $
    if energy lt 8. then estr2='6to8' else if energy lt 12. then estr2='8to12' else $
    if energy lt 20. then estr2='12to20' else estr2='20to79'
;;;  PSF deviants determined at specific angles only ;;;
devpa=[42,77,170]

testpa=devpa
if pa lt devpa[0] then begin
    pp=where(devpa gt 270.)
    if pp[0] ne -1 then testpa[pp]-=360.
endif
;blah=min(abs(testpa-pa),pp)
;strpa='pa'+str(devpa[pp])
;fits_read,getenv('NUSKYBGD_AUXIL')+'/deviant'+ab+'en'+estr+strpa+'.fits',dev
bigger=testpa-pa
smaller=pa-testpa
bb=where(bigger ge 0.)
ss=where(smaller ge 0.)
if ss[0] eq -1 then palist=min(devpa[bb]) else if bb[0] eq -1 then $
      palist=max(devpa[ss]) else palist=[max(devpa[ss]),min(devpa[bb])]
fits_read,getenv('NUSKYBGD_AUXIL')+'/deviant'+ab+'en'+estr+'pa'+$
      str(palist[0])+'.fits',dev1
if n_elements(palist) eq 2 then begin
    fits_read,getenv('NUSKYBGD_AUXIL')+'/deviant'+ab+'en'+estr+'pa'+$
          str(palist[1])+'.fits',dev2
    wpa=(total(abs(pa-palist))-abs(pa-palist))/total(abs(pa-palist))
    dev=dev1*wpa[0]+dev2*wpa[1]
endif else dev=dev1

fits_read,getenv('NUSKYBGD_AUXIL')+'/psfasym'+ab+'en'+estr+'.fits',asym
for i=0,17 do psfall[*,*,i]-=asym
;asym[*,*]=-0.0
ndev=rot(dev,-pa)
dev=fltarr(xpsf,xpsf)
dev[xpsf/2-120:xpsf/2+120,xpsf/2-120:xpsf/2+120]= $
      ndev[xpsf/2-120:xpsf/2+120,xpsf/2-120:xpsf/2+120]
angmask=radialmask(getenv('NUSKYBGD_AUXIL')+'/psfasym'+ab+'en'+estr+'.fits',$
      [xpsf/2,xpsf/2],5.)
x=fltarr(xpsf,xpsf)
y=fltarr(xpsf,xpsf)
for i=0,xpsf-1 do x[*,i]=findgen(xpsf)-xpsf/2
for i=0,xpsf-1 do y[i,*]=findgen(xpsf)-xpsf/2
r=sqrt(x*x+y*y)

x=findgen(xpsf)*2.46/60.
readcol,getenv('NUSKYBGD_AUXIL')+'/coeffs'+ab+estr2+'keV.dat',coeff,/silent
corr=fltarr(xpsf,xpsf)
for f=0,n_elements(coeff)-1 do corr+=coeff[f]*(r*2.46/60.)^f
ii=where(r*2.46/60. gt 6.0)
blah=min(abs(r*2.46/60.-6.0),jj)
corr[ii]=corr[jj]
;corr[*,*]=1.0

asym=rot(asym,120.-pa)
compositepsf=fltarr(xpsf,xpsf)
;if ab eq 'A' then psf0o=psfall[*,*,0] else psf0o=(psfall[*,*,4]+psfall[*,*,5])/2.
psf0o=psfall[*,*,0]  ;-asym
if not nodev then begin
    psf=psfall[*,*,minoff+offmost]
    psf=rot(psf,minphi+phimost)
    psf0=rot(psf0o,minphi+phimost)
    newdev=fltarr(xpsf,xpsf)
    for ang=1,max(angmask) do begin
        angmask[xpsf/2,xpsf/2]=ang
        ii=where(angmask eq ang)
        ii=ii[sort(r[ii])]
        newr=interpol(r[ii],psf0[ii],psf[ii])
        for k=0,n_elements(newr)-1 do begin
            blah=min(abs(newr[k]-r[ii]),kk)
            newdev[ii[k]]=dev[ii[kk]]
        endfor
    endfor
    newdev=gblur(newdev,1.,3.)
endif
for i=0,noff-1 do for j=0,nphi-1 do begin
    jj=where(offaxis eq minoff+i and phi eq minphi+j)
    if jj[0] ne -1 then begin
        psf=psfall[*,*,minoff+i]
;        peak=max(psf,ii)
;        ii2d=array_indices(psf,ii)
;        psf[ii2d[0]-1:ii2d[0]+1,ii2d[1]-1:ii2d[1]+1]= $
;              psf[ii2d[0]-1,ii2d[1]-2]
        psf=rot(psf,minphi+j)
;        psf0=rot(psf0o,minphi+j)
;;        newdev=dev*psf/rot(psf0o,minphi+j)
;        newdev=fltarr(xpsf,xpsf)
;        for ang=1,max(angmask) do begin
;            angmask[xpsf/2,xpsf/2]=ang
;            ii=where(angmask eq ang)
;            ii=ii[sort(r[ii])]
;            newr=interpol(r[ii],psf0[ii],psf[ii])
;            for k=0,n_elements(newr)-1 do begin
;                blah=min(abs(newr[k]-r[ii]),kk)
;                newdev[ii[k]]=dev[ii[kk]]
;            endfor
;        endfor
;        newdev=gblur(newdev,1.,3.)
        if not nodev then psfim=(psf+newdev+asym)*corr else psfim=(psf+asym)*corr
        ii=where(psfim lt 0.)
        if ii[0] ne -1 then psfim[ii]=0.
        compositepsf+=psfim
;        if not nodev then compositepsf+= $
;              (psf+newdev+asym)*corr/total((psf+newdev+asym)*corr)*total(fopt[jj]) $
;        else compositepsf+=(psf+asym)*corr/total((psf+asym)*corr)*total(fopt[jj])
;        compositepsf+=(psf+newdev)/total(psf+newdev)*total(fopt[jj])
    endif
endfor
ii=where(r gt 160)
compositepsf[ii]=0.


;compositepsf=fltarr(n_elements(psfall[*,0,0]),n_elements(psfall[0,*,0]))
;undefine,phiall,offaxisall
;for k=0,n_elements(iiopt)-1 do begin
;    ii2d=array_indices(optax,iiopt[k])
;    offaxis=fix(sqrt((xpos-(ii2d[0]+xoptmin))^2+ $
;          (ypos-(ii2d[1]+yoptmin))^2)*plscl*2.)
;    if offaxis gt 20 then print,"WARNING: position outside PSF library (r<9')"
;    if offaxis ge 18 then offaxis=17
;    phi=-atan(ypos-(ii2d[1]+yoptmin),xpos-(ii2d[0]+xoptmin))*180./!pi
;    push,phiall,phi
;    push,offaxisall,offaxis
;    temppsf=psfall[*,*,offaxis]
;    compositepsf+=rot(temppsf,phi)*fopt[iiopt[k]]
;endfor
;
;phi=mean(phiall)
;offaxis=median(offaxisall)
;if ab eq 'A' then devoff=0 else devoff=5
;if offaxis ne devoff then begin
;    if ab eq 'A' then psf0=rot(psfall[*,*,0],phi) else $
;          psf0=rot((psfall[*,*,4]+psfall[*,*,5])/2.,phi)
;    psf=rot(psfall[*,*,offaxis],phi)
;    fits_read,getenv('NUSKYBGD_AUXIL')+'/deviant'+ab+'en'+estr+'.fits',dev
;    ndev=rot(dev,-pa)
;    dev=fltarr(xpsf,xpsf)
;    dev[xpsf/2-70:xpsf/2+70,xpsf/2-70:xpsf/2+70]=ndev[xpsf/2-70:xpsf/2+70,xpsf/2-70:xpsf/2+70]
;    angmask=radialmask(getenv('NUSKYBGD_AUXIL')+'/deviant'+ab+'en'+estr+'.fits',$
;          [xpsf/2,xpsf/2],1.)
;    x=fltarr(xpsf,xpsf)
;    y=fltarr(xpsf,xpsf)
;    for i=0,xpsf-1 do x[*,i]=findgen(xpsf)-xpsf/2.
;    for i=0,xpsf-1 do y[i,*]=findgen(xpsf)-xpsf/2.
;    r=sqrt(x*x+y*y)
;    newdev=fltarr(xpsf,xpsf)
;    for ang=1,360 do begin
;        ii=where(angmask eq ang)
;        ii=ii[sort(r[ii])]
;        newr=interpol(r[ii],psf0[ii],psf[ii])
;        for i=0,n_elements(newr)-1 do begin
;            blah=min(abs(newr[i]-r[ii]),jj)
;            newdev[ii[i]]=dev[ii[jj]]
;        endfor
;    endfor
;    newdev=gblur(newdev,1.,3.)
;endif else begin
;    fits_read,getenv('NUSKYBGD_AUXIL')+'/deviant'+ab+'en'+estr+'.fits',dev
;    ndev=rot(dev,-pa)
;    dev=fltarr(xpsf,xpsf)
;    dev[xpsf/2-70:xpsf/2+70,xpsf/2-70:xpsf/2+70]=ndev[xpsf/2-70:xpsf/2+70,xpsf/2-70:xpsf/2+70]
;    newdev=gblur(dev,1.,3.)
;endelse
;dev=fltarr(xpsf,xpsf)
;dev[xpsf/2-70:xpsf/2+70,xpsf/2-70:xpsf/2+70]=newdev[xpsf/2-70:xpsf/2+70,xpsf/2-70:xpsf/2+70]
;newdev=dev
;
;;print,offaxis/2.,total(newdev),total(compositepsf),total(newdev)/total(compositepsf),format='(F5.2,E10.2,E10.2,F8.4)'
;;i=0
;;while file_test('80002092_SN2014J/80002092002/psf'+ab+str(i)+'.fits') do i++
;;fits_write,'80002092_SN2014J/80002092002/psf'+ab+str(i)+'.fits',newdev
;
;fits_write,'temp'+ab+'no.fits',compositepsf
;compositepsf=compositepsf+newdev
;fits_write,'temp'+ab+'.fits',compositepsf

;peak=max(compositepsf,ii)
;ii2d=array_indices(compositepsf,ii)
;compositepsf[ii2d[0]-1:ii2d[0]+1,ii2d[1]-1:ii2d[1]+1]= $
;      compositepsf[ii2d[0]-1,ii2d[1]-1]

if keyword_set(psfsmooth) then $
      compositepsf=gblur(compositepsf,psfsmooth,psfsmooth*3)
compositepsf/=total(compositepsf)

return,compositepsf

end
