pro plotimgs,channum,imname1,imname2,imname3,minscl=minscl,maxscl=maxscl,$
      gsm=gsm,levmin=levmin,levmax=levmax,nlev=nlev,residscl=residscl,xy=xy,$
      nobox=nobox,xyprint=xyprint,xypos=xypos,regularimg=regularimg,$
      xycolor=xycolor,screenscl=screenscl

!P.Color = '000000'xL
!P.Background = 'FFFFFF'xL

if size(screenscl,/type) eq 0 then screenscl=1.0
dim=float(get_screen_size())
fscl=(sqrt(dim[0]^2+dim[1]^2)/1686.56)^0.58*screenscl
fscl=round(fscl*10.)/10.
if size(channum,/type) eq 7 then fscl=1
xcmarg=2./90.
ycmarg=0.1

colmodel=255
coldata=0
coldata2=220

if size(gsm,/type) eq 0 then gsm=2.
if size(levmin,/type) eq 0 then levmin=0.6
if size(levmax,/type) eq 0 then levmax=28.
if size(nlev,/type) eq 0 then nlev=5
if size(maxscl,/type) eq 0 then maxscl=-1
if size(minscl,/type) eq 0 then minscl=0.
if size(residscl,/type) eq 0 then residscl=-1
if size(nobox,/type) eq 0 then nobox=0
xsize=fix(900*fscl)
ysize=fix(400*fscl)
xyim=fix(300*fscl)
psres=1
inches=0

if (size(imname1))[0] eq 2 then im1=imname1 else fits_read,imname1,im1,h1
if (size(imname2))[0] eq 2 then im2=imname2 else fits_read,imname2,im2,h2
if (size(imname3))[0] eq 2 then im3=imname3 else fits_read,imname3,im3,h3
im1=gblurt(im1,gsm,gsm*3)
im1con=congrid(im1,xyim,xyim,/interp)
im2=gblurt(im2,gsm,gsm*3)
im2con=congrid(im2,xyim,xyim,/interp)
im3=gblurt(im3,gsm,gsm*3)
if keyword_set(xy) then begin
    imy=intarr(n_elements(im1[*,0]),n_elements(im1[0,*]))
    imxy[xy[0]:xy[1],xy[2]:xy[3]]=1
    imxycon=congrid(imxy,xyim,xyim)
    imxycon=shift(imxycon,1,1)+shift(imxycon,-1,1)+shift(imxycon,1,-1)+$
          shift(imxycon,-1,-1)
    bb=where(imxycon gt 0 and imxycon lt 4)
endif

xsize0=xsize
ysize0=ysize
if screenscl gt 1.5 then csize=3.5 else csize=1.8
if screenscl gt 1.5 then cth=3 else cth=2
if size(channum,/type) eq 7 then begin
    xsize=7.5
    ysize=3.333
    inches=1
    psres=xsize/900.
    psopen,channum,/color,xsi=xsize,$
          ysize=ysize,/inches,bits_per_pixel=8
    csize=0.7
    cth=3
endif else chan,channum,xsize,ysize
erase

if maxscl eq -1 then sclmax=max(im1-minscl) else sclmax=maxscl
colbar=fltarr(256,40)
for i=0,n_elements(colbar[0,*])-1 do colbar[*,i]=findgen(256)
colbar=congrid(colbar,560*fscl,40*fscl,/interp)
ctcustom
;loadct,39,/silent
;;reverse_ct
tv,colbar,xcmarg*xsize0*psres,ycmarg*ysize0*psres,$
      inches=inches,xsize=560*psres,ysize=40*psres
tek_color
plot,[0],/nodata,xra=[minscl,sclmax],/xst,yra=[0,1],/yst,/noerase,$
      yticknam=replicate(' ',60),xtit='Counts/pixel',yticklen=1e-5,$
      position=[xcmarg,ycmarg,2./3.-xcmarg,2*ycmarg],charsi=csize,charth=cth,$
      xthi=3,ythi=3

im1=(im1-minscl)/sclmax*255.
im2=(im2-minscl)/sclmax*255.
ii=where(im1 gt 255.)
if ii[0] ne -1 then im1[ii]=255.
ii=where(im2 gt 255.)
if ii[0] ne -1 then im2[ii]=255.
ii=where(im1 lt 0.)
if ii[0] ne -1 then im1[ii]=0.
ii=where(im2 lt 0.)
if ii[0] ne -1 then im2[ii]=0.

if not keyword_set(regularimg) then begin

if residscl eq -1 then sclmax=5.*stddev(im3)*2. else sclmax=residscl*2.
minscl=-sclmax/2.
ctplusminus
;loadct,39,/silent
;;reverse_ct
colbar=fltarr(256,40)
for i=0,n_elements(colbar[0,*])-1 do colbar[*,i]=findgen(256)
colbar=congrid(colbar,260*fscl,40*fscl,/interp)
tv,colbar,(2./3.+xcmarg)*xsize0*psres,ycmarg*ysize0*psres,$
      inches=inches,xsize=260*psres,ysize=40*psres
tek_color
plot,[0],/nodata,xra=[minscl,sclmax+minscl],/xst,yra=[0,1],/yst,/noerase,$
      yticknam=replicate(' ',60),xtit='Counts/pixel',yticklen=1e-5,$
      position=[2./3.+xcmarg,ycmarg,1.-xcmarg,2*ycmarg],$
      charsi=csize,charth=cth,xth=3,yth=3

im3=(im3-minscl)/sclmax*255.
ii=where(im3 gt 255.)
if ii[0] ne -1 then im3[ii]=255.
ii=where(im3 lt 0.)
if ii[0] ne -1 then im3[ii]=0.

endif else begin

if residscl eq -1 then sclmax=max(im3-minscl) else sclmax=residscl
colbar=fltarr(256,40)
for i=0,n_elements(colbar[0,*])-1 do colbar[*,i]=findgen(256)
colbar=congrid(colbar,560*fscl,40*fscl,/interp)
ctplusminus
;loadct,39,/silent
;;reverse_ct
tv,colbar,620*psres,40*psres,inches=inches,xsize=260*psres,ysize=40*psres
tek_color
plot,[0],/nodata,xra=[minscl,sclmax],/xst,yra=[0,1],/yst,/noerase,$
      yticknam=replicate(' ',60),xtit='Counts/pixel',yticklen=1e-5,$
      position=[620./xsize0,40./ysize0,880./xsize0,80./ysize0],charsi=csize,charth=cth

im3=(im3-minscl)/sclmax*255.
ii=where(im3 gt 255.)
if ii[0] ne -1 then im3[ii]=255.
ii=where(im3 lt 0.)
if ii[0] ne -1 then im3[ii]=0.

endelse

ctcustom
;loadct,39,/silent
;;reverse_ct

cont=intarr(xyim,xyim)
contd=intarr(xyim,xyim)
lev=(findgen(nlev)*sqrt(levmax-levmin)/(nlev-1.))^2+levmin
for i=0,nlev-1 do begin
    c1=intarr(xyim,xyim)
    c2=intarr(xyim,xyim)
    c3=intarr(xyim,xyim)
    ii=where(im2con ge lev[i])
    if ii[0] ne -1 then begin
        c1[ii]=1
        ii=where(shift(im2con,1,0) ge lev[i])
        if ii[0] ne -1 then c2[ii]=1
        ii=where(c1+c2 eq 1)
        if ii[0] ne -1 then cont[ii]=1
        ii=where(shift(im2con,0,1) ge lev[i])
        if ii[0] ne -1 then c3[ii]=1
        ii=where(c1+c3 eq 1)
        if ii[0] ne -1 then cont[ii]=1
    endif
    c1=intarr(xyim,xyim)
    c2=intarr(xyim,xyim)
    c3=intarr(xyim,xyim)
    ii=where(im1con ge lev[i])
    if ii[0] ne -1 then begin
        c1[ii]=1
        ii=where(shift(im1con,1,0) ge lev[i])
        if ii[0] ne -1 then c2[ii]=1
        ii=where(c1+c2 eq 1)
        if ii[0] ne -1 then contd[ii]=1
        ii=where(shift(im1con,0,1) ge lev[i])
        if ii[0] ne -1 then c3[ii]=1
        ii=where(c1+c3 eq 1)
        if ii[0] ne -1 then contd[ii]=1
    endif
endfor
cont[0,*]=0 & cont[-1,*]=0 & cont[*,0]=0 & cont[*,-1]=0
contd[0,*]=0 & contd[-1,*]=0 & contd[*,0]=0 & contd[*,-1]=0

ii=where(cont eq 1)
dd=where(contd eq 1)
disp1=congrid(im1,xyim,xyim)
disp1[dd]=coldata
disp1[ii]=colmodel
if size(bb,/type) ne 0 then disp1[bb]=225
disp2=congrid(im2,xyim,xyim)
disp2[ii]=colmodel
if size(bb,/type) ne 0 then disp2[bb]=225
disp3=congrid(im3,xyim,xyim)
disp3[dd]=coldata2
disp3[ii]=colmodel
if size(bb,/type) ne 0 then disp3[bb]=225

disp2[0:1,*]=colmodel
disp3[0:1,*]=colmodel

tv,disp1,0*psres,xyim/3*psres,inches=inches,xsize=xyim*psres,ysize=xyim*psres
tv,disp2,xyim*psres,xyim/3*psres,inches=inches,xsize=xyim*psres,ysize=xyim*psres
ctplusminus
tv,disp3,xyim*2*psres,xyim/3*psres,inches=inches,xsize=xyim*psres,ysize=xyim*psres

tek_color

if keyword_set(xyprint) and keyword_set(xypos) then begin
    if not keyword_set(xycolor) then xycolor=255
    for i=0,n_elements(xyprint)-1 do $
          xyouts,xypos[0,i],xypos[1,i],xyprint[i],charth=cth,/normal,$
          charsi=csize,color=xycolor[i]
endif

if size(channum,/type) eq 7 then begin
    device,bits_per_pixel=4
    psclose
endif

!P.Color = 'FFFFFF'xL
!P.Background = '000000'xL

end
