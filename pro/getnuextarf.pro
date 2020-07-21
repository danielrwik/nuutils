pro getnuextarf,cldir,obsid,ab,reg,srcimfile,outdir,outcore=outcore, $
      flat=flat,ptsrc=ptsrc,savfile=savfile,boxsize=boxsize

apstopcorr=0

hpa=headfits(cldir+'nu'+obsid+ab+'01_cl.evt',exten=1,/silent)
obsstr=sxpar(hpa,'DATE-OBS')
t=strsplit(obsstr,'T',/extract)
tt=strsplit(t[0],'-',/extract)
ttt=strsplit(t[1],':',/extract)
juldate,[fix(tt),fix(ttt)],mjd

if not keyword_set(outcore) then begin
    blah=strsplit(reg,'.',/extract)
    blah2=strsplit(blah[0],'/',/extract)
    outcore=blah2[n_elements(blah2)-1]
endif
outarf=outcore+ab+'.arf'
if size(outdir,/type) eq 0 then outdir=cldir
if keyword_set(savfile) then begin
    if file_test(savfile) then restore,savfile
endif
if not keyword_set(boxsize) then boxsize=19
if boxsize mod 2 eq 0 then boxsize--
halfbox=boxsize/2  ; integer division on purpose

mask=reg2mask(srcimfile,reg)
ii=where(mask gt 0.5)
ii2d=array_indices(mask,ii)

if keyword_set(ptsrc) then begin
    xpos=mean(ii2d[0,*])
    ypos=mean(ii2d[1,*])
    weight=1.0
endif else begin
    fits_read,srcimfile,srcim
    srcim=float(srcim)
    if keyword_set(flat) then srcim[*,*]=1.0
    temp=srcim
    srcim[*,*]=0.0
    srcim[ii]=temp[ii]/total(temp[ii])
    undefine,xpos,ypos,weight
    i0=min(ii2d[0,*]) & i1=max(ii2d[0,*])
    j0=min(ii2d[1,*]) & j1=max(ii2d[1,*])
    for i=i0,i1 do for j=j0,j1 do $
          if (i-1-i0) mod boxsize eq 0 and (j-1-j0) mod boxsize eq 0 then begin
        npix=total(mask[i-halfbox:i+halfbox,j-halfbox:j+halfbox])
        if npix gt 0.5 then begin
            push,xpos,i
            push,ypos,j
            push,weight,total(srcim[i-halfbox:i+halfbox,j-halfbox:j+halfbox])
        endif
    endif
endelse

arfstr=mrdfits(getcaldbfile('arf',ab,refmjd=mjd),1,harf,/silent)
arf=fltarr(n_elements(arfstr.specresp))
for i=0,n_elements(xpos)-1 do begin
    onearf=getnuarf(cldir,obsid,ab,xpos[i],ypos[i],$
          optax,xoptmin,yoptmin,xoptax,yoptax,xapstop,yapstop,pa,$
          apstopcorr=apstopcorr)
    arf+=onearf*weight[i]
endfor

arfstr.specresp=arf
mwrfits,arfstr,outdir+'/'+outarf,harf,/silent,/create

if keyword_set(savfile) then $
      save,optax,xoptmin,yoptmin,xoptax,yoptax,xapstop,yapstop,pa,$
            filename=savfile

end
