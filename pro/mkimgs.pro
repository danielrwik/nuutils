pro mkimgs,cldir,obsid,ab,en1,en2,t1,t2,outname=outname,det=det,$
      colfilt=colfilt,fexp=fexp,suffevt=suffevt

if float(fix(en1)) eq en1 then e1name=str(fix(en1)) else $
      e1name=str(en1,format='(F'+str(floor(alog10(en1))+3)+'.1)')
if float(fix(en2)) eq en2 then e2name=str(fix(en2)) else $
      e2name=str(en2,format='(F'+str(floor(alog10(en2))+3)+'.1)')
;if size(en1,/type) eq 2 then e1name=str(en1) else $
;      e1name=str(string(en1,format='(F6.2)'))
;if size(en2,/type) eq 2 then e2name=str(en2) else $
;      e2name=str(string(en2,format='(F6.2)'))
if not keyword_set(outname) then outname=cldir+'/im'+ab+e1name+'to'+e2name+'keV.fits'

if keyword_set(suffevt) then evtend='01_cl'+suffevt+'.evt' else evtend='01_cl.evt'
if file_test(cldir,/directory) then evtname=cldir+'/nu'+obsid+ab+evtend $
      else if file_test(cldir) then evtname=cldir $
      else stop,'MKIMGS: Failed miserably attempting to find event file'
evts=mrdfits(evtname,1,header,/silent)
head2=headfits(evtname,exten=0,/silent)
if not keyword_set(det) then begin
    x=evts.x
    y=evts.y
    xtt='X'
    ytt='Y'
endif else begin
    x=evts.det1x
    y=evts.det1y
    xtt='DET1X'
    ytt='DET1Y'
endelse

e1=round((en1-1.6)/0.04+1)
e2=round((en2-1.6)/0.04)
gti=mrdfits(evtname,2,/silent)
exptime=0.
if size(t1,/type) eq 0 then begin
    t1=gti[0].start
    t2=gti[n_elements(gti)-1].stop
endif
if size(colfilt,/type) eq 0 then colfilt=''
if size(colfilt,/type) ne 7 then if colfilt ne '' then colfilt=' and ('+colfilt+')'

undefine,ii
for t=0,n_elements(t1)-1 do begin
    check=execute('thisii=where(evts.pi ge e1 and evts.pi lt e2 and '+$
          'evts.time ge t1[t] and evts.time lt t2[t] and '+$
          'evts.grade ge 0 and evts.grade le 26 and x ge 0 and y ge 0 '+colfilt+')')
    if check eq 0 then stop,'MKIMGS: Error with execute command.'
    push,ii,thisii
    tt=where(gti.start ge t1[t] and gti.stop lt t2[t])
    for i=0,n_elements(tt)-1 do $
          if gti[tt[i]].stop gt t1[t] and gti[tt[i]].start lt t1[t] then $
              exptime+=gti[tt[i]].stop-t1[t] $
          else if gti[tt[i]].start lt t2[t] and gti[tt[i]].stop gt t2[t] then $
              exptime+=t2[t]-gti[tt[i]].start else $
          if t1[t] lt gti[tt[i]].start and t2[t] gt gti[tt[i]].stop then $
              exptime+=gti[tt[i]].stop-gti[tt[i]].start
endfor
deadc=sxpar(head2,'DEADC')
exptime*=deadc
sxaddpar,head2,'EXPOSURE',exptime
sxaddpar,head2,'LIVETIME',exptime

i=1
check=1
while i lt 100 and check ne 0 do begin
    if str(sxpar(header,'TTYPE'+str(i))) eq xtt and $
          str(sxpar(header,'TTYPE'+str(i+1))) eq ytt then check=0
    i++
endwhile
xx=str(i-1)
yy=str(i)

x0=fix(sxpar(header,'TLMIN'+xx))
y0=fix(sxpar(header,'TLMIN'+yy))
im=intarr(fix(sxpar(header,'TLMAX'+xx))-x0+1,$
      fix(sxpar(header,'TLMAX'+yy))-y0+1)
for i=0,n_elements(ii)-1 do im[x[ii[i]]-x0,y[ii[i]]-y0]++

sxaddpar,head2,'BITPIX',32
sxaddpar,head2,'NAXIS',2
sxaddpar,head2,'NAXIS1',fix(sxpar(header,'TLMAX'+xx))-fix(sxpar(header,'TLMIN'+xx))+1
sxaddpar,head2,'NAXIS2',fix(sxpar(header,'TLMAX'+yy))-fix(sxpar(header,'TLMIN'+yy))+1
sxaddpar,head2,'CTYPE1',sxpar(header,'TCTYP'+xx)
sxaddpar,head2,'CTYPE2',sxpar(header,'TCTYP'+yy)
sxaddpar,head2,'CRPIX1',sxpar(header,'TCRPX'+xx)
sxaddpar,head2,'CRPIX2',sxpar(header,'TCRPX'+yy)
sxaddpar,head2,'CRVAL1',sxpar(header,'TCRVL'+xx)
sxaddpar,head2,'CRVAL2',sxpar(header,'TCRVL'+yy)
sxaddpar,head2,'CDELT1',sxpar(header,'TCDLT'+xx)
sxaddpar,head2,'CDELT2',sxpar(header,'TCDLT'+yy)
if keyword_set(fexp) then begin
    sxaddpar,head2,'EXPOSURE',sxpar(header,'EXPOSURE')*fexp
    sxaddpar,head2,'LIVETIME',sxpar(header,'EXPOSURE')*fexp
endif

fits_write,outname,im,head2

end
