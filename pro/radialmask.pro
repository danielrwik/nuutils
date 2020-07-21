function radialmask,imfile,center,angle,rot=rot

fits_read,imfile,im,header
xim=im
yim=im
x=findgen(n_elements(im[*,0]))+1
y=findgen(n_elements(im[0,*]))+1
mask=intarr(n_elements(x),n_elements(y))
for i=0,n_elements(y)-1 do xim[*,i]=x
for i=0,n_elements(x)-1 do yim[i,*]=y

xim=xim-center[0]
yim=yim-center[1]
if size(rot,/type) ne 0 then begin
    x0=xim
    xim=x0*cos(!pi*rot/180.)-yim*sin(!pi*rot/180.)
    yim=x0*sin(!pi*rot/180.)+yim*cos(!pi*rot/180.)
endif
angim=atan(yim,xim)*180./!pi
ii=where(angim lt 0.)
angim[ii]=angim[ii]+360.

nang=fix(360./angle)
for i=0,nang-1 do begin
    ii=where(angim ge i*angle and angim lt (i+1)*angle)
    mask[ii]=i+1
endfor

return,mask

end
