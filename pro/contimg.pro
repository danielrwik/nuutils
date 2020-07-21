function contimg,im,lev

xim=n_elements(im[*,0])
yim=n_elements(im[0,*])
cont=intarr(xim,yim)
c1=intarr(xim,yim)
c2=intarr(xim,yim)
c3=intarr(xim,yim)
ii=where(im ge lev)
if ii[0] ne -1 then begin
    c1[ii]=1
    ii=where(shift(im,1,0) ge lev)
    if ii[0] ne -1 then c2[ii]=1
    ii=where(c1+c2 eq 1)
    if ii[0] ne -1 then cont[ii]=1
    ii=where(shift(im,0,1) ge lev)
    if ii[0] ne -1 then c3[ii]=1
    ii=where(c1+c3 eq 1)
    if ii[0] ne -1 then cont[ii]=1
endif
cont[0,*]=0 & cont[-1,*]=0 & cont[*,0]=0 & cont[*,-1]=0

ii=where(cont gt 0.5)
ii2d=array_indices(cont,ii)

return,cont

end
