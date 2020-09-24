pro trimspecresp,fullbase,fulloutbase,nreg,nsrc,elow,ehigh,cmprmf=cmprmf

lo=round((elow-1.6)/0.04)
hi=round((ehigh-1.6)/0.04)
nchan=hi-lo

ii=strsplit(fullbase,'/')
if n_elements(ii) eq 1 then begin
    base=fullbase
    dir=''
endif else begin
    base=strmid(fullbase,ii[n_elements(ii)-1],strlen(fullbase)-ii[n_elements(ii)-1])
    dir=strmid(fullbase,0,ii[n_elements(ii)-1])
endelse
ii=strsplit(fulloutbase,'/')
if n_elements(ii) eq 1 then begin
    outbase=fulloutbase
    outdir=''
endif else begin
    outbase=strmid(fulloutbase,ii[n_elements(ii)-1], $
          strlen(fulloutbase)-ii[n_elements(ii)-1])
    outdir=strmid(fulloutbase,0,ii[n_elements(ii)-1])
endelse

for r=0,nreg-1 do begin
    specstr=mrdfits(fullbase+str(r+1)+'.pha',1,spech,/silent)
    tspecstr=specstr[lo:hi-1]
    sxaddpar,spech,'NAXIS2',nchan
    sxaddpar,spech,'DETCHANS',nchan
    sxaddpar,spech,'TLMAX1',nchan-1
    sxaddpar,spech,'BACKFILE','bgd'+outbase+str(r+1)+'.pha'
    mwrfits,tspecstr,fulloutbase+str(r+1)+'.pha',spech,/create,/silent

    specstr=mrdfits(dir+'bgd'+base+str(r+1)+'.pha',1,spech,/silent)
    tspecstr=specstr[lo:hi-1]
    sxaddpar,spech,'NAXIS2',nchan
    sxaddpar,spech,'DETCHANS',nchan
    sxaddpar,spech,'TLMAX1',nchan-1
    mwrfits,tspecstr,outdir+'bgd'+outbase+str(r+1)+'.pha',spech,$
          /create,/silent

    rmfstr=mrdfits(fullbase+str(r+1)+'.rmf',2,rh,/silent)
    r2str=mrdfits(fullbase+str(r+1)+'.rmf',1,r2h,/silent)
    tr2str=r2str[lo:hi-1]
    sxaddpar,r2h,'NAXIS2',nchan
    sxaddpar,r2h,'DETCHANS',nchan
    sxaddpar,r2h,'TLMAX1',nchan-1

;    mwrfits,tr2str,fulloutbase+str(r+1)+'.rmf',r2h,/silent,/create
    mwrfits,tr2str,'temp.rmf',r2h,/silent,/create
    trmfstr={trmfstr,ENERG_LO:0.,ENERG_HI:0.,N_GRP:0,F_CHAN:0,$
          N_CHAN:0,MATRIX:fltarr(nchan)}
    trmfstr=replicate(trmfstr,nchan)
    trmfstr.energ_lo=rmfstr[lo:hi-1].energ_lo
    trmfstr.energ_hi=rmfstr[lo:hi-1].energ_hi
    trmfstr.n_grp=rmfstr[lo:hi-1].n_grp
    trmfstr.f_chan=rmfstr[lo:hi-1].f_chan
;    trmfstr.n_chan=rmfstr[lo:hi-1].n_chan
    trmfstr.n_chan=replicate(nchan,n_elements(trmfstr.n_chan))
    trmfstr.matrix=rmfstr[lo:hi-1].matrix[lo:hi-1]
    sxaddpar,rh,'NAXIS2',nchan
;    sxaddpar,rh,'NAXIS1',8206
    sxaddpar,rh,'NAXIS1',34
    sxaddpar,rh,'TFORM6',str(nchan)+'E'
    sxaddpar,rh,'DETCHANS',nchan
    sxaddpar,rh,'TLMAX4',nchan-1
;    mwrfits,trmfstr,fulloutbase+str(r+1)+'.rmf',rh,/silent
    mwrfits,trmfstr,'temp.rmf',rh,/silent

    if keyword_set(cmprmf) then spawn,'cmprmf temp.rmf '+ $
          fulloutbase+str(r+1)+'.rmf 1e-6 clobber=yes' $
        else spawn,'cp -f temp.rmf '+fulloutbase+str(r+1)+'.rmf'
    spawn,'rm -f temp.rmf'

    for s=0,nsrc-1 do begin
        if nsrc eq 1 then regn='' else regn='_'+str(r+1)
        arfstr=mrdfits(fullbase+str(s+1)+regn+'.arf',1,harf,/silent)
        tarfstr=arfstr[lo:hi-1]
        sxaddpar,harf,'NAXIS2',nchan
        mwrfits,tarfstr,fulloutbase+str(s+1)+regn+'.arf',harf,$
              /silent,/create
    endfor

endfor

end
