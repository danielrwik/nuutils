pro clickn,xarr,yarr,ncl

undefine,xarr
undefine,yarr
for i=0,ncl-1 do begin
  chan,0
  cursor,curx,cury
;  print,curx,cury
  wait,0.5
  push,xarr,curx
  push,yarr,cury
endfor


end
