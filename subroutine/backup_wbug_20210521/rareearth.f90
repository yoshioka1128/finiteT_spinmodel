subroutine rareearth(ratom,n4f,nJex,imax,JM,S1,L1,J1)
implicit none
integer(4),parameter :: mm=100
integer(4) :: n4f,nJex,icount,i,jj,imax
real(8) :: S1,L1,J1,JM(mm,2),Jz,J
character(15) ratom

if(ratom.eq."Ce") n4f=1
if(ratom.eq."Pr") n4f=2
if(ratom.eq."Nd") n4f=3
if(ratom.eq."Pm") n4f=4
if(ratom.eq."Sm") n4f=5
if(ratom.eq."Eu") n4f=6
if(ratom.eq."Gd") n4f=7
if(ratom.eq."Tb") n4f=8
if(ratom.eq."Dy") n4f=9
if(ratom.eq."Ho") n4f=10
if(ratom.eq."Er") n4f=11
if(ratom.eq."Tm") n4f=12
if(ratom.eq."Yb") n4f=13
if(nJex.ne.0.and.n4f.ge.8) then
   write(6,*) "error: heavy rare-earth"
   stop
end if
if(n4f.lt.7) then
   S1=dble(n4f)/2.0d0
   L1=dble((2*3-n4f+1)*n4f/2)
   J1=abs(L1-S1)
else
   S1=dble(14-n4f)/2.0d0
   L1=(2*3-(n4f-7)+1)*(n4f-7)/2
   J1=abs(L1+S1)
end if
icount=0
JM=0.0d0
do jj=0,nJex
   J=J1+dble(jj)
   do i=0,int(J*2.0d0)
      icount=icount+1
      Jz=-J+dble(i)
      JM(icount,1)=J
      JM(icount,2)=Jz
   end do
end do
imax=icount

if(imax.gt.mm) then
   write(6,*) "error: imax > mm"
   stop
end if

end subroutine rareearth
