subroutine rareearth_SmL3(ratom,n4f,nJex,imax,JM,S1,L1,J1)
implicit none
integer(4),parameter :: mm=100
integer(4) :: n4f,nJex,icount,i,jj,imax
real(8) :: S1,L1,J1,JM(mm,2),Jz,J
character(15) ratom

if(ratom.eq."Sm") then
   n4f=5
else
   write(6,*) "error R is not Sm"
end if

if(nJex.ne.0.and.n4f.ge.8) then
   write(6,*) "error: heavy rare-earth"
   stop
end if
S1=dble(n4f)/2.0d0
L1=dble(3)
J1=abs(L1-S1) ! L1=3, S1=2.5, J1=0.5, excited J1=1.5, 2.5, 3.5, 4.5, 5.5

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

write(6,*) "JM list"
do i=1,imax
   write(6,*) i,JM(i,1),JM(i,2)
end do

end subroutine rareearth_SmL3
