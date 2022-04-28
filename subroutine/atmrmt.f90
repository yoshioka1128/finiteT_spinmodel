subroutine atmrmt(dir,ratom,iatommax,r0,rmt)
integer(4),parameter :: iatmax=50
integer(4) :: iatommax,i,j,iat,mult
real(8) :: r0(iatommax),rmt(iatommax)
character(32) :: ach,bch,cch,dch,ech,fch
character(len=255) :: dir,ratom(iatmax)

open(120,file=''//trim(dir)//'.struct')
read(120,"(A5)") ach
read(120,"(A27,I3)") bch,iatommax
do i=1,2
   read(120,*) 
end do
do iat=1,iatommax
   read(120,"(A4,I4)") ach,i
   read(120,"(A15,I2)") ach,mult
   if(mult.ne.1) then
      do i=1,mult-1
         read(120,*)
      end do
   end if
   read(120,"(A2,A13,I5,A5,f10.9,A5,f10.4)") ratom(iat),ach,npt,bch,r0(iat),cch,rmt(iat)
   do i=1,3
      read(120,*)
   end do
end do

return
end subroutine atmrmt
