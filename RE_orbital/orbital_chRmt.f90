program order
real(8) :: a,b,dx1,dx2,r0=0.00001d0,c,rmtin,rmtout,x1(0:781),x2(0:781),dr1(781+1),dr2(781+1),abc,r4fin(0:781),r4fout(781)
integer(4) :: i1,i2,imax=781,l
character(32) rmtoutch,rmtinch,atom,file

write(6,*) "atom name ?"
read(5,*) atom

write(6,*) "file name ?"
read(5,*) file

write(6,*) "input Rmt ?"
read(5,*) rmtin
write(6,*) "output Rmt ?"
read(5,*) rmtout

write(rmtoutch,'(f5.1)') rmtout
write(rmtinch,'(f5.1)') rmtin
!if(rmt.ne.3.2d0) then
!   write(6,*) "rmt.ne.3.2 error"
!   stop
!end if

dx1=(1.0d0/dble(imax-1))*log(rmtin/r0)
dx2=(1.0d0/dble(imax-1))*log(rmtout/r0)
open(18,file=''//trim(adjustl(atom))//''//trim(adjustl(rmtoutch))//'_R4f.dat')
open(20,file=''//trim(adjustl(atom))//''//trim(adjustl(rmtinch))//'_R4f.dat')
open(19,file='./source/'//trim(adjustl(file))//'')
do i1=1,imax
   read(19,*) r4fin(i1)
end do
r4fin(0)=0
x1(0)=0

do i2=1,imax ! output pints
   x2(i2)=r0*exp(dx2*dble(i2-1))
   dr2(i2)=r0*(dexp(dx2*dble(i2-1))-dexp(dx2*dble(i2-2)))
   do i1=1,imax ! input points
      x1(i1)=r0*exp(dx1*dble(i1-1))
      dr1(i1)=r0*(dexp(dx1*dble(i1-1))-dexp(dx1*dble(i1-2)))
      if(x1(i1).gt.x2(i2)) then
         r4fout(i2)=(r4fin(i1)-r4fin(i1-1))/dr1(i1)*(x2(i2)-x1(i1-1))+r4fin(i1-1)
         write(18,*) x2(i2),r4fout(i2)
         exit
      end if
   end do
end do

dr2(1)=r0
dr2(782)=0.0d0


do l=2,8,2
   a=0.0d0
   do i=1,781
      a=a+x2(i)**(l-2)*r4fout(i)*(dr2(i)+dr2(i+1))/2.0d0
   end do
   write(6,*) "l=",l-2,"<r^l>",a
end do

end program order
