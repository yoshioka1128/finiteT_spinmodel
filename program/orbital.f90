program order
real(8) :: a,b,dx,r0=0.00001d0,c,rmt,x(0:781),dr(781+1),abc,r4f(781)
integer(4) :: i,imax=781,l
character(32) rmtch,atom,file

write(6,*) "atom name ?"
read(5,*) atom

write(6,*) "file name ?"
read(5,*) file

write(6,*) "Rmt ?"
read(5,*) rmt
write(rmtch,'(f5.1)') rmt
if(rmt.ne.3.2d0) then
   write(6,*) "rmt.ne.3.2 error"
   stop
end if

dx=(1.0d0/dble(imax-1))*log(rmt/r0)
open(18,file=''//trim(adjustl(atom))//''//trim(adjustl(rmtch))//'_R4f.dat')
open(19,file='./source/'//trim(adjustl(file))//'')
do i=1,imax
   x(i)=r0*exp(dx*dble(i-1))
   dr(i)=r0*(dexp(dx*dble(i-1))-dexp(dx*dble(i-2)))
   read(19,*) r4f(i)
   write(18,*) x(i),r4f(i)
end do
dr(1)=r0
dr(782)=0.0d0


do l=2,8,2
   a=0.0d0
   do i=1,781
      a=a+x(i)**(l-2)*r4f(i)*(dr(i)+dr(i+1))/2.0d0
   end do
   write(6,*) "l=",l-2,"<r^l>",a
end do

end program order
