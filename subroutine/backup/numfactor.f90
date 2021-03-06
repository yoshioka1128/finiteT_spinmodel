subroutine numfactor(cclm)
real(8) :: cclm(6,-6:6)
real(8),parameter :: pi=dacos(-1.0d0)
integer(4) :: l,m

! set clm table
cclm=0.0d0
!cclm(0,0)=dsqrt(1.0d0/pi)/2.0d0
cclm(2,0)=dsqrt(5.0d0/pi)/4.0d0
cclm(2,1)=dsqrt(15.0d0/pi)/2.0d0
cclm(2,2)=dsqrt(15.0d0/pi)/4.0d0
cclm(4,0)=dsqrt(1.0d0/pi)*3.0d0/16.0d0
cclm(4,1)=dsqrt(10.0d0/pi)*3.0d0/8.0d0
cclm(4,2)=dsqrt(5.0d0/pi)*3.0d0/8.0d0
cclm(4,3)=dsqrt(70.0d0/pi)*3.0d0/8.0d0
cclm(4,4)=dsqrt(35.0d0/pi)*3.0d0/16.0d0
cclm(6,0)=dsqrt(13.0d0/pi)/32.0d0
cclm(6,1)=dsqrt(273.0d0/pi)/16.0d0
cclm(6,2)=dsqrt(2730d0/pi)/64.0d0
cclm(6,3)=dsqrt(2730d0/pi)/32.0d0
cclm(6,4)=dsqrt(13.0d0/7.0d0/pi)*21.0d0/32.0d0
cclm(6,5)=dsqrt(2002.0d0/pi)*3.0d0/32.0d0
cclm(6,6)=dsqrt(26.0d0/231.0d0/pi)*231.0d0/64.0d0

do l=2,6,2
   do m=0,6
      cclm(l,-m)=cclm(l,m)
   end do
end do

return
end subroutine numfactor
