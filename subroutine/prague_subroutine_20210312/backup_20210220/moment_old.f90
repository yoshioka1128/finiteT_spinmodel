subroutine moment(imax,J,zmoment)
implicit none
real(8) :: J,Jz,Jzst,Ost,Ost2,abc,theta
integer(4) :: imax
integer(4) :: l,m,i,ii
complex(kind(0d0)) :: zmoment(imax,imax,3),imag=(0.0d0,1.0d0)

zmoment=0.0d0
Jz=-J-1.0d0
do ii=1,imax
   Jz=Jz+1.0d0
   zmoment(ii,ii,3)=Jz
end do

Jzst=-J-1.0d0
do ii=1,imax-1
   Jzst=Jzst+1.0d0
   Jz=Jzst
   Ost=dsqrt((J-Jz)*(J+Jz+1.0d0))
   zmoment(ii+1,ii,1)=Ost/2.0d0
   zmoment(ii+1,ii,2)=-imag*Ost/2.0d0
end do
do i=1,2
   do ii=1,imax-1
      zmoment(ii,ii+1,i)=conjg(zmoment(ii+1,ii,i))
   end do
end do

return
end subroutine moment


