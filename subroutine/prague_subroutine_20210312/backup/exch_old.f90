subroutine exch(imax,J,zmag,theta,phi)
implicit none
real(8) :: J,Jz,Jzst,Ost,Ost2,abc,theta,phi
integer(4) :: imax
integer(4) :: l,m,i,ii
complex(kind(0d0)) :: zmag(imax,imax),imag=(0.0d0,1.0d0)

zmag=0.0d0
Jz=-J-1.0d0
do ii=1,imax
   Jz=Jz+1.0d0
   zmag(ii,ii)=Jz*dcos(theta)
end do

Jzst=-J-1.0d0
do ii=1,imax-1
   Jzst=Jzst+1.0d0
   Jz=Jzst
   Ost=dsqrt((J-Jz)*(J+Jz+1.0d0))
   zmag(ii+1,ii)=0.5d0*(dcos(phi)-imag*dsin(phi))*Ost*dsin(theta)
!   zmag(ii+1,ii)=dsqrt(2.0d0)/4.0d0*(1.0d0-imag)*Ost*dsin(theta)
end do

return
end subroutine exch
