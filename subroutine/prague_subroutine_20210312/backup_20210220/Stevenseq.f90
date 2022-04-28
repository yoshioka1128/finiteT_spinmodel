subroutine Stevenseq(imax,J,ZOst)
implicit none
real(8) :: J,Jz,Jzst,Ost,Ost2,abc
integer(4) :: imax
integer(4) :: l,m,i,ii
complex(kind(0d0)) :: ZOst(imax,imax,6,-6:6),imag=(0.0d0,1.0d0)

ZOst=0.0d0
! Stevens eq operator
l=2
m=0
!write(6,*) "l,m",l,m
Jz=-J-1.0d0
do ii=1,imax
   Jz=Jz+1.0d0
   ZOst(ii+m,ii,l,m)=3.0d0*Jz**2-J*(J+1)
end do

l=2
m=2
!write(6,*) "l,m",l,m
Jzst=-J-1.0d0
do ii=1,imax-m
   Jzst=Jzst+1.0d0
   Jz=Jzst
   Ost=1.0d0
   do i=1,m
      Ost=Ost*dsqrt((J-Jz)*(J+Jz+1.0d0))
      Jz=Jz+1.0d0
   end do
   ZOst(ii+m,ii,l,m)=-0.5d0*Ost*imag
end do

l=4
m=0
Jz=-J-1.0d0
do ii=1,imax
   Jz=Jz+1.0d0
   ZOst(ii+m,ii,l,m)=35.0d0*Jz**4-30.0d0*J*(J+1.0d0)*Jz**2+25.0d0*Jz**2-6.0d0*J*(J+1.0d0)+3.0d0*J**2*(J+1.0d0)**2
end do

l=4
m=2
Jzst=-J-1.0d0
do ii=1,imax-m
   Jzst=Jzst+1.0d0
   Jz=Jzst
   Ost=1.0d0
   Ost2=(7.0d0*Jz**2-J*(J+1.0d0)-5.0d0)
   do i=1,m
      Ost=Ost*dsqrt((J-Jz)*(J+Jz+1.0d0))
      Ost2=Ost2*dsqrt((J-Jz)*(J+Jz+1.0d0)) 
      Jz=Jz+1.0d0
   end do
   ZOst(ii+m,ii,l,m)=-imag*0.25*(Ost*(7.0d0*Jz**2-J*(J+1.0d0)-5.0d0)+Ost2)
end do

l=4
m=4
Jzst=-J-1.0d0
do ii=1,imax-m
   Jzst=Jzst+1.0d0
   Jz=Jzst
   Ost=1.0d0
   do i=1,m
      Ost=Ost*dsqrt((J-Jz)*(J+Jz+1))
      Jz=Jz+1.0d0
   end do
   ZOst(ii+m,ii,l,m)=Ost*0.5d0
end do

l=6
m=0
Jz=-J-1.0d0
do ii=1,imax
   Jz=Jz+1.0d0
   ZOst(ii+m,ii,l,m)=231.0d0*Jz**6-315.0d0*J*(J+1)*Jz**4+735.0d0*Jz**4+105.0d0*J**2*(J+1.0d0)**2*Jz**2-525.0d0*J*(J+1.0d0)*Jz**2+294.0d0*Jz**2-5.0d0*J**3*(J+1.0d0)**3+40.0d0*J**2*(J+1.0d0)**2-60.0d0*J*(J+1.0d0)
end do

l=6
m=2
Jzst=-J-1.0d0
do ii=1,imax-m
   Jzst=Jzst+1.0d0
   Jz=Jzst
   Ost=1.0d0
   Ost2=(33.0d0*Jz**4-18.0d0*Jz**2*J*(J+1.0d0)-123.0d0*Jz**2+J**2*(J+1.0d0)**2+10.0d0*J*(J+1.0d0)+102.0d0)
   do i=1,m
      Ost=Ost*dsqrt((J-Jz)*(J+Jz+1.0d0))
      Ost2=Ost2*dsqrt((J-Jz)*(J+Jz+1.0d0)) 
      Jz=Jz+1.0d0
   end do
   ZOst(ii+m,ii,l,m)=-imag*0.25*(Ost*(33.0d0*Jz**4-18.0d0*Jz**2*J*(J+1.0d0)-123.0d0*Jz**2+J**2*(J+1.0d0)**2+10.0d0*J*(J+1.0d0)+102.0d0)+Ost2)
end do


l=6
m=4
Jzst=-J-1.0d0
do ii=1,imax-m
   Jzst=Jzst+1.0d0
   Jz=Jzst
   Ost=1.0d0
   Ost2=(11.0d0*Jz**2-J*(J+1.0d0)-38.0d0)
   do i=1,m
      Ost=Ost*dsqrt((J-Jz)*(J+Jz+1))
      Ost2=Ost2*dsqrt((J-Jz)*(J+Jz+1))
      Jz=Jz+1.0d0
   end do
   abc=(11.0d0*Jz**2-J*(J+1.0d0)-38.0d0)*Ost+Ost2
   ZOst(ii+m,ii,l,m)=abc/4.0d0
end do


l=6
m=6
Jzst=-J-1.0d0
do ii=1,imax-m
   Jzst=Jzst+1.0d0
   Jz=Jzst
   Ost=1.0d0
   do i=1,m
      Ost=Ost*dsqrt((J-Jz)*(J+Jz+1.0d0))
      Jz=Jz+1.0d0
   end do
   ZOst(ii+m,ii,l,m)=-0.5d0*Ost*imag
end do

return
end subroutine Stevenseq
