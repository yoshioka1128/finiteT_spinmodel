subroutine moment(imax,L1,S1,zmoment,zmomentS,JM)
implicit none
real(8) :: J,Jz,Jzst,Ost,Ost2,abc,theta,JM(100,2),J1,J2,M1,M2,f,L1,S1
integer(4) :: l,m,i,ii,jj,imax,i1,i2
complex(kind(0d0)) :: zmoment(imax,imax,3),imag=(0.0d0,1.0d0),zmomentS(imax,imax,3)

zmoment=0.0d0
zmomentS=0.0d0
do i2=1,imax
   J2=JM(i2,1)
   M2=JM(i2,2)
   do i1=1,imax
      J1=JM(i1,1)
      M1=JM(i1,2)
      if(J1.eq.J2) then
         if(M1.eq.M2) then
            zmoment(i1,i1,3)=JM(i1,2)
         else if(M2.eq.M1+1.0d0) then
            Ost=dsqrt((J1-M1)*(J1+M1+1.0d0))
            zmoment(i2,i1,1)=Ost/2.0d0
            zmoment(i2,i1,2)=-imag*Ost/2.0d0
         end if
      else if(J1+1.0d0.eq.J2) then
         f=dsqrt((J1+L1+S1+2.0d0)*(-J1+S1+L1)*(J1+S1-L1+1.0d0)*(J1+L1-S1+1.0d0)*(J1+M1+1.0d0)*(J1-M1+1.0d0)/&
              (4.0d0*(J1+1)**2*(2.0d0*J1+1.0d0)*(2.0d0*J1+3.0d0)))
         if(M1.eq.M2) then
            zmomentS(i2,i1,3)=f
         else if(M1+1.0d0.eq.M2) then
            zmomentS(i2,i1,1)=-0.5d0*f*dsqrt((J1+M1+2.0d0)/(J1-M1+1))
            zmomentS(i2,i1,2)=0.5d0*imag*f*dsqrt((J1+M1+2.0d0)/(J1-M1+1))
         else if(M1-1.0d0.eq.M2) then
            zmomentS(i2,i1,1)=0.5d0*f*dsqrt((J1-M1+2.0d0)/(J1+M1+1))
            zmomentS(i2,i1,2)=0.5d0*imag*f*dsqrt((J1-M1+2.0d0)/(J1+M1+1))
         end if
      end if
   end do
end do

do i=1,3
   do ii=1,imax
      do jj=ii,imax
         zmoment(ii,jj,i)=conjg(zmoment(jj,ii,i))
         zmomentS(ii,jj,i)=conjg(zmomentS(jj,ii,i))
      end do
   end do
end do

return
end subroutine moment


