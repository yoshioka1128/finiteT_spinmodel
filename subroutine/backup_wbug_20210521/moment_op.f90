subroutine moment_op(imax,L1,S1,zspin,zorbit,JM,Vst,Ust)
implicit none
integer(4),parameter :: mm=100
real(8) :: J,Jz,Jzst,Ost,Ost2,abc,theta,JM(mm,2),J1,J2,M1,M2,f,L1,S1,rdcoffS,rdcoffL,Vst(imax,imax,-1:1)
real(8) :: Ust(imax,imax,6,-6:6)
integer(4) :: l,m,i,ii,jj,imax,i1,i2,isign
complex(kind(0d0)) :: im=(0.0d0,1.0d0),zspin(imax,imax,3),zorbit(imax,imax,3)

! reduced matrix element of 1st tensor
rdcoffS=dsqrt(S1*(S1+1.0d0)*(2.0d0*S1+1.0d0))
rdcoffL=dsqrt(L1*(L1+1.0d0)*(2.0d0*L1+1.0d0))

zspin=0.0d0
do i2=1,imax
   M2=JM(i2,2)
   do i1=1,imax
      M1=JM(i1,2)
      if(M1.eq.M2) then
         zspin(i1,i2,3)=Vst(i1,i2,0)*rdcoffS
      else if(M1-1.0d0.eq.M2) then
         zspin(i1,i2,1)=-Vst(i1,i2,1)*rdcoffS/dsqrt(2.0d0)
         zspin(i1,i2,2)=im*Vst(i1,i2,1)*rdcoffS/dsqrt(2.0d0)
      else if(M1+1.0d0.eq.M2) then
         zspin(i1,i2,1)=Vst(i1,i2,-1)*rdcoffS/dsqrt(2.0d0)
         zspin(i1,i2,2)=im*Vst(i1,i2,-1)*rdcoffS/dsqrt(2.0d0)
      end if
   end do
end do

zorbit=0.0d0
do i2=1,imax
   M2=JM(i2,2)
   do i1=1,imax
      M1=JM(i1,2)
      if(M1.eq.M2) then
         zorbit(i1,i2,3)=Ust(i1,i2,1,0)*rdcoffL
      else if(M1-1.0d0.eq.M2) then
         zorbit(i1,i2,1)=-Ust(i1,i2,1,1)*rdcoffL/dsqrt(2.0d0)
         zorbit(i1,i2,2)=im*Ust(i1,i2,1,1)*rdcoffL/dsqrt(2.0d0)
      else if(M1+1.0d0.eq.M2) then
         zorbit(i1,i2,1)=Ust(i1,i2,1,-1)*rdcoffL/dsqrt(2.0d0)
         zorbit(i1,i2,2)=im*Ust(i1,i2,1,-1)*rdcoffL/dsqrt(2.0d0)
      end if
   end do
end do

return
end subroutine moment_op

