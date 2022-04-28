subroutine Hexch(theta,phi,imax,zmatrixr,zHr,zspin,Hex)
implicit none
integer(4) :: i,j,imax
real(8) :: Hex,theta,phi
complex(kind(0d0)) :: zmatrixr(imax,imax),zHr(imax,imax),zspin(imax,imax,3)

do j=1,imax
   do i=1,imax
      zmatrixr(i,j)=zHr(i,j)&
           +(zspin(i,j,1)*dcos(phi)*dsin(theta)&
           +zspin(i,j,2)*dsin(phi)*dsin(theta)&
           +zspin(i,j,3)*dcos(theta))*2.0d0*Hex
   end do
end do

return
end subroutine Hexch
