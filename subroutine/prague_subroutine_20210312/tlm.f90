subroutine tlm(coffmag,theta,phi)
implicit none
real(8) :: theta,phi,coffmag(6,-6:6)

coffmag=0.0d0
! m=0
coffmag(2,0)=(3.0d0*dcos(theta)**2-1.0d0)
coffmag(4,0)=(35.0d0*dcos(theta)**4-30.0d0*dcos(theta)**2+3.0d0) ! error in second term as 30*cos(theeta)
coffmag(6,0)=(231.0d0*dcos(theta)**6-315.0d0*dcos(theta)**4+105.d0*dcos(theta)**2-5.0d0)
! m=2
coffmag(2,2)=dcos(2.0d0*phi)*dsin(theta)**2
coffmag(4,2)=dcos(2.0d0*phi)*dsin(theta)**2*(7.0d0*dcos(theta)**2-1.0d0)
coffmag(6,2)=dcos(2.0d0*phi)*dsin(theta)**2*(33.0d0*dcos(theta)**4-18.0d0*dcos(theta)**2+1.0d0)
! m=4
coffmag(4,4)=dcos(4.0d0*phi)*dsin(theta)**4
coffmag(6,4)=dcos(4.0d0*phi)*dsin(theta)**4*(11.0d0*dcos(theta)**2-1.0d0)

return
end subroutine tlm
