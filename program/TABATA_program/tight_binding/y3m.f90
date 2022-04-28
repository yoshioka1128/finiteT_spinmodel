function y3m(m,theta,phi)
  implicit none
  integer :: i,j,k
  integer :: m                          !magnetic quantum number
  double precision :: theta,phi         !angular
  double complex :: y3m
  complex*16,parameter :: imag=cmplx(0d0,1d0)
  if(m==3) then
     y3m=-sqrt(5d0/16d0)*(sin(theta)**3d0)*(cos(3d0*phi)-imag*sin(3d0*phi))
  else if(m==2) then
     y3m=sqrt(15d0/8d0)*(sin(theta)**2d0)*cos(theta)*(cos(2d0*phi)-imag*sin(2d0*phi))
  else if(m==1) then
     y3m=-sqrt(3d0/16d0)*sin(theta)*(5d0*cos(theta)*cos(theta)-1d0)*(cos(phi)-imag*sin(phi))
  else if(m==0) then
     y3m=sqrt(1d0/4d0)*(5d0*(cos(theta)**3d0)-3d0*cos(theta))
  else if(m==-1) then
     y3m=sqrt(3d0/16d0)*sin(theta)*(5d0*cos(theta)*cos(theta)-1d0)*(cos(phi)+imag*sin(phi))
  else if(m==-2) then
     y3m=sqrt(15d0/8d0)*(sin(theta)**2d0)*cos(theta)*(cos(2d0*phi)+imag*sin(2d0*phi))
  else if(m==-3) then
     y3m=sqrt(5d0/16d0)*(sin(theta)**3d0)*(cos(3d0*phi)+imag*sin(3d0*phi))
  else
     y3m=0
  end if
  
  return 
end function y3m
