subroutine vector_to_angle(vector,theta,phi,abs)
real(8),parameter :: pi=dacos(-1.0d0)
real(8) :: vector(3),theta,phi,abs
integer(4) :: i

abs=dsqrt(vector(1)**2+vector(2)**2+vector(3)**2)
if(abs.eq.0.0d0) then
   theta=0.0d0
   phi=0.0d0
else
   theta=dacos(vector(3)/abs)
   if(theta.eq.0) then
      phi=0.0d0
   else
      phi=dacos(vector(1)/dsqrt(vector(1)**2+vector(2)**2))
      if(vector(2).le.0.0d0) then
         phi=2.0d0*pi-phi
      end if
   end if
end if

end subroutine vector_to_angle
