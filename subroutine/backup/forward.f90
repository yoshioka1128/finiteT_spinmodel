!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine forward(v,v0,nr,nunit,Gt,f,ik,abc,nrtot) ! forward difference
implicit none
integer(4) :: nr,nunit,ir,idim,ik,nrtot,ir0
real(8) :: v(nr*nunit,3),v0(nr*nunit,3),f(nr*nunit,3,4),Gt,abc

!do ir=1,nr*nunit
do ir=1,nrtot
   do idim=1,3
      v(ir,idim)=v0(ir,idim)+abc*f(ir,idim,ik)*Gt
   end do
end do

end subroutine forward
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine normalize(v,nr,nunit,nrtot) ! normalization of v
implicit none
integer(4) :: nr,nunit,ir,idim,nrtot,ir0
real(8) :: v(nr*nunit,3),abs

do ir=1,nrtot
   abs=0.0d0
   do idim=1,3
      abs=abs+v(ir,idim)**2
   end do
   abs=dsqrt(abs)
   if(abs.eq.0.0d0) then
      do idim=1,3
         v(ir,idim)=0.0d0 ! for surface
      end do
   else
      do idim=1,3
         v(ir,idim)=v(ir,idim)/abs
      end do
   end if
end do

end subroutine normalize
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


