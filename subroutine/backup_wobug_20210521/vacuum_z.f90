subroutine vacuum(nunitx,nunity,nunitz,nunit,nr,nfe,vr0,vfe0,irpt,ifept)
integer(4) :: nunitx,nunity,nunitz,nr,nfe,ix,iy,i,nunit
integer(4) :: irpt(nr,nunitx,nunity,nunitz),ifept(nfe,nunitx,nunity,nunitz)
real(8) :: vr0(nr*nunit,3),vfe0(nfe*nunit,3)

! surface correction
ifept0=1
do iy=1,nunity ! vacuum
   do ix=1,nunitx
      do i=1,nr
         do idim=1,3
            vr0(irpt(i,ix,iy,1),idim)=0.0d0
            vr0(irpt(i,ix,iy,nunitz),idim)=0.0d0
         end do
      end do
      do i=1,nfe
         do idim=1,3
            vfe0(ifept(i,ix,iy,1),idim)=0.0d0
            vfe0(ifept(i,ix,iy,nunitz),idim)=0.0d0
         end do
      end do
   end do
end do

end subroutine vacuum
