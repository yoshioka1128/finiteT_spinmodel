subroutine surfaceKave_xgb(Klcr,Ku,percent,nr,nunit,nunitx,nunity,nunitz,irpt,irvx)
implicit none
integer(4) :: nr,nunit,nunitx,nunity,nunitz,irvx
integer(4) :: iy,iz,irpt(nr,nunitx,nunity,nunitz)
real(8) :: Klcr(9,nr*nunit),Ku,percent

! averaged
do iz=1,nunitz
   do iy=1,nunity
         Klcr(1,irpt(4,irvx,iy,iz))=Ku*(1.0d0-percent)
         Klcr(1,irpt(5,irvx,iy,iz))=Ku*(1.0d0-percent)
         Klcr(1,irpt(2,irvx,iy,iz))=Ku*(1.0d0-percent)
         Klcr(1,irpt(7,irvx,iy,iz))=Ku*(1.0d0-percent)
         Klcr(1,irpt(1,irvx,iy,iz))=Ku*(1.0d0-percent)
         Klcr(1,irpt(8,irvx,iy,iz))=Ku*(1.0d0-percent)
         Klcr(1,irpt(3,irvx,iy,iz))=Ku*(1.0d0-percent)
         Klcr(1,irpt(6,irvx,iy,iz))=Ku*(1.0d0-percent)
   end do
end do

end subroutine surfaceKave_xgb
