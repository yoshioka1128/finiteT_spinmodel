subroutine defectR_clusterx4(Klcr,Ku,ran,defectrate,irpt,nr,nunit,nunitx,nunity,nunitz,icount,nunitxbd)
integer(4) :: nunit,nunitx,nunity,nunitz,nr,nunitxbd(2)
integer(4) :: irpt(nr,nunitx,nunity,nunitz),ix,iy,iz
integer,allocatable :: seed(:)
integer(4) :: i1,i2,i3,i4,icount,seedsize,ic
real(8) :: Klcr(9,nr*nunit),Ku,ran(nunity*nunitz*4),defectrate
real(8),parameter :: mev=1.60218d-21 ! (Merg/mev) 
real(8),parameter :: kb=8.6173303d-2*mev ! (Merg/K)

do iz=nunitz/2,nunitz/2+1
   do iy=nunity/2,nunity/2+1
      ix=nunitxbd(2)
      do i=1,nr
         Klcr(1,irpt(i,ix,iy,iz))=-Ku
         do ik=2,9
            Klcr(ik,irpt(i,ix,iy,iz))=0.0d0
         end do
      end do
   end do
end do

open(93,file="R_defect_K1_in_gb_phase.txt")
icount=0
do iz=1,nunitz
   do iy=1,nunity
      do ix=1,nunitx
         do i=1,nr
            if(Klcr(1,irpt(i,ix,iy,iz)).lt.0.0d0) then
               write(93,*) i,ix,iy,iz,Klcr(1,irpt(i,ix,iy,iz))/kb
               icount=icount+1
            end if
         end do
      end do
   end do
end do

end subroutine defectR_clusterx4
