subroutine defectR_x(Klcr,Ku,ran,defectrate,irpt,nr,nunit,nunitx,nunity,nunitz,icount,nunitxbd)
integer(4) :: nunit,nunitx,nunity,nunitz,nr
integer(4) :: irpt(nr,nunitx,nunity,nunitz),ix,iy,iz,nunitxbd(2)
integer,allocatable :: seed(:)
integer(4) :: i1,i2,i3,i4,icount,seedsize,ic
real(8) :: Klcr(9,nr*nunit),Ku,ran(nunity*nunitz*4),defectrate
real(8),parameter :: mev=1.60218d-21 ! (Merg/mev) 
real(8),parameter :: kb=8.6173303d-2*mev ! (Merg/K)

! random
call system_clock(count=ic)
call random_seed(size=seedsize)
allocate(seed(seedsize))
call random_seed(get=seed)
seed=492983
call random_seed(put=seed)
call random_number(ran) ! random number 0.0 to 1.0

! reset Klcr
!do iz=1,nunitz
!   do iy=1,nunity
!      ix=idef
!      do ir=1,nr
!         do ik=1,9
!            Klcr(ik,irpt(ir,ix,iy,iz))=0.0d0
!         end do
!      end do
!   end do
!end do

if(defectrate.le.0.5d0) then
   do iz=1,nunitz
      do iy=1,nunity
         ix=nunitxbd(2)
         i1=iy+(iz-1)*nunity
         i2=iy+(iz-1)*nunity+nunitz*nunity
         i3=iy+(iz-1)*nunity+2*nunitz*nunity
         i4=iy+(iz-1)*nunity+3*nunitz*nunity
         if(ran(i1)/2.0d0.lt.defectrate) then
            Klcr(1,irpt(4,ix,iy,iz))=-Ku
         end if
         if(ran(i2)/2.0d0.lt.defectrate) then
            Klcr(1,irpt(5,ix,iy,iz))=-Ku
         end if
         if(ran(i3)/2.0d0.lt.defectrate) then
            Klcr(1,irpt(2,ix,iy,iz))=-Ku
         end if
         if(ran(i4)/2.0d0.lt.defectrate) then
            Klcr(1,irpt(7,ix,iy,iz))=-Ku
         end if
      end do
   end do
else
   do iz=1,nunitz
      do iy=1,nunity
         ix=nunitxbd(2)
         Klcr(1,irpt(4,ix,iy,iz))=-Ku
         Klcr(1,irpt(5,ix,iy,iz))=-Ku
         Klcr(1,irpt(2,ix,iy,iz))=-Ku
         Klcr(1,irpt(7,ix,iy,iz))=-Ku
      end do
   end do

   do iz=1,nunitz
      do iy=1,nunity
         ix=nunitxbd(2)
         i1=iy+(iz-1)*nunity
         i2=iy+(iz-1)*nunity+nunitz*nunity
         i3=iy+(iz-1)*nunity+2*nunitz*nunity
         i4=iy+(iz-1)*nunity+3*nunitz*nunity
         if(ran(i1)/2.0d0.lt.defectrate-0.5d0) then
            Klcr(1,irpt(1,ix,iy,iz))=-Ku
         end if
         if(ran(i2)/2.0d0.lt.defectrate-0.5d0) then
            Klcr(1,irpt(3,ix,iy,iz))=-Ku
         end if
         if(ran(i3)/2.0d0.lt.defectrate-0.5d0) then
            Klcr(1,irpt(6,ix,iy,iz))=-Ku
         end if
         if(ran(i4)/2.0d0.lt.defectrate-0.5d0) then
            Klcr(1,irpt(8,ix,iy,iz))=-Ku
         end if
      end do
   end do
end if

open(93,file="R_defect_K1_in_gb_phase.txt")
icount=0
do ix=1,nunitx
   do iy=1,nunity
      do iz=1,nunitz
         do i=1,nr
            if(Klcr(1,irpt(i,ix,iy,iz)).lt.0.0d0) then
               write(93,*) i,ix,iy,iz,Klcr(1,irpt(i,ix,iy,iz))/kb
               icount=icount+1
            end if
         end do
      end do
   end do
end do

end subroutine defectR_x
