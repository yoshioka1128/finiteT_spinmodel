subroutine grainboundary_x(nunit,nunitx,nunity,nunitz,nr,nfe,magR,magFe,Klcfe,Klcr,irpt,ifept,igb,gbrate,icountgb,ran,rdgbmag)
integer(4) :: nunit,nunitx,nunity,nunitz,nr,nfe,ix,iy,iz,i,ik,idim,ife,ixx
integer(4) :: irpt(nr,nunitx,nunity,nunitz),ifept(nfe,nunitx,nunity,nunitz)
real(8) :: delta(1000),xx,eta=0.2d0
real(8) :: pi=dacos(-1.0d0),Klcfe(9,nfe*nunit),Klcr(9,nr*nunit)
real(8) :: ran(nunity*nunitz*5*nfe),gbrate,magR(nr*nunit),magFe(nfe*nunit),rdgbmag
integer(4) :: icount,seedsize,ic,i1,i2,i3,i4,irv0,irvx0,irvy0,igb(2),icountgb,icountrd
real(8),parameter :: muB=9.274d-24 ! (Merg/kOe) or (kerg/Oe)
integer,allocatable :: seed(:)

! remove Fe randomly
call system_clock(count=ic)
call random_seed(size=seedsize)
allocate(seed(seedsize))
call random_seed(get=seed)
seed=3984598
call random_seed(put=seed)
call random_number(ran)
open(41,file="Fe_moments_in_gb_phase.txt")
write(41,*) "moments in gb phase"
write(41,*) "i,ix,iy,iz,mag"
icountrd=0
! grain boundary correction
icount=0
do iz=1,nunitz
   do iy=1,nunity
      do ix=igb(1),igb(2)
         do i=1,nr
            magR(irpt(i,ix,iy,iz))=0.0d0
            do ik=1,9
               Klcr(ik,irpt(i,ix,iy,iz))=0.0d0
            end do
         end do
         do i=1,nfe
            ife=ifept(i,ix,iy,iz)
            icount=icount+1
            Klcfe(1,ife)=Klcfe(1,ife)*gbrate
            if(ran(icount).lt.1.0d0-gbrate) then
               icountrd=icountrd+1
               magFe(ife)=0.0d0
            else
               magFe(ife)=magFe(ife)*rdgbmag
               write(41,*) i,ix,iy,iz,magFe(ife)/muB
            end if
         end do
      end do
   end do
end do

icountgb=dble(igb(2)-igb(1)+1)*dble(nunity*nunitz*nfe)-icountrd

end subroutine grainboundary_x
