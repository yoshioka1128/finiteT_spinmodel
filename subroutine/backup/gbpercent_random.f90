subroutine grainboundary(nunit,nunitx,nunity,nunitz,nr,nfe,vr0,vfe0,irpt,ifept,rr,rfe,igb,gbpercent,icountb)
integer(4) :: nunit,nunitx,nunity,nunitz,nr,nfe,ix,iy,iz,i,idim,ife,ixx
integer(4) :: irpt(nr,nunitx,nunity,nunitz),ifept(nfe,nunitx,nunity,nunitz)
real(8) :: vr0(nr*nunit,3),vfe0(nfe*nunit,3),rfe(nfe,3),rr(nr,3),delta(1000),xx,eta=0.2d0
real(8) :: pi=dacos(-1.0d0)
real(8) :: ran(nunity*nunitz*5*nfe),gbpercent
integer(4) :: icount,seedsize,ic,i1,i2,i3,i4,irv0,irvx0,irvy0,igb(2),icountgb,icountrd
integer,allocatable :: seed(:)

gbpercent=gbpercent/100.0d0
! remove Fe randomly
call system_clock(count=ic)
call random_seed(size=seedsize)
allocate(seed(seedsize))
call random_seed(get=seed)
seed=3984598
call random_seed(put=seed)
call random_number(ran)

icountrd=0
! grain boundary correction
icount=0
do iz=1,nunitz
   do iy=1,nunity
      do ix=igb(1),igb(2)
         do i=1,nr
            do idim=1,3
               vr0(irpt(i,ix,iy,iz),idim)=0.0d0
            end do
         end do
         do i=1,nfe
            icount=icount+1
            if(ran(icount).lt.1.0d0-gbpercent) then
               icountrd=icountrd+1
               do idim=1,3
                  vfe0(ifept(i,ix,iy,iz),idim)=0.0d0
               end do
            end if
         end do
      end do
   end do
end do

icountgb=dble(igb(2)-igb(1)+1)*dble(nunity*nunitz*nfe)-icountrd
!write(6,*) icountrd,icountgb
write(6,*) "gb ",100.0d0*dble(icountgb)/dble(igb(2)-igb(1)+1)/dble(nunity*nunitz*nfe),"% Fe exist"

end subroutine grainboundary
