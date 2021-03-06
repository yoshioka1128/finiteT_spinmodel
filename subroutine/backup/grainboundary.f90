subroutine grainboundary_edge_SmFe12(nunit,nunitx,nunity,nunitz,nr,nfe,irpt,ifept,igb,gbrate,icountgb,icountfe,ran,rfept,gbv,ax,pbcx,pbcy,pbcz)
integer(4) :: nunit,nunitx,nunity,nunitz,nr,nfe,ix,iy,iz,i,ik,idim,ife,ifegb,ixx
integer(4) :: irpt(nr,nunitx,nunity,nunitz),ifept(nfe,nunitx,nunity,nunitz)
real(8) :: delta(1000),xx,eta=0.2d0
real(8) :: rfept(nfe,0:nunitx+1,0:nunity+1,0:nunitz+1,3),rfe,rfegb,drfe
real(8) :: pi=dacos(-1.0d0)
real(8) :: ran(nunit*nfe),gbrate,gbv(6)
integer(4) :: icountfe,seedsize,ic,i1,i2,i3,i4,irv0,irvx0,irvy0,igb(2,3),icountgb,icountrd
real(8),parameter :: muB=9.274d-24 ! (Merg/kOe) or (kerg/Oe)
integer,allocatable :: seed(:)
character(10) :: pbcx,pbcy,pbcz

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
write(41,*) "x-direction",igb(1,1),igb(2,1)
write(41,*) "y-direction",igb(1,2),igb(2,2)
write(41,*) "z-direction",igb(1,3),igb(2,3)
write(41,*) "i,ix,iy,iz,mag"

icountrd=0
! grain boundary correction
icountfe=0
do iz=1,nunitz
   do iy=1,nunity
      do ix=1,nunitx

         if(pbcx.eq."off".and.ix.eq.nunitx) then
            do i=1,nr
               if(i.ne.1) then
                  irpt(i,ix,iy,iz)=0
               end if
            end do
            do i=1,nfe
               if(i.ne.11.and.i.ne.12.and.i.ne.23.and.i.ne.24) then
                  ifept(i,ix,iy,iz)=0
               end if
            end do
         end if

         if(pbcy.eq."off".and.iy.eq.nunity) then
            do i=1,nr
               if(i.ne.1) then
                  irpt(i,ix,iy,iz)=0
               end if
            end do
            do i=1,nfe
               if(i.ne.9.and.i.ne.10.and.i.ne.21.and.i.ne.22) then
                  ifept(i,ix,iy,iz)=0
               end if
            end do
         end if

         if(pbcz.eq."off".and.iz.eq.nunitz) then
            do i=1,nr
               if(i.ne.1) then
                  irpt(i,ix,iy,iz)=0
               end if
            end do
            do i=1,nfe
               if(i.ne.9.and.i.ne.10.and.i.ne.11.and.i.ne.12.and.i.ne.17.and.i.ne.18.and.i.ne.19.and.i.ne.20) then
                  ifept(i,ix,iy,iz)=0
               end if
            end do
         end if

      end do
   end do
end do
icountgb=icountfe-icountrd
!write(6,*) icountgb,icountfe-icountrd

end subroutine grainboundary_edge_SmFe12


subroutine grainboundary(nunit,nunitx,nunity,nunitz,nr,nfe,irpt,ifept,igb,gbrate,icountgb,icountfe,ran)
integer(4) :: nunit,nunitx,nunity,nunitz,nr,nfe,ix,iy,iz,i,ik,idim,ife,ixx
integer(4) :: irpt(nr,nunitx,nunity,nunitz),ifept(nfe,nunitx,nunity,nunitz)
real(8) :: delta(1000),xx,eta=0.2d0
real(8) :: pi=dacos(-1.0d0)
real(8) :: ran(nunit*nfe),gbrate
integer(4) :: icountfe,seedsize,ic,i1,i2,i3,i4,irv0,irvx0,irvy0,igb(2,3),icountgb,icountrd
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
write(41,*) "x-direction",igb(1,1),igb(2,1)
write(41,*) "y-direction",igb(1,2),igb(2,2)
write(41,*) "z-direction",igb(1,3),igb(2,3)
write(41,*) "i,ix,iy,iz,mag"

icountrd=0
! grain boundary correction
icountfe=0
do iz=1,nunitz
   do iy=1,nunity
      do ix=1,nunitx
         if((ix.ge.igb(1,1).and.ix.le.igb(2,1)).or.(iy.ge.igb(1,2).and.iy.le.igb(2,2)).or.(iz.ge.igb(1,3).and.iz.le.igb(2,3))) then
            do i=1,nr
               irpt(i,ix,iy,iz)=0
            end do
            do i=1,nfe
               icountfe=icountfe+1
               if(ran(icountfe).lt.1.0d0-gbrate) then
                  icountrd=icountrd+1
                  ifept(i,ix,iy,iz)=0
               end if
            end do
         end if
      end do
   end do
end do
icountgb=icountfe-icountrd
!write(6,*) icountgb,icountfe-icountrd

end subroutine grainboundary




subroutine grainboundary_valley(nunit,nunitx,nunity,nunitz,nr,nfe,irpt,ifept,igb,gbrate,icountgb,icountfe,ran,rfept,gbv,ax)
integer(4) :: nunit,nunitx,nunity,nunitz,nr,nfe,ix,iy,iz,i,ik,idim,ife,ifegb,ixx
integer(4) :: irpt(nr,nunitx,nunity,nunitz),ifept(nfe,nunitx,nunity,nunitz)
real(8) :: delta(1000),xx,eta=0.2d0,ax
real(8) :: rfept(nfe,0:nunitx+1,0:nunity+1,0:nunitz+1,3),rfe,rfegb,drfe
real(8) :: pi=dacos(-1.0d0)
real(8) :: ran(nunit*nfe),gbrate,gbv(6)
integer(4) :: icountfe,seedsize,ic,i1,i2,i3,i4,irv0,irvx0,irvy0,igb(2,3),icountgb,icountrd
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
write(41,*) "x-direction",igb(1,1),igb(2,1)
write(41,*) "y-direction",igb(1,2),igb(2,2)
write(41,*) "z-direction",igb(1,3),igb(2,3)
write(41,*) "i,ix,iy,iz,mag"

icountrd=0
! grain boundary correction
icountfe=0
do iz=1,nunitz
   do iy=1,nunity
      do ix=1,nunitx

         if((ix.ge.igb(1,1).and.ix.le.igb(2,1)).or.(iy.ge.igb(1,2).and.iy.le.igb(2,2)).or.(iz.ge.igb(1,3).and.iz.le.igb(2,3))) then
            do i=1,nr
               irpt(i,ix,iy,iz)=0
            end do
            do i=1,nfe
               drfe=rfept(i,ix,iy,iz,1)-rfept(51,igb(1,1),iy,iz,1)
               icountfe=icountfe+1
!               if(iz.eq.nunitz/2) write(6,*) drfe/ax
!               if(rfept(i,ix,iy,iz).ge.igb)
!               write(6,*) ifept(i,ix,iy,iz),Klcfe(1,ifept(i,ix,iy,iz))
               if(drfe.le.0.5d0*ax) then
                  if(ran(icountfe).lt.1.0d0-gbv(1)) then
                     icountrd=icountrd+1
                     ifept(i,ix,iy,iz)=0
                  end if
               else if(drfe.le.1.0d0*ax) then
                  if(ran(icountfe).lt.1.0d0-gbv(2)) then
                     icountrd=icountrd+1
                     ifept(i,ix,iy,iz)=0
                  end if
               else if(drfe.le.1.5d0*ax) then
                  if(ran(icountfe).lt.1.0d0-gbv(3)) then
                     icountrd=icountrd+1
                     ifept(i,ix,iy,iz)=0
                  end if
               else if(drfe.le.2.0d0*ax) then
                  if(ran(icountfe).lt.1.0d0-gbv(4)) then
                     icountrd=icountrd+1
                     ifept(i,ix,iy,iz)=0
                  end if
               else if(drfe.le.2.5d0*ax) then
                  if(ran(icountfe).lt.1.0d0-gbv(5)) then
                     icountrd=icountrd+1
                     ifept(i,ix,iy,iz)=0
                  end if
               else if(drfe.le.3.0d0*ax) then
                  if(ran(icountfe).lt.1.0d0-gbv(6)) then
                     icountrd=icountrd+1
                     ifept(i,ix,iy,iz)=0
                  end if
               end if
            end do
         end if

      end do
   end do
end do
icountgb=icountfe-icountrd
!write(6,*) icountgb,icountfe-icountrd

end subroutine grainboundary_valley



subroutine grainboundary_edge_Nd2Fe14B(nunit,nunitx,nunity,nunitz,nr,nfe,irpt,ifept,igb,gbrate,icountgb,icountfe,ran,rfept,gbv,ax,pbcx,pbcy,pbcz)
integer(4) :: nunit,nunitx,nunity,nunitz,nr,nfe,ix,iy,iz,i,ik,idim,ife,ifegb,ixx
integer(4) :: irpt(nr,nunitx,nunity,nunitz),ifept(nfe,nunitx,nunity,nunitz)
real(8) :: delta(1000),xx,eta=0.2d0
real(8) :: rfept(nfe,0:nunitx+1,0:nunity+1,0:nunitz+1,3),rfe,rfegb,drfe
real(8) :: pi=dacos(-1.0d0)
real(8) :: ran(nunit*nfe),gbrate,gbv(6)
integer(4) :: icountfe,seedsize,ic,i1,i2,i3,i4,irv0,irvx0,irvy0,igb(2,3),icountgb,icountrd
real(8),parameter :: muB=9.274d-24 ! (Merg/kOe) or (kerg/Oe)
integer,allocatable :: seed(:)
character(10) :: pbcx,pbcy,pbcz

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
write(41,*) "x-direction",igb(1,1),igb(2,1)
write(41,*) "y-direction",igb(1,2),igb(2,2)
write(41,*) "z-direction",igb(1,3),igb(2,3)
write(41,*) "i,ix,iy,iz,mag"

icountrd=0
! grain boundary correction
icountfe=0
do iz=1,nunitz
   do iy=1,nunity
      do ix=1,nunitx

         if(pbcx.eq."off".and.ix.eq.nunitx) then
            do i=1,nr
               irpt(i,ix,iy,iz)=0
            end do
            do i=1,nfe
               if(i.ne.51.and.i.ne.52.and.i.ne.53.and.i.ne.54) then
                  ifept(i,ix,iy,iz)=0
               end if
            end do
         end if

         if(pbcy.eq."off".and.iy.eq.nunity) then
            do i=1,nr
               irpt(i,ix,iy,iz)=0
            end do
            do i=1,nfe
               if(i.ne.51.and.i.ne.52.and.i.ne.55.and.i.ne.56) then
                  ifept(i,ix,iy,iz)=0
               end if
            end do
         end if

         if(pbcz.eq."off".and.iz.eq.1) then
            do i=1,nr
               if(i.ne.1.and.i.ne.2.and.i.ne.3.and.i.ne.4) then
                  irpt(i,ix,iy,iz)=0
               end if
            end do
            do i=1,nfe
               if(i.ne.53.and.i.ne.56) then
                  ifept(i,ix,iy,iz)=0
               end if
            end do
         end if

      end do
   end do
end do
icountgb=icountfe-icountrd
!write(6,*) icountgb,icountfe-icountrd

end subroutine grainboundary_edge_Nd2Fe14B
