!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine position_RFe12_allsite(irpt,irpt0,ifept,ifept0,nunit,nunitx,nunity,nunitz,nr,nfe,rfe,rr,rfept,rrpt,ax,cx,&
     ipmx,ipmy,ipmz,iunit,irsub0,irunit0,ifesub0,ifeunit0)
implicit none
integer(4) :: nunit,nr,nfe,nunitx,nunity,nunitz,irpbc(nr,0:nunitx+1,0:nunity+1,0:nunitz+1),ifepbc(nfe,0:nunitx+1,0:nunity+1,0:nunitz+1)
integer(4) :: i,i2,ir,ir2,ife,ife2,ix,iy,iz,ix2,iy2,iz2
integer(4) :: ipmx(nunitx,-1:1),ipmy(nunity,-1:1),ipmz(nunitz,-1:1)
integer(4) :: irpt(nr,nunitx,nunity,nunitz),ifept(nfe,nunitx,nunity,nunitz)
integer(4) :: irpt0(nr,nunitx,nunity,nunitz),ifept0(nfe,nunitx,nunity,nunitz)
integer(4) :: ipt1,ipt2,idim,idim2
real(8) :: rfe(nfe,3),rr(nr,3),ax,cx,abc,rrpt(nr,0:nunitx+1,0:nunity+1,0:nunitz+1,3),rfept(nfe,0:nunitx+1,0:nunity+1,0:nunitz+1,3)
real(8) :: dis,dr(3),dfe(3),drfe(3),dfer(3)
integer(4) :: iunit(nunitx,nunity,nunitz)
integer(4) :: irunit0(nr*nunit,3),ifeunit0(nfe*nunit,3),irsub0(nr*nunit),ifesub0(nfe*nunit)

! number of unitcell
i=0
do iz=1,nunitz
   do iy=1,nunity
      do ix=1,nunitx
         i=i+1
         iunit(ix,iy,iz)=i
      end do
   end do
end do

! number of site
do iz=1,nunitz
   do iy=1,nunity
      do ix=1,nunitx
         do i=1,nr
            irpt(i,ix,iy,iz)=i+(ix-1)*nr+(iy-1)*nr*nunitx+(iz-1)*nr*nunitx*nunity
            irpt0(i,ix,iy,iz)=i+(ix-1)*nr+(iy-1)*nr*nunitx+(iz-1)*nr*nunitx*nunity
         end do
         do i=1,nfe
            ifept(i,ix,iy,iz)=i+(ix-1)*nfe+(iy-1)*nfe*nunitx+(iz-1)*nfe*nunitx*nunity
            ifept0(i,ix,iy,iz)=i+(ix-1)*nfe+(iy-1)*nfe*nunitx+(iz-1)*nfe*nunitx*nunity
         end do
      end do
   end do
end do

! number of unitcell
i=0
do iz=1,nunitz
   do iy=1,nunity
      do ix=1,nunitx
         i=i+1
         iunit(ix,iy,iz)=i
      end do
   end do
end do

! shift unitcell
do ix=1,nunitx
   ipmx(ix,-1)=ix-1
   ipmx(ix,0)=ix
   ipmx(ix,1)=ix+1
end do
do iy=1,nunity
   ipmy(iy,-1)=iy-1
   ipmy(iy,0)=iy
   ipmy(iy,1)=iy+1
end do
do iz=1,nunitz
   ipmz(iz,-1)=iz-1
   ipmz(iz,0)=iz
   ipmz(iz,1)=iz+1
end do
ipmx(nunitx,1)=1
ipmx(1,-1)=nunitx
ipmy(nunity,1)=1
ipmy(1,-1)=nunity
ipmz(nunitz,1)=1
ipmz(1,-1)=nunitz

rr=0.0d0
rfe=0.0d0
! RE
 ! RE1
rr(1,1)=0.00000000d0*ax
rr(1,2)=0.00000000d0*ax
rr(1,3)=0.00000000d0*cx
 ! RE2
rr(2,1)=0.50000000d0*ax
rr(2,2)=0.50000000d0*ax
rr(2,3)=0.50000000d0*cx
! Fe
 ! Fe1
rfe(1,1)=0.25000000d0*ax
rfe(1,2)=0.25000000d0*ax
rfe(1,3)=0.25000000d0*cx
 ! Fe2
rfe(2,1)=0.75000000d0*ax
rfe(2,2)=0.25000000d0*ax
rfe(2,3)=0.25000000d0*cx
 ! Fe3
rfe(3,1)=0.75000000d0*ax
rfe(3,2)=0.75000000d0*ax
rfe(3,3)=0.25000000d0*cx
 ! Fe4
rfe(4,1)=0.25000000d0*ax
rfe(4,2)=0.75000000d0*ax
rfe(4,3)=0.25000000d0*cx
 ! Fe5
rfe(5,1)=0.75000000d0*ax
rfe(5,2)=0.75000000d0*ax
rfe(5,3)=0.75000000d0*cx
 ! Fe6
rfe(6,1)=0.25000000d0*ax
rfe(6,2)=0.75000000d0*ax
rfe(6,3)=0.75000000d0*cx
 ! Fe7
rfe(7,1)=0.25000000d0*ax
rfe(7,2)=0.25000000d0*ax
rfe(7,3)=0.75000000d0*cx
 ! Fe8
rfe(8,1)=0.75000000d0*ax
rfe(8,2)=0.25000000d0*ax
rfe(8,3)=0.75000000d0*cx
 ! Fe9
rfe(9,1)=0.35900000d0*ax
rfe(9,2)=0.00000000d0*ax
rfe(9,3)=0.00000000d0*cx
 ! Fe10
rfe(10,1)=0.64100000d0*ax
rfe(10,2)=0.00000000d0*ax
rfe(10,3)=0.00000000d0*cx
 ! Fe11
rfe(11,1)=0.00000000d0*ax
rfe(11,2)=0.35900000d0*ax
rfe(11,3)=0.00000000d0*cx
 ! Fe12
rfe(12,1)=0.00000000d0*ax
rfe(12,2)=0.64100000d0*ax
rfe(12,3)=0.00000000d0*cx
 ! Fe13
rfe(13,1)=0.85900000d0*ax
rfe(13,2)=0.50000000d0*ax
rfe(13,3)=0.50000000d0*cx
 ! Fe14
rfe(14,1)=0.14100000d0*ax
rfe(14,2)=0.50000000d0*ax
rfe(14,3)=0.50000000d0*cx
 ! Fe15
rfe(15,1)=0.50000000d0*ax
rfe(15,2)=0.85900000d0*ax
rfe(15,3)=0.50000000d0*cx
 ! Fe16
rfe(16,1)=0.50000000d0*ax
rfe(16,2)=0.14100000d0*ax
rfe(16,3)=0.50000000d0*cx
 ! Fe17
rfe(17,1)=0.27000000d0*ax
rfe(17,2)=0.50000000d0*ax
rfe(17,3)=0.00000000d0*cx
 ! Fe18
rfe(18,1)=0.73000000d0*ax
rfe(18,2)=0.50000000d0*ax
rfe(18,3)=0.00000000d0*cx
 ! Fe19
rfe(19,1)=0.50000000d0*ax
rfe(19,2)=0.27000000d0*ax
rfe(19,3)=0.00000000d0*cx
 ! Fe20
rfe(20,1)=0.50000000d0*ax
rfe(20,2)=0.73000000d0*ax
rfe(20,3)=0.00000000d0*cx
 ! Fe21
rfe(21,1)=0.77000000d0*ax
rfe(21,2)=0.00000000d0*ax
rfe(21,3)=0.50000000d0*cx
 ! Fe22
rfe(22,1)=0.23000000d0*ax
rfe(22,2)=0.00000000d0*ax
rfe(22,3)=0.50000000d0*cx
 ! Fe23
rfe(23,1)=0.00000000d0*ax
rfe(23,2)=0.77000000d0*ax
rfe(23,3)=0.50000000d0*cx
 ! Fe24
rfe(24,1)=0.00000000d0*ax
rfe(24,2)=0.23000000d0*ax
rfe(24,3)=0.50000000d0*cx

! position of all sites
rrpt=0.0d0
rfept=0.0d0
do iz=0,nunitz+1
   do iy=0,nunity+1
      do ix=0,nunitx+1
         do i=1,nr
            rrpt(i,ix,iy,iz,1)=rr(i,1)+dble(ix-1)*ax
            rrpt(i,ix,iy,iz,2)=rr(i,2)+dble(iy-1)*ax
            rrpt(i,ix,iy,iz,3)=rr(i,3)+dble(iz-1)*cx
         end do
         do i=1,nfe
            rfept(i,ix,iy,iz,1)=rfe(i,1)+dble(ix-1)*ax
            rfept(i,ix,iy,iz,2)=rfe(i,2)+dble(iy-1)*ax
            rfept(i,ix,iy,iz,3)=rfe(i,3)+dble(iz-1)*cx
         end do
      end do
   end do
end do

! site number
do iz=1,nunitz
   do iy=1,nunity
      do ix=1,nunitx
         do i=1,nr
            ir=irpt0(i,ix,iy,iz)
            irunit0(ir,1)=ix
            irunit0(ir,2)=iy
            irunit0(ir,3)=iz
            irsub0(ir)=i
         end do
         do i=1,nfe
            ife=ifept0(i,ix,iy,iz)
            ifeunit0(ife,1)=ix
            ifeunit0(ife,2)=iy
            ifeunit0(ife,3)=iz
            ifesub0(ife)=i
         end do
      end do
   end do
end do

end subroutine position_RFe12_allsite
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
