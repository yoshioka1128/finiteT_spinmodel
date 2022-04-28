!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine reset_position(irpt,irpt0,irpbc,ifept,ifept0,ifepbc,nunit,nunitx,nunity,nunitz,nr,nfe,&
     nrtot,nfetot,irsub,irunit,ifesub,ifeunit,ipmx,ipmy,ipmz,pbcx,pbcy,pbcz)
implicit none
integer(4) :: nunit,nr,nfe,nunitx,nunity,nunitz,nrtot,nfetot,icount,idx,idy,idz
integer(4) :: irpbc(nr,0:nunitx+1,0:nunity+1,0:nunitz+1),ifepbc(nfe,0:nunitx+1,0:nunity+1,0:nunitz+1)
integer(4) :: irpt(nr,nunitx,nunity,nunitz),ifept(nfe,nunitx,nunity,nunitz)
integer(4) :: irpt0(nr,nunitx,nunity,nunitz),ifept0(nfe,nunitx,nunity,nunitz)
integer(4) :: i,i2,ir,ir2,ife,ife2,ix,iy,iz,ix2,iy2,iz2
integer(4) :: ipt1,ipt2,idim,idim2
integer(4) :: ipmx(nunitx,-1:1),ipmy(nunity,-1:1),ipmz(nunitz,-1:1)
real(8) :: abc,dis,dr(3),dfe(3),drfe(3),dfer(3)
integer(4) :: irunit(nr*nunit,3),ifeunit(nfe*nunit,3),irsub(nr*nunit),ifesub(nfe*nunit)
character(10) :: pbcx,pbcy,pbcz ! "on" or "off"

icount=0
irunit=0
irsub=0
do iz=1,nunitz
   do iy=1,nunity
      do ix=1,nunitx
         do i=1,nr
            ir=irpt(i,ix,iy,iz)
            if(ir.ne.0) then
               icount=icount+1
               irunit(icount,1)=ix
               irunit(icount,2)=iy
               irunit(icount,3)=iz
               irsub(icount)=i
               irpt(i,ix,iy,iz)=icount
            end if
         end do
      end do
   end do
end do
nrtot=icount

icount=0
ifeunit=0
ifesub=0
do iz=1,nunitz
   do iy=1,nunity
      do ix=1,nunitx
         do i=1,nfe
            ife=ifept(i,ix,iy,iz)
            if(ife.ne.0) then
               icount=icount+1
               ifeunit(icount,1)=ix
               ifeunit(icount,2)=iy
               ifeunit(icount,3)=iz
               ifesub(icount)=i
               ifept(i,ix,iy,iz)=icount
            end if
         end do

      end do
   end do
end do
nfetot=icount

! set pbc
do iz=1,nunitz
do iy=1,nunity
do ix=1,nunitx
   do idz=-1,1
      do idy=-1,1
         do idx=-1,1
            do i=1,nr
               irpbc(i,ix+idx,iy+idy,iz+idz)=irpt(i,ipmx(ix,idx),ipmy(iy,idy),ipmz(iz,idz))
            end do
            do i=1,nfe
               ifepbc(i,ix+idx,iy+idy,iz+idz)=ifept(i,ipmx(ix,idx),ipmy(iy,idy),ipmz(iz,idz))
            end do
         end do
      end do
   end do
end do
end do
end do

if(pbcx.eq."off") then
   do iz=0,nunitz+1
      do iy=0,nunity+1
         do i=1,nr
            irpbc(i,0,iy,iz)=0
            irpbc(i,nunitx+1,iy,iz)=0
         end do
         do i=1,nfe
            ifepbc(i,0,iy,iz)=0
            ifepbc(i,nunitx+1,iy,iz)=0
         end do
      end do
   end do
end if
if(pbcy.eq."off") then
   do iz=0,nunitz+1
      do ix=0,nunitx+1
         do i=1,nr
            irpbc(i,ix,0,iz)=0
            irpbc(i,ix,nunity+1,iz)=0
         end do
         do i=1,nfe
            ifepbc(i,ix,0,iz)=0
            ifepbc(i,ix,nunity+1,iz)=0
         end do
      end do
   end do
end if
if(pbcz.eq."off") then
   do iy=0,nunity+1
      do ix=0,nunitx+1
         do i=1,nr
            irpbc(i,ix,iy,0)=0
            irpbc(i,ix,iy,nunitz+1)=0
         end do
         do i=1,nfe
            ifepbc(i,ix,iy,0)=0
            ifepbc(i,ix,iy,nunitz+1)=0
         end do
      end do
   end do
end if
!stop

end subroutine reset_position
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
