!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Intexch(ax,cx,range,nunit,nunitx,nunity,nunitz,numint0,ipmx,ipmy,ipmz,iIntr,numintr,&
     magR,irpt,irptpbc,irpt0,irunit,irsub,nr,nrtot,rr,magFe,ifept,ifeptpbc,ifept0,ifeunit,ifesub,nfe,nfetot,rfe)
implicit none
integer(4) :: numint0,nunit,nr,nfe,iIntr(nr*nunit,numint0),nunitx,nunity,nunitz,irpt(nr,nunitx,nunity,nunitz),ifept(nfe,nunitx,nunity,nunitz)
integer(4) :: i,i2,ir,ife,ipmx(nunitx,-1:1),ipmy(nunity,-1:1),ipmz(nunitz,-1:1)
integer(4) :: irptpbc(nr,0:nunitx+1,0:nunity+1,0:nunitz+1),ifeptpbc(nfe,0:nunitx+1,0:nunity+1,0:nunitz+1)
integer(4) :: ix,iy,iz,ix2,iy2,iz2,ipt1,ipt2,icount,icount0,idx,idy,idz,numintr(nr*nunit),idim,ir0,ife0,nrtot,nfetot
real(8) :: ax,cx,abc,rr(nr,3),rfe(nfe,3),range,magFe(nfe*nunit),magR(nr*nunit),dis
integer(4) :: irunit(nr*nunit,3),ifeunit(nfe*nunit,3),irsub(nr*nunit),ifesub(nfe*nunit),irpt0(nr*nunit),ifept0(nfe*nunit)

iIntr=0.0d0
numintr=0
do ir0=1,nrtot
   ir=irpt0(ir0)
   iz=irunit(ir,3)
   iy=irunit(ir,2)
   ix=irunit(ir,1)
   i=irsub(ir)
   if(magR(ir).eq.0.0d0) cycle
   
   icount=0
   
   do idz=-1,1 ! ifeunit-irunit
      do idy=-1,1
         do idx=-1,1
            do i2=1,nfe
!               ife=ifeptpbc(i2,ix+idx,iy+idy,iz+idz)
               ife=ifeptpbc(i2,ipmx(ix,idx),ipmy(iy,idy),ipmz(iz,idz))
               dis=0.0d0
               dis=dis+(rr(i,1)-(rfe(i2,1)+dble(idx)*ax))**2
               dis=dis+(rr(i,2)-(rfe(i2,2)+dble(idy)*ax))**2
               dis=dis+(rr(i,3)-(rfe(i2,3)+dble(idz)*cx))**2
               dis=dsqrt(dis)
               
               if((dis.lt.range).and.(dis.ne.0.0d0).and.ife.ne.0) then
!               if((dis.lt.range).and.(dis.ne.0.0d0).and.magFe(ife).ne.0) then
                  icount=icount+1
                  iIntr(ir,icount)=ife
                  numintr(ir)=icount
                  if(icount.gt.numint0) then
                     write(6,*) "the number of itneracting pair is over limit"
                     stop
                  end if
               end if
               
            end do
         end do
      end do
      
   end do
end do

end subroutine Intexch
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
