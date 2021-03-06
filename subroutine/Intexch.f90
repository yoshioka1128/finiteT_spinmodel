!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Intexch(ax,cx,range,nunit,nunitx,nunity,nunitz,numint0,ipmx,ipmy,ipmz,iIntr,numintr,&
     irpbc,irunit,irsub,nr,nrtot,rrpt,ifepbc,ifeunit,ifesub,nfe,nfetot,rfept)
implicit none
integer(4) :: numint0,nunit,nr,nfe,iIntr(nr*nunit,numint0),nunitx,nunity,nunitz
integer(4) :: i,i2,ir,ife,ipmx(nunitx,-1:1),ipmy(nunity,-1:1),ipmz(nunitz,-1:1)
integer(4) :: irpbc(nr,0:nunitx+1,0:nunity+1,0:nunitz+1),ifepbc(nfe,0:nunitx+1,0:nunity+1,0:nunitz+1)
integer(4) :: ix,iy,iz,ix2,iy2,iz2,ipt1,ipt2,icount,icount0,idx,idy,idz,numintr(nr*nunit),idim,ir0,ife0,nrtot,nfetot
real(8) :: ax,cx,abc,range,dis
integer(4) :: irunit(nr*nunit,3),ifeunit(nfe*nunit,3),irsub(nr*nunit),ifesub(nfe*nunit)
real(8) :: rrpt(nr,0:nunitx+1,0:nunity+1,0:nunitz+1,3),rfept(nfe,0:nunitx+1,0:nunity+1,0:nunitz+1,3)

iIntr=0.0d0
numintr=0
do ir=1,nrtot
   iz=irunit(ir,3)
   iy=irunit(ir,2)
   ix=irunit(ir,1)
   i=irsub(ir)
   
   icount=0
   do idz=-1,1
      do idy=-1,1
         do idx=-1,1
            do i2=1,nfe
               ife=ifepbc(i2,ix+idx,iy+idy,iz+idz)

               dis=0.0d0
               do idim=1,3
                  dis=dis+(rrpt(i,ix,iy,iz,idim)-rfept(i2,ix+idx,iy+idy,iz+idz,idim))**2
               end do
               dis=dsqrt(dis)

               if((dis.lt.range).and.(dis.ne.0.0d0).and.ife.ne.0) then
                  icount=icount+1
                  iIntr(ir,icount)=ife
                  numintr(ir)=icount
                  if(icount.gt.numint0) then
                     write(6,*) "the number of interacting pair is over limit"
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
