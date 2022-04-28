subroutine numeng(theta,phi,nfe,Femag,temp,imax,KuFe,Hmag,nr,zH,eigen,RWORK,WORK,free,dmult)
implicit none
integer(4),parameter :: nrmax=30,itempmax=1000
real(8),parameter :: pi=dacos(-1.0d0)
real(8) :: theta,phi,free,eigen(imax),RWORK(3*imax-2),freeloc(nrmax),part(nrmax),temp,thetamag,phimag,Hmag,Femag
real(8) :: KuFe,dmult
integer(4) :: ir,nr,jj,ii,LWORK,imax,INFO,nfe,itemp
complex(kind(0d0)) :: zH(imax,imax,nrmax),zH0(imax,imax,nr),WORK(2*imax-1)

part=0.0d0
free=0.0d0
do ir=1,nr
   LWORK=2*imax-1
   call zheev('V','L',imax,zH(1,1,ir),imax,eigen,WORK,LWORK,RWORK,INFO)
   if(temp.eq.0.0d0) then
      freeloc(ir)=eigen(1) 
   else
      do  ii=1,imax ! state sum
         part(ir)=part(ir)+dexp(-(eigen(ii)-eigen(1))/temp)
      end do
      freeloc(ir)=eigen(1)-temp*dlog(part(ir))
   end if
   free=free+freeloc(ir)
end do
! free eng per unitcell
free=dmult*free-Hmag*Femag*(dsin(thetamag)*dsin(theta)*dcos(phimag-phi)+dcos(thetamag)*dcos(theta))&
     +KuFe*dsin(theta)**2

end subroutine numeng
