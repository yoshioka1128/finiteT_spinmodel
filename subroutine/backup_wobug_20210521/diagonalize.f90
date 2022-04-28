subroutine diagonalize(it,ip,nfe,Femag,temp,imax,KuFe,dtheta0,Hmag,nr,zH,zspin,eigen,RWORK,WORK,zunit0,part0,free0,freeloc0,theta0,phi0,eigen0,dmult)
implicit none
integer(4),parameter :: nrmax=30,itempmax=1000
real(8),parameter :: pi=dacos(-1.0d0)
real(8) :: theta,phi,free,eigen(imax),RWORK(3*imax-2),freeloc(nrmax),part(nrmax),temp,thetamag,phimag,Hmag,Femag
real(8) :: freeloc0(0:4,0:4,nrmax),part0(nrmax),free0,theta0,phi0,eigen0(imax),KuFe,freet0loc(nrmax),dtheta0,dmult
integer(4) :: ir,nr,jj,ii,LWORK,imax,INFO,it,ip,nfe,itemp
complex(kind(0d0)) :: zH(imax,imax,nrmax),zH0(imax,imax,nr),zspin(imax,imax,3),WORK(2*imax-1),zunit0(imax,imax,nr)

phi=pi*dble(ip)/8.0d0
if(ip.eq.3) phi=pi/12.0d0 ! for K3 tilda
theta=dtheta0*dble(it)

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
   if(it.eq.0) freet0loc(ir)=freeloc(ir)
   if(it.le.3.and.ip.le.4)  freeloc0(it,ip,ir)=freeloc(ir)-freet0loc(ir) 
   free=free+freeloc(ir)
end do
! free eng per unitcell
free=dmult*free-Hmag*Femag*(dsin(thetamag)*dsin(theta)*dcos(phimag-phi)+dcos(thetamag)*dcos(theta))&
     +KuFe*dsin(theta)**2
if(free.lt.free0) then
   zunit0=zH
   part0=part
   free0=free
   theta0=theta
   phi0=phi
   eigen0=eigen
end if

end subroutine diagonalize
