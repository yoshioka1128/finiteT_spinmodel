subroutine moment_exp(temp,nr,imax,zunit0,zspin,zorbit,thetamag,phimag,theta0,phi0,Femag,dmult,&
     part0,eigen0,stspin0,stspinS0,stspinL0,magtot,Rmag0,Femag0)
implicit none
integer(4),parameter :: nrmax=30
real(8),parameter :: mu0=4.0d0*dacos(-1.0d0)*1.0d-7 ! (NA^-2)
real(8),parameter :: muB=9.2740100783d-24 ! (Merg/kOe) or (kerg/Oe) or (J/T)
real(8),parameter :: pi=dacos(-1.0d0)
integer(4) :: ir,nr,ll,m,imax,jj,ii
real(8) :: stspin0(3,nrmax),stspinS0(3,nrmax),stspinL0(3,nrmax),eigen0(imax),part0(nrmax),temp,dmult
real(8) :: def,ghi,Rmag0,Femag,Femag0,magtot(3),thetamag,phimag,theta0,phi0
complex(kind(0d0)) :: zunit0(imax,imax,nr),zspin(imax,imax,3),zorbit(imax,imax,3)

! angular momentum
stspin0=0.0d0
stspinS0=0.0d0
stspinL0=0.0d0
do ir=1,nr
   do ll=1,3
      do m=1,imax
         def=0.0d0
         ghi=0.0d0
         do jj=1,imax
            do ii=1,imax
               def=def-2.0d0*dreal(conjg(zunit0(ii,m,ir))*zspin(ii,jj,ll)*zunit0(jj,m,ir))
               ghi=ghi-dreal(conjg(zunit0(ii,m,ir))*zorbit(ii,jj,ll)*zunit0(jj,m,ir))
            end do
         end do
         if(temp.eq.0.0d0) then
            stspinS0(ll,ir)=def
            stspinL0(ll,ir)=ghi
            exit
         else
            stspinS0(ll,ir)=stspinS0(ll,ir)+def*dexp(-(eigen0(m)-eigen0(1))/temp)/part0(ir)                        
            stspinL0(ll,ir)=stspinL0(ll,ir)+ghi*dexp(-(eigen0(m)-eigen0(1))/temp)/part0(ir)                        
         end if
      end do ! m energy level
      stspin0(ll,ir)=stspinS0(ll,ir)+stspinL0(ll,ir)
   end do ! for spin component
end do ! for ir=1,2

! magnetization of Fe and RE along the applied field direction
Femag0=Femag*(dsin(thetamag)*dsin(theta0)*dcos(phimag-phi0)+dcos(thetamag)*dcos(theta0))
Rmag0=0.0d0
do ir=1,nr
   Rmag0=Rmag0+(dsin(thetamag)*(stspin0(1,ir)*dcos(phimag)+stspin0(2,ir)*dsin(phimag))+stspin0(3,ir)*dcos(thetamag))
end do
Rmag0=Rmag0*dmult

! total magnetization
magtot=0.0d0
magtot(1)=Femag*dsin(theta0)*dcos(phi0)
magtot(2)=Femag*dsin(theta0)*dsin(phi0)
magtot(3)=Femag*dcos(theta0)
do ir=1,nr
   do ll=1,3
      magtot(ll)=magtot(ll)+dmult*stspin0(ll,ir)
   end do
end do

end subroutine moment_exp
