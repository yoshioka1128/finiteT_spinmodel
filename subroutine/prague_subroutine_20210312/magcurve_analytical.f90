subroutine magcurve_analytical(temp,sys,method,wwoK1Fe,volume,lambda,Alm,rdch,KuFe,Hex,Ms,ratom,nr,dmult,direction)
implicit none
integer(4),parameter :: nrmax=30,itempmax=1000,ianmax=2*10+1
real(8),parameter :: pi=dacos(-1.0d0)
real(8),parameter :: kb=1.3806488d-23 ! (J/K)
real(8),parameter :: muB=9.2740100783d-24 ! (Merg/kOe) or (kerg/Oe) or (J/T)
real(8),parameter :: mu0=4.0d0*pi*1.0d-7 ! (NA^-2)
real(8),parameter :: joule=10.0d0 ! (Merg/J)
real(8),parameter :: dtheta=0.001d0,dphi0=0.015d0
integer(4) :: iHmag,nr,n4f,itemp,m,i,j,k,l,imax,ii,it,ip,ll,ir,itmax,ipmax
real(8) :: temp,S1,L1,J1,Ms,Hmag,thetamag,phimag,theta,phi,dcoff,Dso,Dex(nrmax),abc,KuFe,g1,volume,dmult,dphi
real(8) :: JM(100,2),Stevens(6),extBJ(-1:8,nrmax),rdcoff(0:6),Alm(6,-6:6,nrmax),lambda
real(8) :: dmagRT0(3),thetaan(3),phian(3),Tmix1(6,nrmax),Tmix2(6,nrmax),free0(3),freeR(3),freeFe,gibbs(3),Hex(nrmax)
real(8) :: engmixcf1(nrmax),engmixcf2(nrmax),engcf(nrmax),coffmag(6,-6:6),xi(6),sc(3)
character(50) :: sys,ratom,method,wwoK1Fe,rdch,Tch,direction

dcoff=kb/volume*1.0d-6 ! convert from [K] to [MJ/m^3]
write(Tch,'(f5.1)') temp
Tch=adjustl(Tch)
if(sys.eq."SmFe12_L3") then
   call rareearth_SmL3(ratom,n4f,0,imax,JM,S1,L1,J1)
else 
   call rareearth(ratom,n4f,0,imax,JM,S1,L1,J1)
end if
g1=1.0d0+(J1*(J1+1.0d0)+S1*(S1+1.0d0)-L1*(L1+1.0d0))/(2.0d0*J1*(J1+1.0d0))

! matrix element of spherical tensor operators
if(sys.eq."SmFe12_L3") then
   call rdmatrix_SmL3(n4f,rdcoff,L1)
else
   call rdmatrix(n4f,rdcoff,L1) ! reduced matrix element
end if
call stfactor(Stevens,L1,J1,S1,rdcoff) ! Stevens factor
call xifactor(xi,L1,J1,S1,rdcoff) ! Stevens factor

! Spherical Tensor operator @ T=0
Dso=lambda*(J1+1.0d0)
do ir=1,nr
   Dex(ir)=abs(-2.0d0*S1/(J1+1.0d0)*Hex(ir))
   call extendedBJ2(extBJ(:,ir),Hex(ir),temp,J1,JM,S1)
   call TVJl(Dex(ir),Dso,J1,L1,S1,Tmix1(:,ir),Tmix2(:,ir),Hex(ir),extBJ(:,ir))
end do

dmagRT0=0.0d0
do ir=1,nr
   abc=(Dex(ir)/Dso)*(2.0d0*(L1+1.0d0))/(3.0d0*(J1+1.0d0))
   dmagRT0(1)=dmagRT0(1)+g1*extBJ(1,ir)
   dmagRT0(2)=dmagRT0(2)+g1*extBJ(1,ir)-abc*Tmix1(1,ir)
   dmagRT0(3)=dmagRT0(3)+g1*extBJ(1,ir)-abc*Tmix1(1,ir)-abc*Tmix2(1,ir)
end do

! direction of magnetic field
if(direction.eq."001") then
   thetamag=0.0d0
   phimag=0.0d0
else if(direction.eq."100") then
   thetamag=pi/2.0d0
   phimag=0.0d0
else if(direction.eq."110") then
   thetamag=pi/2.0d0
   phimag=pi/4.0d0
end if
! output files
open(799,file=""//trim(sys)//"_"//trim(method)//"_"//trim(wwoK1Fe)//"_magcurve_T"//trim(Tch)//"_"//trim(direction)//"_lowestJ_rd"//trim(rdch)//".txt")
open(800,file=""//trim(sys)//"_"//trim(method)//"_"//trim(wwoK1Fe)//"_magcurve_T"//trim(Tch)//"_"//trim(direction)//"_mixexcf1_rd"//trim(rdch)//".txt")
open(801,file=""//trim(sys)//"_"//trim(method)//"_"//trim(wwoK1Fe)//"_magcurve_T"//trim(Tch)//"_"//trim(direction)//"_mixexcf2_rd"//trim(rdch)//".txt")
do iHmag=0,400 ! to 400 kOe
   Hmag=dble(iHmag)*muB/joule/kb
   abc=0.0d0
   thetaan=0.0d0
   phian=0.0d0
   gibbs=0.0d0
   itmax=int(0.5d0*pi/dtheta)
   do it=0,itmax
      theta=0.5d0*pi*dble(it)/dble(itmax)
      if(it.ne.0) then
         dphi=dphi0/dsin(theta)
         ipmax=int(0.5d0*pi/dphi)+1
      else
         ipmax=1
      end if
      do ip=0,ipmax
         phi=0.5d0*pi*dble(ip)/dble(ipmax)
         call tlm(coffmag,theta,phi)
         engcf=0.0d0
         engmixcf1=0.0d0
         engmixcf2=0.0d0
         do ir=1,nr
            do l=2,6,2
               do m=-l,l
                  engcf(ir)=engcf(ir)+coffmag(l,m)*Alm(l,m,ir)*Stevens(l)*extBJ(l,ir)
                  abc=xi(l)*dble(l*(l+1))/(dble(2*l+1))*Dex(ir)/Dso*coffmag(l,m)*Alm(l,m,1)
                  engmixcf1(ir)=engmixcf1(ir)+abc*Tmix1(l,ir)
                  engmixcf2(ir)=engmixcf2(ir)+abc*Tmix2(l,ir)
               end do
            end do
         end do
         
         freeR=0.0d0
         do ir=1,nr
            freeR(1)=freeR(1)+engcf(1)
            freeR(2)=freeR(2)+engcf(1)+engmixcf1(1)
            freeR(3)=freeR(3)+engcf(1)+engmixcf1(1)+engmixcf2(1)
         end do
         freeFe=-Hmag*(Ms+dmagRT0(3)*dmult)*(dsin(thetamag)*dsin(theta)*dcos(phimag-phi)+dcos(thetamag)*dcos(theta))&
              +KuFe*dsin(theta)**2
         if(it.eq.0) then ! set energy level
            free0(1)=freeR(1)*dmult+freeFe
            free0(2)=freeR(2)*dmult+freeFe
            free0(3)=freeR(3)*dmult+freeFe
         end if
         do i=1,3
            if((freeR(i)*dmult+freeFe-free0(i))*dcoff.le.gibbs(i)) then
               gibbs(i)=(freeR(i)*dmult+freeFe-free0(i))*dcoff
               thetaan(i)=theta
               phian(i)=phi
            end if
         end do
      end do  ! ip theta
   end do   ! it theta
   write(6,*) "magnetic field [T]",Hmag*kb/muB
   write(6,*) "gradient & y-intercept",Hmag*(Ms+dmagRT0(3)*dmult)*dcoff,gibbs(3)
   write(6,*) "sin(theta),Mxy [muB], Gibbs [MJ/m^3]"
   write(6,*) dsin(thetaan(3)),dsin(thetaan(3))*(Ms+dmagRT0(3)*dmult)*muB/volume*mu0,gibbs(3)
   write(6,*) 
   do i=1,3
      sc(i)=dsin(thetaan(i))*dsin(thetamag)*dcos(phian(i)-phimag)+dcos(thetaan(i))*dcos(thetamag)
   end do
   write(799,*) Hmag*kb*joule/muB/10.0d0,sc(1)*(Ms+dmagRT0(1)*dmult)*muB/volume*mu0,gibbs(1)
   write(800,*) Hmag*kb*joule/muB/10.0d0,sc(2)*(Ms+dmagRT0(2)*dmult)*muB/volume*mu0,gibbs(2)
   write(801,*) Hmag*kb*joule/muB/10.0d0,sc(3)*(Ms+dmagRT0(3)*dmult)*muB/volume*mu0,gibbs(3)
end do ! iHmag
return
   
end subroutine magcurve_analytical
