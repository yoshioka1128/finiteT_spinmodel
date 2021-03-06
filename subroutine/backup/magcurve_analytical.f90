subroutine magcurve_analytical(temp,sys,method,wwoK1TM,volume,lambda,Alm,rdch,KuTM,Hex,Ms,ratom,nr,dmult,direction,Hmax,dH)
implicit none
integer(4),parameter :: nrmax=30,itempmax=1000,mm=100
real(8),parameter :: pi=dacos(-1.0d0)
real(8),parameter :: kb=1.3806488d-23 ! (J/K)
real(8),parameter :: muB=9.2740100783d-24 ! (Merg/kOe) or (kerg/Oe) or (J/T)
real(8),parameter :: mu0=4.0d0*pi*1.0d-7 ! (NA^-2)
real(8),parameter :: joule=10.0d0 ! (Merg/J)
real(8),parameter :: dtheta0=0.001d0,dphi0=0.015d0
integer(4) :: iHmag,nr,n4f,itemp,m,i,j,k,l,imax,ii,it,ip,ll,ir,itmax,ipmax
real(8) :: temp,S1,L1,J1,Ms,Hmag,thetamag,phimag,theta,phi,dcoff,Dso,Dex(nrmax)
real(8) :: abc,KuTM,g1,volume,dmult,dphi,freeZ(3),gibbs0(3),gibbs(3),Hex(nrmax)
real(8) :: JM(mm,2),Stevens(6),extBJ(-1:8,nrmax),rdcoff(0:6),Alm(6,-6:6,nrmax),lambda
real(8) :: dmagRT0(3),thetaan(3),phian(3),Tmix1(6,nrmax),Tmix2(6,nrmax),freet0(3),freeR(3)
real(8) :: engmixcf1(nrmax),engmixcf2(nrmax),engcf(nrmax),coffmag(6,-6:6),xi(6),sc(3),Hmax,dH
character(50) :: sys,ratom,method,wwoK1TM,rdch,Tch,direction


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
   abc=(2.0d0*(L1+1.0d0))/(3.0d0*(J1+1.0d0))
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
open(799,file=""//trim(sys)//"_"//trim(method)//"_"//trim(wwoK1TM)//"_&
     magcurve_T"//trim(Tch)//"_"//trim(direction)//"_lowestJ_rd"//trim(rdch)//".txt")
open(800,file=""//trim(sys)//"_"//trim(method)//"_"//trim(wwoK1TM)//"_&
     magcurve_T"//trim(Tch)//"_"//trim(direction)//"_mixexcf1_rd"//trim(rdch)//".txt")
open(801,file=""//trim(sys)//"_"//trim(method)//"_"//trim(wwoK1TM)//"_&
     magcurve_T"//trim(Tch)//"_"//trim(direction)//"_mixexcf2_rd"//trim(rdch)//".txt")
do iHmag=0,int(Hmax/dH) ! to 400 kOe
   Hmag=dble(iHmag)*dH*muB/joule/kb
   abc=0.0d0
   thetaan=0.0d0
   phian=0.0d0
   gibbs0=0.0d0
   itmax=int(0.5d0*pi/dtheta0)
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
                  engcf(ir)=engcf(ir)+coffmag(l,m)*Stevens(l)*Alm(l,m,ir)*extBJ(l,ir)
                  abc=xi(l)*dble(l*(l+1))/(dble(2*l+1))*coffmag(l,m)*Alm(l,m,ir)
                  engmixcf1(ir)=engmixcf1(ir)+abc*Tmix1(l,ir)
                  engmixcf2(ir)=engmixcf2(ir)+abc*Tmix2(l,ir)
               end do
            end do
         end do
         
         freeR=0.0d0
         do ir=1,nr
            freeR(1)=freeR(1)+engcf(ir)
            freeR(2)=freeR(2)+engcf(ir)+engmixcf1(ir)
            freeR(3)=freeR(3)+engcf(ir)+engmixcf1(ir)+engmixcf2(ir)
         end do
         do i=1,3
            freeZ(i)=-Hmag*(Ms+dmagRT0(i)*dmult)*&
                 (dsin(thetamag)*dsin(theta)*dcos(phimag-phi)+dcos(thetamag)*dcos(theta))&
                 +KuTM*dsin(theta)**2
         end do
         if(it.eq.0) then ! set energy level
            do i=1,3
               freet0(i)=freeR(i)*dmult+freeZ(i)
            end do
         end if
         do i=1,3
            gibbs(i)=(freeR(i)*dmult+freeZ(i)-freet0(i))*dcoff
            if(gibbs(i).le.gibbs0(i)) then
               gibbs0(i)=gibbs(i)
               thetaan(i)=theta
               phian(i)=phi
            end if
         end do
      end do  ! ip theta
   end do   ! it theta
   do i=1,3
      sc(i)=dsin(thetaan(i))*dsin(thetamag)*dcos(phian(i)-phimag)+dcos(thetaan(i))*dcos(thetamag)
   end do
   write(799,*) dble(iHmag)*dH/10.0d0,sc(1)*(Ms+dmagRT0(1)*dmult)*muB/volume*mu0,gibbs0(1)
   write(800,*) dble(iHmag)*dH/10.0d0,sc(2)*(Ms+dmagRT0(2)*dmult)*muB/volume*mu0,gibbs0(2)
   write(801,*) dble(iHmag)*dH/10.0d0,sc(3)*(Ms+dmagRT0(3)*dmult)*muB/volume*mu0,gibbs0(3)
end do ! iHmag
return
   
end subroutine magcurve_analytical
