program Ki_statistical
implicit none
integer(4),parameter :: nrmax=30,itempmax=1000,nJexmax=10
real(8),parameter :: pi=dacos(-1.0d0)
real(8),parameter :: kb=1.3806488d-23 ! (J/K)
real(8),parameter :: muB=9.2740100783d-24 ! (Merg/kOe) or (kerg/Oe) or (J/T)
real(8),parameter :: mu0=4.0d0*pi*1.0d-7 ! (NA^-2) 
real(8),parameter :: joule=10.0d0 ! (Merg/J)
integer(4) :: nJex,nfe,nr,n4f,itemp,ncfp,m,n,i,j,k,l,imax,ii,jj,INFO,it,ip,ll,ir,LWORK,ir0
real(8) :: J0,S1,L1,J1,lambda,Hmax,dH,Tmin,Tmax,dT,p,Femag,Tc,s,mFe0,MsT_Kuzmin
real(8) :: JM(100,2),SO(100),Clm(6,-6:6,nrmax),g(100),cclm(6,-6:6),WJ(0:nJexmax)
real(8) :: rdcoff(0:6),Alm(6,-6:6,nrmax),rl(6),Kloc(11,nrmax),dK1(2),dK2(3),dK3(4)
real(8) :: phij0(nrmax),thetaj0(nrmax),absj0(nrmax),phiS0(nrmax),thetaS0(nrmax),absS0(nrmax)
real(8) :: phiL0(nrmax),thetaL0(nrmax),absL0(nrmax)
real(8) :: Hex0(nrmax),thetatot,phitot,abstot,theta,phi,dcoff
real(8) :: abc,def,theta0,phi0,freeloc(nrmax),freeloc0(0:4,0:4,nrmax),freet0(nrmax)
real(8) :: part(nrmax),part0(nrmax),temp,free,free0,dtheta,dtheta0,Rmag0,Femag0
real(8) :: mFeT(itempmax),HexT(nrmax,0:itempmax),K1FeT(itempmax),magtot(3)
real(8) :: aleng,bleng,cleng,gmm,volume,dmult,rd,K1Fe,rK2Fe,rK3Fe,KuFe
real(8) :: stspin0(3,nrmax),stspinS0(3,nrmax),stspinL0(3,nrmax)
character(50) :: sys,nJexch,ratom,tatom,method,sch,wwoK1Fe,rdch
real(8),allocatable :: eigen(:),eigen0(:),RWORK(:),Ust(:,:,:,:),Vst(:,:,:)
complex(kind(0d0)),allocatable :: ZUst(:,:,:,:),WORK(:),zU(:,:,:),zH0(:,:,:),zspin(:,:,:),zorbit(:,:,:),zH(:,:,:)

! set paremeter
rewind(5)
read(5,*) sys
read(5,*) nJex ! 0:lowestJ, 1:LowestJ+secondJ
read(5,*) method
read(5,*) wwoK1Fe
read(5,*) Tmin,Tmax,dT
read(5,*) dtheta
write(6,*) "system=",sys
write(6,*) "number of excited states",nJex
write(6,*) "method=",method
write(6,*) "wwoK1Fe",wwoK1Fe
write(6,*) "temperature",Tmin,Tmax,dT,"(K)"
write(6,*) "increment of \theta",dtheta,"(rad)",dtheta*180.0d0/pi,"(deg.)"
write(nJexch,"(I2)") nJex
nJexch=adjustl(nJexch)

! input files
open(33,file="~/research/CFPs/Sm_lambda411K/cfp"//trim(sys)//"_"//trim(method)//".txt")
read(33,*) ratom,nr,tatom,nfe
read(33,*) dmult
read(33,*) aleng,bleng,cleng,gmm
read(33,*) lambda
read(33,*) mFe0,Tc,s,p
read(33,*) K1Fe,rK2Fe,rK3Fe
write(6,*) "system=",ratom,nr,tatom,nfe
write(6,*) "multiplicity",dmult
write(6,*) "(a,b,c)=",aleng,bleng,cleng,"(angstrom)"
write(6,*) "\lambda=",lambda,"(K)"
write(6,*) "Ms, Tc, s, p for Kuzmin fitting"
write(6,*) mFe0, Tc, s, p ! mFe0 [muB per #nfe]
write(6,*) "K1Fe0 [K/unitcell],K2Fe0/K1Fe0,K3Fe0/K1Fe0"
write(6,*) K1Fe,rK2Fe,rK3Fe
write(sch,"(f6.3)") s
write(6,*) "Kuz'min parameter",s
Alm=0.0d0
do ir=1,nr
   read(33,*)
   read(33,*) Hex0(ir),rd
   write(6,*) Hex0(ir),rd
   write(rdch,"(f5.1)") rd
   Hex0(ir)=rd*Hex0(ir)
   read(33,*) rl(2),rl(4),rl(6),ncfp
   write(6,*) "ncfp=",ncfp
   do i=1,ncfp
      read(33,*) abc,l,m
      if(l.lt.0) then
         write(6,*) "error: negative l"
         stop
      end if
      Alm(l,m,ir)=abc*rl(l)
   end do
end do
if(nr.gt.nrmax) then
   write(6,*) "error: nr > nrmax"
   stop
end if
rdch=adjustl(rdch)
volume=aleng*bleng*cleng*dsin(pi*gmm/180.0d0)*1.0d-30
dcoff=kb/volume*1.0d-6 ! convert from [K] to [MJ/m^3]

! Fe sublattice
K1FeT=0.0d0
open(287,file=""//trim(sys)//"_MsT_Kuzmin_s"//trim(sch)//".txt")
open(288,file=""//trim(sys)//"_K1FeT_Kuzmin_s"//trim(sch)//".txt")
do itemp=0,int((Tmax-Tmin)/dT)
   temp=Tmin+dble(itemp)*dT
   mFeT(itemp)=MsT_Kuzmin(temp,mFe0,Tc,s,p) ! in unit Tesla
   abc=mFeT(itemp)/mFe0
   def=abc**3+&
        8.0d0/7.0d0*rK2Fe*(abc**3-abc**10)+&
        8.0d0/7.0d0*rK3Fe*(abc**3-18.0d0/11.0d0*abc**10+7.0d0/11.0d0*abc**21)
   K1FeT(itemp)=K1Fe*def
   if(sys.eq."Sm2Fe17N3tri") K1FeT(itemp)=0.0d0
   write(287,*) temp,mu0*mFeT(itemp)*muB/volume ! in unit Tesla
   write(288,*) temp,K1FeT(itemp)/(aleng**2*cleng*1.0d-30)*kb*1.0d-6 ! [MJ/m^3]
   do ir=1,nr
      HexT(ir,itemp)=-Hex0(ir)*abc
   end do
end do

! set quantum number
if(sys.eq."SmFe12_L3") then
   call rareearth_SmL3(ratom,n4f,nJex,imax,JM,S1,L1,J1)
else 
   call rareearth(ratom,n4f,nJex,imax,JM,S1,L1,J1)
end if
if(imax.gt.100) then
   write(6,*) "imax too large"
   stop
end if

! spin orbit coupling
do i=1,imax
   SO(i)=0.5d0*lambda*(JM(i,1)*(JM(i,1)+1.0d0)-L1*(L1+1.0d0)-S1*(S1+1.0d0))
   g(i)=1.0d0+(JM(i,1)*(JM(i,1)+1.0d0)+S1*(S1+1.0d0)-L1*(L1+1.0d0))/(2.0d0*JM(i,1)*(JM(i,1)+1.0d0))
end do

! analytical calculations
call Ki_analytical(sys,method,wwoK1Fe,Tmin,Tmax,dT,volume,lambda,Alm,rdch,K1FeT,HexT,mFeT,Tc,ratom,nr,dmult)
ir0=1
call Tmix_analytical(ir0,sys,method,Tmin,Tmax,dT,SO,HexT(ir0,:),ratom,Tc)

! statistical calculations
allocate(eigen(imax),eigen0(imax),RWORK(3*imax-2),Ust(imax,imax,6,-6:6),Vst(imax,imax,-1:1))
allocate(zH0(imax,imax,nr),ZUst(imax,imax,6,-6:6),WORK(2*imax-1),zspin(imax,imax,3),zorbit(imax,imax,3),&
     zU(imax,imax,nr),zH(imax,imax,nr))
LWORK=2*imax-1

! matrix element of spherical tensor operators
if(sys.eq."SmFe12_L3") then
   call rdmatrix_SmL3(n4f,rdcoff,L1)
else
   call rdmatrix(n4f,rdcoff,L1) ! reduced matrix element
end if
call numfactor(cclm)
call CkqL(L1,S1,J1,Ust,imax,nJex) ! matrix element of tensor operator
call CkqS(L1,S1,J1,Vst,imax,nJex) ! matrix element of vector operator
call OkqL(ZUst,Ust,imax)
call moment_op(imax,L1,S1,zspin,zorbit,JM,Vst,Ust) ! matrix element of magnetic moment

! crystal field parameters for tensor operators
Clm=0.0d0
do ir=1,nr
   do l=2,6,2
      do m=-l,l
         def=1.0d0/cclm(l,m)*dsqrt((2.0d0*dble(l)+1.0d0)/(4.0d0*pi))*rdcoff(l)
         if(m.ne.0) then
            def=def/dsqrt(2.0d0)
         end if
         Clm(l,m,ir)=Alm(l,m,ir)*def
      end do
   end do
end do

! Hamiltonian
zH0=(0.0d0,0.0d0)
do ir=1,nr
   do jj=1,imax
      zH0(jj,jj,ir)=zH0(jj,jj,ir)+SO(jj)
      do ii=1,imax
         do l=2,6,2
            do m=-l,l
               zH0(ii,jj,ir)=zH0(ii,jj,ir)+ZUst(ii,jj,l,m)*Clm(l,m,ir) ! crystal field term
            end do
         end do
      end do
   end do
end do

! temperature dependence of Ki
open(26,file=""//trim(sys)//"_"//trim(method)//"_RE1-2_K1_nJex"//trim(nJexch)//"_rd"//trim(rdch)//".txt") 
open(27,file=""//trim(sys)//"_"//trim(method)//"_RE1-2_K2_nJex"//trim(nJexch)//"_rd"//trim(rdch)//".txt") 
open(28,file=""//trim(sys)//"_"//trim(method)//"_RE1-2_K3_nJex"//trim(nJexch)//"_rd"//trim(rdch)//".txt") 
open(31,file=""//trim(sys)//"_"//trim(method)//"_"//trim(wwoK1Fe)//"_K1tot_nJex"//trim(nJexch)//"_rd"//trim(rdch)//".txt")
open(32,file=""//trim(sys)//"_"//trim(method)//"_"//trim(wwoK1Fe)//"_K2tot_nJex"//trim(nJexch)//"_rd"//trim(rdch)//".txt") 
open(34,file=""//trim(sys)//"_"//trim(method)//"_"//trim(wwoK1Fe)//"_K3tot_nJex"//trim(nJexch)//"_rd"//trim(rdch)//".txt") 
open(131,file=""//trim(sys)//"_"//trim(method)//"_K1TM_nJex"//trim(nJexch)//"_rd"//trim(rdch)//".txt") !"# temperature (K), K1, K1' (
open(169,file=""//trim(sys)//"_"//trim(method)//"_"//trim(wwoK1Fe)//"_MsT_nJex"//trim(nJexch)//"_rd"//trim(rdch)//".txt") 
open(171,file=""//trim(sys)//"_"//trim(method)//"_"//trim(wwoK1Fe)//"_Fe_theta-phi_nJex"//trim(nJexch)//"_rd"//trim(rdch)//".txt") 
open(172,file=""//trim(sys)//"_"//trim(method)//"_"//trim(wwoK1Fe)//"_RE1-2_J_theta-phi_nJex"//trim(nJexch)//"_rd"//trim(rdch)//".txt") 
open(173,file=""//trim(sys)//"_"//trim(method)//"_"//trim(wwoK1Fe)//"_RE1-2_S_theta-phi_nJex"//trim(nJexch)//"_rd"//trim(rdch)//".txt") 
open(174,file=""//trim(sys)//"_"//trim(method)//"_"//trim(wwoK1Fe)//"_RE1-2_L_theta-phi_nJex"//trim(nJexch)//"_rd"//trim(rdch)//".txt") 
open(175,file=""//trim(sys)//"_"//trim(method)//"_"//trim(wwoK1Fe)//"_Mtot_theta-phi_nJex"//trim(nJexch)//"_rd"//trim(rdch)//".txt")
write(169,*) "# temperature (K), magnetization (per unit cell)"
write(171,*) "# temperature (K), theta (degree), phi (degree), abs"
write(172,*)"# temperature (K), theta (degree), phi (degree), abs for RE 1-2"
write(173,*) "# temperature (K), theta (degree), phi (degree), abs for RE 1-2"
write(174,*) "# temperature (K), theta (degree), phi (degree), abs for RE 1-2"
write(175,*) "# temperature (K), theta (degree), phi (degree), abs for magnetization"
if(nr.gt.2) then
   open(182,file=""//trim(sys)//"_"//trim(method)//"_"//trim(wwoK1Fe)//"_RE3-4_J_theta-phi_nJex"//trim(nJexch)//"_rd"//trim(rdch)//".txt") 
   open(183,file=""//trim(sys)//"_"//trim(method)//"_"//trim(wwoK1Fe)//"_RE3-4_S_theta-phi_nJex"//trim(nJexch)//"_rd"//trim(rdch)//".txt") 
   open(184,file=""//trim(sys)//"_"//trim(method)//"_"//trim(wwoK1Fe)//"_RE3-4_L_theta-phi_nJex"//trim(nJexch)//"_rd"//trim(rdch)//".txt") 
   write(182,*) "# temperature (K), theta (degree), phi (degree), abs for RE 3-4"
   write(183,*) "# temperature (K), theta (degree), phi (degree), abs for RE 3-4"
   write(184,*) "# temperature (K), theta (degree), phi (degree), abs for RE 3-4"
end if
if(int((Tmax-Tmin)/dT).gt.itempmax) then
   write(6,*) "error: itempmax over the number of T plot"
end if
! temperature dependence
do itemp=0,int((Tmax-Tmin)/dT)  ! loop for Temperature 
   temp=Tmin+dble(itemp)*dT
   Femag=mFeT(itemp) ! [muB/unitcell]
   if(wwoK1Fe.eq."wK1Fe") then
      KuFe=K1FeT(itemp) ! [K/unitcell]
   else if(wwoK1Fe.eq."woK1Fe") then
      KuFe=0.0d0
   else
      write(6,*) "read wwoK1Fe error",wwoK1Fe
      stop
   end if

! free energy at theta 0
   part=0.0d0
   do ir=1,nr
      call Hexch(0.0d0,0.0d0,imax,zH(:,:,ir),zH0(:,:,ir),zspin,HexT(ir,itemp))
      call zheev('V','L',imax,zH(:,:,ir),imax,eigen,WORK,LWORK,RWORK,INFO)
      if(temp.eq.0.0d0) then
         freet0(ir)=eigen(1) 
      else
         do  ii=1,imax ! state sum
            part(ir)=part(ir)+dexp(-(eigen(ii)-eigen(1))/temp)
         end do
         freet0(ir)=eigen(1)-temp*dlog(part(ir))
      end if
   end do ! ir   

! angular dpendent free energy
   free0=1.0d10
   part0=0.0d0
   theta0=0.0d0
   phi0=0.0d0
   dtheta0=0.5d0*pi/dble(int(0.5d0*pi/dtheta))
   do ip=0,4 ! \phi=pi*dble(ip)/8.0d0
      phi=pi*dble(ip)/8.0d0
      if(ip.eq.3) phi=pi/12.0d0 ! for K3 tilda
      do it=0,int(0.5d0*pi/dtheta) ! take three point
         theta=dtheta0*dble(it)
         part=0.0d0
         free=0.0d0
         do ir=1,nr
            call Hexch(theta,phi,imax,zH(:,:,ir),zH0(:,:,ir),zspin,HexT(ir,itemp))
            call zheev('V','L',imax,zH(:,:,ir),imax,eigen,WORK,LWORK,RWORK,INFO)
            if(temp.eq.0.0d0) then
               freeloc(ir)=eigen(1)-freet0(ir)
            else
               do  ii=1,imax ! state sum
                  part(ir)=part(ir)+dexp(-(eigen(ii)-eigen(1))/temp)
               end do
               freeloc(ir)=eigen(1)-temp*dlog(part(ir))-freet0(ir)
            end if
            if(it.le.3.and.ip.le.4)  freeloc0(it,ip,ir)=freeloc(ir)
            free=free+freeloc(ir)
         end do ! do ir 
         free=dmult*free+KuFe*dsin(theta)**2
         if(free.lt.free0) then
            zU=zH
            part0=part
            free0=free
            theta0=theta
            phi0=phi
            eigen0=eigen
         end if
      end do ! for \theta it
   end do ! for \phi ip
   
   do ir=1,nr
      call momentR_exp(temp,imax,zU(:,:,ir),zspin,zorbit,part0(ir),eigen0,stspin0(:,ir),stspinS0(:,ir),stspinL0(:,ir))
      call vector_to_angle(stspin0(:,ir),thetaj0(ir),phij0(ir),absj0(ir))
      call vector_to_angle(stspinS0(:,ir),thetaS0(ir),phiS0(ir),absS0(ir))
      call vector_to_angle(stspinL0(:,ir),thetaL0(ir),phiL0(ir),absL0(ir))
   end do
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
   call vector_to_angle(magtot,thetatot,phitot,abstot)
   call anisotropy_tetra(temp,nr,freeloc0,dtheta0,dK1,dK2,dK3,Kloc)

   ! output
   write(6,*) 
   write(6,*) "temperature [K]",temp
   write(6,*) "energy [K]",free0
   do ir=1,nr
      write(6,"(a20,3f12.5)") "RE theta-phi-abs",thetaj0(ir)/pi*180.0d0,phij0(ir)/pi*180.0d0,absj0(ir)
   end do
   write(6,"(a20,3f12.5)") "Fe theta-phi-abs",theta0/pi*180.0d0,phi0/pi*180.0d0,Femag/(dmult*dble(nfe))
   write(6,"(a20,3f12.5)") "total theta-phi-abs",thetatot/pi*180.0d0,phitot/pi*180.0d0,abstot
   
   write(169,'(F10.2,3F15.5)') temp,magtot(1),magtot(2),magtot(3)
   write(171,"(4f12.5)") temp,theta0/pi*180.0d0,phi0/pi*180.0d0,Femag/(dmult*dble(nfe))
   write(172,"(7f12.5)") temp,thetaj0(1)/pi*180.0d0,phij0(1)/pi*180.0d0,absj0(1),&
        thetaj0(2)/pi*180.0d0,phij0(2)/pi*180.0d0,absj0(2)
   write(173,"(7f12.5)") temp,thetaS0(1)/pi*180.0d0,phiS0(1)/pi*180.0d0,absS0(1),&
        thetaS0(2)/pi*180.0d0,phiS0(2)/pi*180.0d0,absS0(2)
   write(174,"(7f12.5)") temp,thetaL0(1)/pi*180.0d0,phiL0(1)/pi*180.0d0,absL0(1),&
        thetaL0(2)/pi*180.0d0,phiL0(2)/pi*180.0d0,absL0(2)
   if(nr.gt.2) then
      write(182,"(7f12.5)") temp,thetaj0(3)/pi*180.0d0,phij0(3)/pi*180.0d0,absj0(3),&
           thetaj0(4)/pi*180.0d0,phij0(4)/pi*180.0d0,absj0(4)
      write(183,"(7f12.5)") temp,thetaS0(3)/pi*180.0d0,phiS0(3)/pi*180.0d0,absS0(3),&
           thetaS0(4)/pi*180.0d0,phiS0(4)/pi*180.0d0,absS0(4)
      write(184,"(7f12.5)") temp,thetaL0(3)/pi*180.0d0,phiL0(3)/pi*180.0d0,absL0(3),&
           thetaL0(4)/pi*180.0d0,phiL0(4)/pi*180.0d0,absL0(4)
   end if
   write(175,"(4f12.5)") temp,thetatot/pi*180.0d0,phitot/pi*180.0d0,mu0*abstot*muB/volume
   
   write(26,'(F10.2,4F20.5)') temp,dK1(1)/dble(nr)
   write(27,'(F10.2,6F20.5)') temp,dK2(1)/dble(nr)
   write(28,'(F10.2,8F20.5)') temp,Kloc(8,1),Kloc(8,2),Kloc(9,1),Kloc(9,2),Kloc(10,1),Kloc(10,2),Kloc(11,1),Kloc(11,2) ! K3333 
   write(31,'(F10.2,2F20.10)') temp,(dK1(1)*dmult+KuFe)*dcoff,dK1(2)*dmult*dcoff
   write(32,'(F10.2,3F20.10)') temp,dK2(1)*dmult*dcoff,dK2(2)*dmult*dcoff,dK2(3)*dmult*dcoff
   write(34,'(F10.2,4F20.10)') temp,dK3(1)*dmult*dcoff
   write(131,'(F10.2,4F20.10)') temp,KuFe*dcoff

   write(6,*) "K1tot [K, MJ/m^3]",dK1(1)/dble(nr),(dK1(1)*dmult+KuFe)*dcoff
   write(6,*) "K2tot [K, MJ/m^3]",dK2(1)/dble(nr),dK2(1)*dmult*dcoff
   write(6,*) "K3tot [K, MJ/m^3]",dK3(1)/dble(nr),dK3(1)*dmult*dcoff
   write(6,*) "K1TM [K, MJ/m^3]",KuFe/dble(nfe),KuFe*dcoff
end do   ! do itemp
   
end program Ki_statistical


