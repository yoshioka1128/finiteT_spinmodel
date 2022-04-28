program eng_level
implicit none
integer(4),parameter :: nrmax=30,itempmax=1000,nJexmax=10
real(8),parameter :: pi=dacos(-1.0d0)
real(8),parameter :: kb=1.3806488d-23 ! (J/K)
real(8),parameter :: muB=9.2740100783d-24 ! (Merg/kOe) or (kerg/Oe) or (J/T)
real(8),parameter :: mu0=4.0d0*pi*1.0d-7 ! (NA^-2)

integer(4) :: nJex,nfe,nr,n4f,itemp,ncfp,m,i,j,k,l,imax,ii,jj,INFO,it,ip,ll,ir,LWORK
real(8) :: S1,L1,J1,lambda,Hmax,dH,p,Tc,s,mFe0,MsT_Kuzmin,Femag,KuFe,abc,def,dcoff
real(8) :: aleng,bleng,cleng,gmm,volume,dmult,rd,K1Fe,rK2Fe,rK3Fe,Dso,Dex,freeFe
real(8) :: JM(100,2),SO(100),Stevens(6),Clm(6,-6:6,nrmax),g(100),cclm(6,-6:6),rl(6)
real(8) :: extBJ(-1:8,nrmax),rdcoff(0:6),Alm(6,-6:6,nrmax),Tmix1(6,nrmax),Tmix2(6,nrmax)
real(8) :: Hex0(nrmax),theta,phi,mFeT(itempmax),Hex(nrmax),K1FeT(itempmax)
real(8) :: part(nrmax),temp,free,dtheta,dtheta0,dmagRT0(3),thetaan(3),free0(3),freeR0(3),freeR(3)
real(8) :: coffmix1(6,-1:8),coffmix2(6,-1:8),coffmix01(6,-1:8),coffmix02(6,-1:8),coffmag(6,-6:6),xi(6)
real(8) :: engmixcf1(nrmax),engmixcf2(nrmax),engcf(nrmax),WJ(0:nJexmax)
real(8),allocatable :: matrix(:,:),eigen(:),englvcf(:),eigen0(:),RWORK(:),Ust(:,:,:,:),Vst(:,:,:),CJst(:,:,:)
real(8),allocatable :: englvmixex1(:),englvmixex2(:),englvmixcf1(:),englvmixcf2(:)
complex(kind(0d0)),allocatable :: ZUst(:,:,:,:),WORK(:),zH(:,:,:),zH0(:,:,:),zspin(:,:,:),zorbit(:,:,:)
character(50) :: sys,nJexch,ratom,tatom,method,wwoK1Fe,rdch,Tch

! set paremeter
read(5,*) sys
read(5,*) nJex ! 0:lowestJ, 1:LowestJ+secondJ
read(5,*) method
read(5,*) wwoK1Fe
read(5,*) temp
read(5,*) dtheta
write(6,*) "system=",sys
write(6,*) "number of excited states",nJex
write(6,*) "method=",method
write(6,*) "wwoK1Fe",wwoK1Fe
write(6,*) "temperature",temp
write(6,*) "increment of \theta",dtheta,"(rad)",dtheta*180.0d0/pi,"(deg.)"
write(nJexch,"(I1)") nJex
write(Tch,"(F5.1)") temp

! input files
open(33,file="~/research/CFPs/Sm_lambda411K/cfp"//trim(adjustl(sys))//"_"//trim(adjustl(method))//".txt")
read(33,*) ratom,nr,tatom,nfe
write(6,*) "system=",ratom,nr,tatom,nfe
read(33,*) dmult
write(6,*) "multiplicity",dmult
read(33,*) aleng,bleng,cleng,gmm
write(6,*) "(a,b,c)=",aleng,bleng,cleng,"(angstrom)"
read(33,*) lambda
write(6,*) "\lambda=",lambda,"(K)"
read(33,*) mFe0,Tc,s,p
write(6,*) "Ms, Tc, s, p for Kuzmin fitting"
write(6,*) mFe0, Tc, s, p ! mFe0 [muB per #nfe]
read(33,*) K1Fe,rK2Fe,rK3Fe
write(6,*) "K1Fe0 [K/unitcell],K2Fe0/K1Fe0,K3Fe0/K1Fe0"
write(6,*) K1Fe,rK2Fe,rK3Fe
write(6,*) "Kuz'min parameter",s
if(nr.gt.nrmax) then
   write(6,*) "error: nr > nrmax"
   stop
end if

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
rdch=adjustl(rdch)
volume=aleng*bleng*cleng*dsin(pi*gmm/180.0d0)*1.0d-30
dcoff=kb/volume*1.0d-6
write(6,*) "Fe sublattice magnetization @T=0 [T]"
write(6,*) mu0*mFe0*muB/volume

call rareearth(ratom,n4f,nJex,imax,JM,S1,L1,J1)

! spin orbit coupling
do i=1,imax
   SO(i)=0.5d0*lambda*(JM(i,1)*(JM(i,1)+1.0d0)-L1*(L1+1.0d0)-S1*(S1+1.0d0))
   g(i)=1.0d0+(JM(i,1)*(JM(i,1)+1.0d0)+S1*(S1+1.0d0)-L1*(L1+1.0d0))/(2.0d0*JM(i,1)*(JM(i,1)+1.0d0))
end do
abc=J1+1.0d0
write(6,*) g(int(2.0d0*J1)+2),1.0d0+(abc*(abc+1.0d0)+S1*(S1+1.0d0)-L1*(L1+1.0d0))/(2.0d0*abc*(abc+1.0d0))

allocate(matrix(imax,imax),eigen(imax),englvcf(imax),eigen0(imax),RWORK(3*imax-2),Ust(imax,imax,6,-6:6))
allocate(englvmixex1(imax),englvmixex2(imax),englvmixcf1(imax),englvmixcf2(imax))
allocate(CJst(imax,-1:8,-8:8),Vst(imax,imax,-1:1))
allocate(zH0(imax,imax,nr),zH(imax,imax,nr),ZUst(imax,imax,6,-6:6),WORK(2*imax-1),zspin(imax,imax,3),zorbit(imax,imax,3))
LWORK=2*imax-1

! matrix element of spherical tensor operators
call rdmatrix(n4f,rdcoff,L1) ! reduced matrix element
call stfactor(Stevens,L1,J1,S1,rdcoff) ! Stevens factor
call xifactor(xi,L1,J1,S1,rdcoff) ! Stevens factor
call numfactor(cclm)

! crystal field parameters
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
call CkqL(L1,S1,J1,Ust,imax,nJex) ! matrix element of tensor operator
call CkqS(L1,S1,J1,Vst,imax,nJex) ! matrix element of vector operator
call OkqL(ZUst,Ust,imax)
call moment_op(imax,L1,S1,zspin,zorbit,JM,Vst,Ust) ! matrix element of magnetic moment

! TM sublattice
K1FeT=0.0d0
Femag=MsT_Kuzmin(temp,mFe0,Tc,s,p) ! in unit Tesla
abc=Femag/mFe0
def=abc**3+&
     8.0d0/7.0d0*rK2Fe*(abc**3-abc**10)+&
     8.0d0/7.0d0*rK3Fe*(abc**3-18.0d0/11.0d0*abc**10+7.0d0/11.0d0*abc**21)
KuFe=K1Fe*def
do ir=1,nr
   Hex(ir)=-Hex0(ir)*abc
end do
if(wwoK1Fe.eq."woK1Fe") KuFe=0.0d0

! numerical calculations for energy level at T=0 and theta dependent Free energy
! spin-orbit term
zH0=0.0d0
do ir=1,nr
   do ii=1,imax
      zH0(ii,ii,ir)=zH0(ii,ii,ir)+SO(ii) ! spin orbit term
   end do
end do
open(270,file=""//trim(adjustl(sys))//"_"//trim(adjustl(method))//"_englevel_&
     T"//trim(adjustl(Tch))//"_diag_so_nJex"//trim(adjustl(nJexch))//"_rd"//trim(adjustl(rdch))//".txt")
do ii=1,imax
   write(270,*) ii,real(zH0(ii,ii,1))-SO(1),0 ! Hso only
   write(270,*) ii,real(zH0(ii,ii,1))-SO(1),1
   write(270,*)
end do
! spin-orbit and exchange term
open(271,file=""//trim(adjustl(sys))//"_"//trim(adjustl(method))//"_englevel_&
     T"//trim(adjustl(Tch))//"_diag_soex_nJex"//trim(adjustl(nJexch))//"_rd"//trim(adjustl(rdch))//".txt")
dtheta0=0.5d0*pi/dble(int(0.5d0*pi/dtheta))
call Hexch(0.0d0,0.0d0,imax,zH(:,:,1),zH0(:,:,1),zspin,Hex(1))
call zheev('V','L',imax,zH(:,:,1),imax,eigen,WORK,LWORK,RWORK,INFO)
do ii=1,imax
   write(271,*) ii,eigen(ii)-SO(1),0
   write(271,*) ii,eigen(ii)-SO(1),1
   write(271,*)
end do
! spin-orbit, exchange, and CF term
do ir=1,nr
   do jj=1,imax
      do ii=1,imax
         abc=0.0d0
         do l=2,6,2
            do m=-l,l
               abc=abc+ZUst(ii,jj,l,m)*Clm(l,m,ir) ! crystal field term
            end do
         end do
         zH0(ii,jj,ir)=zH0(ii,jj,ir)+abc
      end do
   end do
end do
open(272,file=""//trim(adjustl(sys))//"_"//trim(adjustl(method))//"_englevel_&
     T"//trim(adjustl(Tch))//"_diag_soexcf_nJex"//trim(adjustl(nJexch))//"_rd"//trim(adjustl(rdch))//".txt")
call Hexch(theta,phi,imax,zH(:,:,1),zH0(:,:,1),zspin,Hex(1))
call zheev('V','L',imax,zH(:,:,1),imax,eigen,WORK,LWORK,RWORK,INFO)      
do ii=1,imax
   write(272,*) ii,eigen(ii)-SO(1),0 ! Hso+Hex+Hcf
   write(272,*) ii,eigen(ii)-SO(1),1
   write(272,*)
end do




! angular dpendence of MA energy
open(285,file=""//trim(adjustl(sys))//"_"//trim(adjustl(method))//"_thetadept_&
     phi0_single_eng_T"//trim(adjustl(Tch))//"_&
     nJex"//trim(adjustl(nJexch))//"_rd"//trim(adjustl(rdch))//".txt")
open(286,file=""//trim(adjustl(sys))//"_"//trim(adjustl(method))//"_"//trim(adjustl(wwoK1Fe))//"_thetadept_&
     phi0_tot_eng_T"//trim(adjustl(Tch))//"_&
     nJex"//trim(adjustl(nJexch))//"_rd"//trim(adjustl(rdch))//".txt")
do it=0,int(0.5d0*pi/dtheta)
   theta=dtheta0*dble(it)
   freeFe=KuFe*dsin(theta)**2
   free=0.0d0
   do ir=1,nr
      call Hexch(theta,0.0d0,imax,zH(:,:,ir),zH0(:,:,ir),zspin,Hex(ir))
      call zheev('V','L',imax,zH(:,:,ir),imax,eigen,WORK,LWORK,RWORK,INFO)      
      if(temp.eq.0.0d0) then
         free=free+eigen(1) 
      else
         do  ii=1,imax ! state sum
            part(ir)=part(ir)+dexp(-(eigen(ii)-eigen(1))/temp)
         end do
         free=free+eigen(1)-temp*dlog(part(ir))
      end if
   end do ! do ir
   free=dmult*free+KuFe*dsin(theta)**2
   freeR(1)=(free-freeFe)/dmult/dble(nr)
   if(it.eq.0) then
      free0(1)=free
      freeR0(1)=freeR(1)
   end if
   write(285,*) theta,freeR(1)-freeR0(1)
   write(286,*) theta,(free-free0(1))*dcoff
end do ! for \phi it
open(385,file=""//trim(adjustl(sys))//"_"//trim(adjustl(method))//"_phidept_theta0.5pi_&
     single_eng_T"//trim(adjustl(Tch))//"_nJex"//trim(adjustl(nJexch))//"_rd"//trim(adjustl(rdch))//".txt")
open(386,file=""//trim(adjustl(sys))//"_"//trim(adjustl(method))//"_"//trim(adjustl(wwoK1Fe))//"_phidept_theta0.5pi_&
     tot_eng_T"//trim(adjustl(Tch))//"_nJex"//trim(adjustl(nJexch))//"_rd"//trim(adjustl(rdch))//".txt")
theta=pi/2.0d0
do ip=0,100 ! \phi=pi*dble(ip)/8.0d0
   phi=pi/2.0d0*dble(ip)/100.0d0
   freeFe=KuFe*dsin(theta)**2
   free=0.0d0
   do ir=1,nr
      call Hexch(pi/2.0d0,phi,imax,zH(:,:,ir),zH0(:,:,ir),zspin,Hex(ir))
      call zheev('V','L',imax,zH(1,1,ir),imax,eigen,WORK,LWORK,RWORK,INFO)
      if(temp.eq.0.0d0) then
         free=free+eigen(1) 
      else
         do  ii=1,imax ! state sum
            part(ir)=part(ir)+dexp(-(eigen(ii)-eigen(1))/temp)
         end do
         free=free+eigen(1)-temp*dlog(part(ir))
      end if
   end do ! do ir
   free=dmult*free+KuFe*dsin(theta)**2
   freeR(1)=(free-freeFe)/dmult/dble(nr)
   write(385,*) phi,freeR(1)-freeR0(1)
   write(386,*) phi,(free-free0(1))*def
end do ! for \phi ip






! Kuz'min linear theory
call CkqJ(J1,CJst,JM) ! matrix element of vector operator
Dex=abs(-2.0d0*S1/(J1+1.0d0)*Hex(1))
Dso=lambda*(J1+1.0d0)
write(6,*) "epsilon 1PT",3.0d0*(4.0d0*J1-3.0d0)*Alm(2,0,1)*Stevens(2)/(2.0d0*Dex)
write(6,*) "epsilon mixing 1",12.0d0*300.0d0/(2.0d0*J1-1.0d0)/Dso
write(6,*) "epsilon mixing 2",-(Dex/Dso)**2*6.0d0/5.0d0*(2.0d0*J1+3.0d0)/2.0d0/(J1**2*(2.0d0*J1-1.0d0)*(J1+1.0d0)*(2.0d0*J1+3.0d0)/(60.0d0*J1**3)*(J1*Dex/300.0d0)**2)*(abs(g(7)-g(1))/(g(1)-1.0d0))*J1*(J1+1.0d0)/3.0d0,&
     -12.0d0*(L1+S1+1.0d0)/(Dso**2*(2.0d0*J1-1.0d0)*(J1+2.0d0)*S1)*300.0d0**2
do ir=1,nr
   call extendedBJ2(extBJ(:,ir),Hex(ir),temp,J1,JM,S1)
   call TVJl(Dex,Dso,J1,L1,S1,Tmix1(:,ir),Tmix2(:,ir),Hex(ir),extBJ(:,ir))
end do
coffmix1=0.0d0
coffmix2=0.0d0
do l=1,6
   coffmix1(l,l-1)=(2.0d0*J1+dble(l)+1.0d0)/2.0d0
   coffmix1(l,l+1)=-2.0d0/(2.0d0*J1+dble(l)+2.0d0)
   do j=-1,1,2
      coffmix2(l,l-1+j)=coffmix2(l,l-1+j)&
           -coffmix1(l,l+j)*(dble(l+j)*(2.0d0*J1-dble(l+j)+1.0d0)*(2.0d0*J1+dble(l+j)+1.0d0)/(4.0d0*(2.0d0*dble(l+j)+1.0d0)))
      coffmix2(l,l+1+j)=coffmix2(l,l+1+j)&
           -coffmix1(l,l+j)*(dble(l+j)+1.0d0)/(2.0d0*dble(l+j)+1.0d0)
   end do
end do
coffmix2=coffmix2*(Dex/Dso)*(L1+S1+1.0d0)/((J1+2.0d0)*S1)

do l=2,6,2
   do ll=-1,8
      coffmix01(l,ll)=-xi(l)*dble(l*(l+1))/(dble(2*l+1))*coffmix1(l,ll)
      coffmix02(l,ll)=xi(l)*dble(l*(l+1))/(dble(2*l+1))*coffmix2(l,ll)
   end do
end do
coffmix01=coffmix01*Dex/Dso
coffmix02=coffmix02*Dex/Dso



! theta dependence for analytical free energy
open(282,file=""//trim(adjustl(sys))//"_"//trim(adjustl(method))//"_thetadept_phi0_single_eng_&
     T"//trim(adjustl(Tch))//"_lowstJ_rd"//trim(adjustl(rdch))//".txt")
open(283,file=""//trim(adjustl(sys))//"_"//trim(adjustl(method))//"_thetadept_phi0_single_eng_&
     T"//trim(adjustl(Tch))//"_mixexcf1_rd"//trim(adjustl(rdch))//".txt")
open(284,file=""//trim(adjustl(sys))//"_"//trim(adjustl(method))//"_thetadept_phi0_single_eng_&
     T"//trim(adjustl(Tch))//"_mixexcf2_rd"//trim(adjustl(rdch))//".txt")

open(297,file=""//trim(adjustl(sys))//"_"//trim(adjustl(method))//"_"//trim(adjustl(wwoK1Fe))//&
     "_thetadept_phi0_tot_eng_T"//trim(adjustl(Tch))//"_lowstJ_rd"//trim(adjustl(rdch))//".txt")
open(298,file=""//trim(adjustl(sys))//"_"//trim(adjustl(method))//"_"//trim(adjustl(wwoK1Fe))//&
     "_thetadept_phi0_tot_eng_T"//trim(adjustl(Tch))//"_mixexcf1_rd"//trim(adjustl(rdch))//".txt")
open(299,file=""//trim(adjustl(sys))//"_"//trim(adjustl(method))//"_"//trim(adjustl(wwoK1Fe))//&
     "_thetadept_phi0_tot_eng_T"//trim(adjustl(Tch))//"_mixexcf2_rd"//trim(adjustl(rdch))//".txt")
phi=0.0d0
do it=0,int(0.5d0*pi/dtheta)
   theta=dtheta0*dble(it)
   coffmag=0.0d0
   ! m=0
   coffmag(2,0)=(3.0d0*dcos(theta)**2-1.0d0)
   coffmag(4,0)=(35.0d0*dcos(theta)**4-30.0d0*dcos(theta)**2+3.0d0) ! error in second term 30*cos(theta)
   coffmag(6,0)=(231.0d0*dcos(theta)**6-315.0d0*dcos(theta)**4+105.d0*dcos(theta)**2-5.0d0)
   ! m=4
   coffmag(4,4)=dcos(4.0d0*phi)*dsin(theta)**4
   coffmag(6,4)=dcos(4.0d0*phi)*dsin(theta)**4*(11.0d0*dcos(theta)**2-1.0d0)

   engcf=0.0d0
   engmixcf1=0.0d0
   engmixcf2=0.0d0
   do ir=1,nr
      do l=2,6,2
         do m=-l,l
            engcf(ir)=engcf(ir)+coffmag(l,m)*Alm(l,m,ir)*Stevens(l)*extBJ(l,ir)
            abc=xi(l)*dble(l*(l+1))/(dble(2*l+1))*Dex/Dso*coffmag(l,m)*Alm(l,m,1)
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
   freeFe=KuFe*dsin(theta)**2 ! Helmholtz free energy

   if(it.eq.0) then ! set energy level
      ii=1
      free0(1)=freeR(1)*dmult+freeFe
      free0(2)=freeR(2)*dmult+freeFe
      free0(3)=freeR(3)*dmult+freeFe
      freeR0(1)=freeR(1)
      freeR0(2)=freeR(2)
      freeR0(3)=freeR(3)
   end if

   write(282,*) theta,(freeR(1)-freeR0(1))/dble(nr)
   write(283,*) theta,(freeR(2)-freeR0(2))/dble(nr)
   write(284,*) theta,(freeR(3)-freeR0(3))/dble(nr)
   write(297,*) theta,(freeR(1)*dmult+freeFe-free0(1))*dcoff
   write(298,*) theta,(freeR(2)*dmult+freeFe-free0(2))*dcoff
   write(299,*) theta,(freeR(3)*dmult+freeFe-free0(3))*dcoff
end do ! for it

! phi dependence for analytical free energy
open(382,file=""//trim(adjustl(sys))//"_"//trim(adjustl(method))//"_phidept_theta0.5pi_single_eng_&
     T"//trim(adjustl(Tch))//"_lowstJ_rd"//trim(adjustl(rdch))//".txt")
open(383,file=""//trim(adjustl(sys))//"_"//trim(adjustl(method))//"_phidept_theta0.5pi_single_eng_&
     T"//trim(adjustl(Tch))//"_mixexcf1_rd"//trim(adjustl(rdch))//".txt")
open(384,file=""//trim(adjustl(sys))//"_"//trim(adjustl(method))//"_phidept_theta0.5pi_single_eng_&
     T"//trim(adjustl(Tch))//"_mixexcf2_rd"//trim(adjustl(rdch))//".txt")
open(397,file=""//trim(adjustl(sys))//"_"//trim(adjustl(method))//"_"//trim(adjustl(wwoK1Fe))//&
     "_phidept_theta0.5pi_tot_eng_T"//trim(adjustl(Tch))//"_lowstJ_rd"//trim(adjustl(rdch))//".txt")
open(398,file=""//trim(adjustl(sys))//"_"//trim(adjustl(method))//"_"//trim(adjustl(wwoK1Fe))//&
     "_phidept_theta0.5pi_tot_eng_T"//trim(adjustl(Tch))//"_mixexcf1_rd"//trim(adjustl(rdch))//".txt")
open(399,file=""//trim(adjustl(sys))//"_"//trim(adjustl(method))//"_"//trim(adjustl(wwoK1Fe))//&
     "_phidept_theta0.5pi_tot_eng_T"//trim(adjustl(Tch))//"_mixexcf2_rd"//trim(adjustl(rdch))//".txt")
theta=pi/2.0d0
do ip=0,100
   phi=pi/2.0d0*dble(ip)/100.0d0
   coffmag=0.0d0
   ! m=0
   coffmag(2,0)=(3.0d0*dcos(theta)**2-1.0d0)
   coffmag(4,0)=(35.0d0*dcos(theta)**4-30.0d0*dcos(theta)**2+3.0d0) ! error in second term 30*cos(theta)
   coffmag(6,0)=(231.0d0*dcos(theta)**6-315.0d0*dcos(theta)**4+105.d0*dcos(theta)**2-5.0d0)
   ! m=4
   coffmag(4,4)=dcos(4.0d0*phi)*dsin(theta)**4
   coffmag(6,4)=dcos(4.0d0*phi)*dsin(theta)**4*(11.0d0*dcos(theta)**2-1.0d0)

   engcf=0.0d0
   engmixcf1=0.0d0
   engmixcf2=0.0d0
   do ir=1,nr
      do l=2,6,2
         do m=-l,l
            engcf(ir)=engcf(ir)+coffmag(l,m)*Alm(l,m,ir)*Stevens(l)*extBJ(l,ir)
            abc=xi(l)*dble(l*(l+1))/(dble(2*l+1))*Dex/Dso*coffmag(l,m)*Alm(l,m,1)
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
   freeFe=KuFe*dsin(theta)**2 ! Helmholtz free energy

   write(382,*) phi,(freeR(1)-freeR0(1))/dble(nr)
   write(383,*) phi,(freeR(2)-freeR0(2))/dble(nr)
   write(384,*) phi,(freeR(3)-freeR0(3))/dble(nr)
   
   write(397,*) phi,(freeR(1)*dmult+freeFe-free0(1))*dcoff
   write(398,*) phi,(freeR(2)*dmult+freeFe-free0(2))*dcoff
   write(399,*) phi,(freeR(3)*dmult+freeFe-free0(3))*def
   
end do ! for ip

end program eng_level



