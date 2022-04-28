program phidept 
implicit none
integer(4),parameter :: nrmax=30,itempmax=1000,nJexmax=10
real(8),parameter :: pi=dacos(-1.0d0)
real(8),parameter :: kb=1.3806488d-23 ! (J/K)
real(8),parameter :: muB=9.2740100783d-24 ! (Merg/kOe) or (kerg/Oe) or (J/T)
real(8),parameter :: mu0=4.0d0*pi*1.0d-7 ! (NA^-2) 
real(8),parameter :: joule=10.0d0 ! (Merg/J)

integer(4) :: nJex,iHmag,nfe,nr,n4f,itemp,ncfp,lm,lp
real(8) :: J0,S1,L1,J1,lambda,rl(6),Hmax,dH,Tmin,Tmax,dT,p,f(8),f2(8),WJ(0:nJexmax)
real(8) :: ax,cx,Femag,Femag0,Rmag0,Tc,s,mFe0,MsT_Kuzmin,Jz2exp,Jz4exp,Ms0,Tcomp
real(8) :: JM(100,2),SO(100),Stevens(6),xi(6),Clm(6,-6:6,nrmax),Blm(6,-6:6,nrmax),g(100),cclm(6,-6:6)
real(8) :: extBJ(-1:8,nrmax),extLJ(-1:8,nrmax),CC(7,-7:7,nrmax),disp(7,nrmax),rdcoffL,rdcoffS,rdcoff(0:6)
real(8) :: Hex0(nrmax),Hmag,thetamag,phimag,thetatot,theta,phi,dcoff,coffT(6),extBJ0(-1:8,nrmax)
real(8) :: abc,def,ghi,jkl,theta0,phi0,freeloc(nrmax),freeloc0(0:4,0:100,nrmax),freet0loc(nrmax)
real(8) :: part(nrmax),part0(nrmax),temp,free,free0,dtheta,dtheta0,dl,dlm,dlp
real(8) :: stspin0(3,nrmax),stspinS0(3,nrmax),stspinL0(3,nrmax),Alm(6,-6:6,nrmax)
real(8) :: absj0(nrmax),abss0(nrmax),absl0(nrmax),thetaj0(nrmax),phij0(nrmax)
real(8) :: thetaS0(nrmax),phiS0(nrmax),thetaL0(nrmax),phiL0(nrmax)

real(8) :: mFeT(itempmax),HexT(nrmax,0:itempmax),K1FeT(itempmax),KuFe
real(8) :: Kloc(11,nrmax),dK1(2),dK2(3),dK3(4),dKi1st(9),dKimix1(9),dKimix21(9),dKi2nd(9),dKidisp(9)
real(8) :: dmag1st(2,nrmax),dmagmixcf1(2,nrmax),dmagmixcf2(2,nrmax),dmagmixex1(2,nrmax),dmagmixex2(2,nrmax)
real(8) :: dmagcl1st(2,nrmax),dmagclmixcf1(2,nrmax),dmagclmixcf2(2,nrmax),dmagclmixex1(2,nrmax),dmagclmixex2(2,nrmax)
real(8) :: dclKi(9),dclKimix1(9),dclKimix2(9)
real(8) :: coffKi(9,6,-6:6),coffmix1(6,-1:8),coffmix2(6,-1:8),coffKi2nd(9,6,0:6),coffdisp(9,6),coffmag(9)
real(8) :: Tmix1(6,nrmax),Tmix2(6,nrmax),T0mix1(6,nrmax),T0mix2(6,nrmax),Tclmix1(6,nrmax),Tclmix2(6,nrmax)
real(8) :: d3j(6,nrmax),d3j2(6,nrmax),t3j(6,nrmax),q3j(6,nrmax)

integer(4) :: m,n,i,j,k,l,imax,ii,jj,INFO,it,ip,ll,ir,LWORK
character(50) :: sys,direction,nJexch,Hch,ratom,tatom,method,sch,wwoK1Fe,rdch
real(8) :: aleng,bleng,cleng,gmm,volume,dmult,rd
real(8) :: K1sub(0:4,nrmax),K2sub(0:4,nrmax),K3sub(0:4,nrmax),K1Fe,K2Fe,K3Fe,rK2Fe,rK3Fe
!real(8) :: K1Fe=0.77d0,K2Fe=1.21d0,K3Fe=0.11d0
real(8) :: ri1,ri2,rk1,rk2,rip,rkp,Brillouin,Langevin,dnu,Dex,Dso,w6jsym,xfp

real(8),allocatable :: matrix(:,:),eigen(:),eigen0(:),RWORK(:),Ust(:,:,:,:),Vst(:,:,:)
complex(kind(0d0)),allocatable :: ZUst(:,:,:,:),WORK(:),zunit0(:,:,:),zunit(:,:,:),zH(:,:,:),zH0(:,:,:),zspin(:,:,:),zorbit(:,:,:),zmatrix(:,:)

! set paremeter
read(5,*) sys
write(6,*) "system=",sys
read(5,*) nJex ! 0:lowestJ, 1:LowestJ+secondJ
write(6,*) "number of excited states",nJex
read(5,*) method
write(6,*) "method=",method
read(5,*) wwoK1Fe
write(6,*) "wwoK1Fe",wwoK1Fe
read(5,*) Tmin,Tmax,dT
write(6,*) "temperature",Tmin,Tmax,dT,"(K)"
read(5,*) dtheta
write(6,*) "increment of \theta",dtheta,"(rad)",dtheta*180.0d0/pi,"(deg.)"
write(nJexch,"(I2)") nJex
nJexch=adjustl(nJexch)

! input files
open(33,file="~/research/CFPs/Sm_lambda411K/cfp"//trim(sys)//"_"//trim(method)//".txt")
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
rdch=adjustl(rdch)
volume=aleng*bleng*cleng*dsin(pi*gmm/180.0d0)*1.0d-30
dcoff=kb/volume*1.0d-6 ! convert from [K] to [MJ/m^3]
write(6,*) "Fe sublattice magnetization @T=0 [T]"
write(6,*) mu0*mFe0*muB/volume

if(nr.gt.nrmax) then
   write(6,*) "error: nr > nrmax"
   stop
end if

if(sys.eq."SmFe12_L3") then
   call rareearth_SmL3(ratom,n4f,nJex,imax,JM,S1,L1,J1)
else 
   call rareearth(ratom,n4f,nJex,imax,JM,S1,L1,J1)
end if

! spin orbit coupling
do jj=1,imax
   SO(jj)=0.5d0*lambda*(JM(jj,1)*(JM(jj,1)+1.0d0)-L1*(L1+1.0d0)-S1*(S1+1.0d0))
end do

! lande g factor for excited J multiplet
if(imax.gt.100) then
   write(6,*) "imax too large"
   stop
end if
do i=1,imax
   g(i)=1.0d0+(JM(i,1)*(JM(i,1)+1.0d0)+S1*(S1+1.0d0)-L1*(L1+1.0d0))/(2.0d0*JM(i,1)*(JM(i,1)+1.0d0))
end do

allocate(matrix(imax,imax),eigen(imax),eigen0(imax),RWORK(3*imax-2),Ust(imax,imax,6,-6:6),Vst(imax,imax,-1:1))
allocate(zH0(imax,imax,nr),zH(imax,imax,nr),ZUst(imax,imax,6,-6:6),WORK(2*imax-1),zspin(imax,imax,3),zorbit(imax,imax,3),&
     zunit0(imax,imax,nr),zunit(imax,imax,nr),zmatrix(imax,imax))

! matrix element of spherical tensor operators
if(sys.eq."SmFe12_L3") then
   call rdmatrix_SmL3(n4f,rdcoff,L1)
else
   call rdmatrix(n4f,rdcoff,L1) ! reduced matrix element
end if
call stfactor(Stevens,L1,J1,S1,rdcoff) ! Stevens factor
call xifactor(xi,L1,J1,S1,rdcoff) ! Stevens factor
call numfactor(cclm)
write(6,*) "Stevens check "
do ll=2,6,2
   write(6,*) Stevens(ll),xi(ll)
end do

! crystal field parameters
Blm=0.0d0
Clm=0.0d0
do ir=1,nr
   do l=2,6,2
      do m=-l,l
         def=1.0d0/cclm(l,m)*dsqrt((2.0d0*dble(l)+1.0d0)/(4.0d0*pi))*rdcoff(l)
         if(m.ne.0) then
            def=def/dsqrt(2.0d0)
         end if
         Blm(l,m,ir)=Alm(l,m,ir)*Stevens(l)
         Clm(l,m,ir)=Alm(l,m,ir)*def
      end do
   end do
end do
call CkqL(L1,S1,J1,Ust,imax,nJex) ! matrix element of tensor operator
call CkqS(L1,S1,J1,Vst,imax,nJex) ! matrix element of vector operator
call OkqL(ZUst,Ust,imax)
call moment_op(imax,L1,S1,zspin,zorbit,JM,Vst,Ust) ! matrix element of magnetic moment

! Hamiltonian
zH0=0.0d0
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


! saturated magnetization, molecular field, anisotropy constants for Fe sublattice
if(int((Tmax-Tmin)/dT).gt.itempmax) then
   write(6,*) "error: itempmax over the number of T plot"
end if

! in the case of Ce2Fe14B, K1Fe=0.77d0 [MJ/m^3], K2Fe=1.21 [MJ/m^3], K3Fe= 0.11 [MJ/m^3]
!write(6,*) K1Fe/(aleng**2*cleng*1.0d-30)*kb*1.0d-6,K1Fe*rK2Fe/(aleng**2*cleng*1.0d-30)*kb*1.0d-6,K1Fe*rK3Fe/(aleng**2*cleng*1.0d-30)*kb*1.0d-6
!stop

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
!   K1FeT(itemp)=K1Fe*abc**3+8.0d0/7.0d0*K2Fe*(abc**3-abc**10)+8.0d0/7.0d0*K3Fe*(abc**3-18.0d0/11.0d0*abc**10+7.0d0/11.0d0*abc**21)
   if(sys.eq."Sm2Fe17N3tri") K1FeT(itemp)=0.0d0
   write(287,*) temp,mu0*mFeT(itemp)*muB/volume ! in unit Tesla
   write(288,*) temp,K1FeT(itemp)/(aleng**2*cleng*1.0d-30)*kb*1.0d-6 ! [MJ/m^3]
   do ir=1,nr
      HexT(ir,itemp)=-Hex0(ir)*abc
   end do
end do

Hmag=0.0d0
thetamag=0.0d0
phimag=0.0d0
iHmag=0
Hmag=dble(iHmag)
zH=(0.0d0,0.0d0)
do ir=1,nr
   do jj=1,imax
      do ii=1,imax
         zH(ii,jj,ir)=zH0(ii,jj,ir)+& ! Zeeman term
              ((2.0d0*zspin(ii,jj,3)+zorbit(ii,jj,3))*dcos(thetamag)&
              +(2.0d0*zspin(ii,jj,1)+zorbit(ii,jj,1))*dcos(phimag)*dsin(thetamag)&
              +(2.0d0*zspin(ii,jj,2)+zorbit(ii,jj,2))*dsin(phimag)*dsin(thetamag))*Hmag
      end do
   end do
end do

! temperature dependence
itemp=0  ! loop for Temperature 
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


! angular dpendent free energy
free0=0.0d0
part0=0.0d0
theta0=0.0d0
phi0=0.0d0
dtheta0=0.5d0*pi/dble(int(0.5d0*pi/dtheta))
open(75,file=""//trim(sys)//"_phi_theta0.txt")
open(76,file=""//trim(sys)//"_phi_zoom_theta0.txt")
open(77,file=""//trim(sys)//"_totangle_theta0.txt")
open(78,file=""//trim(sys)//"_SLangle_theta0.txt")
!do ip=0,40 ! \phi=0
!   phi=pi*dble(ip)/40.0d0
do ip=0,100 ! \phi=pi*dble(ip)/8.0d0
   phi=pi*dble(ip)/100.0d0
   if(ip.eq.3) phi=pi/100.0d0 ! for K3 tilda
   do it=0,1 ! take three point
      theta=dtheta0*dble(it)
      part=0.0d0
      free=0.0d0
      do ir=1,nr
         zmatrix=(0.0d0,0.0d0)
         do jj=1,imax
            do ii=1,imax
               zmatrix(ii,jj)=zH(ii,jj,ir)&
                    +(zspin(ii,jj,1)*dcos(phi)*dsin(theta)&
                    +zspin(ii,jj,2)*dsin(phi)*dsin(theta)&
                    +zspin(ii,jj,3)*dcos(theta))*2.0d0*HexT(ir,itemp)
            end do
         end do
         
         LWORK=2*imax-1
         call zheev('V','L',imax,zmatrix,imax,eigen,WORK,LWORK,RWORK,INFO)
         freeloc(ir)=eigen(1) 
         if(it.eq.0) freet0loc(ir)=freeloc(ir)
         if(it.le.3.and.ip.le.100)  freeloc0(it,ip,ir)=freeloc(ir)-freet0loc(ir) 
         free=free+freeloc(ir)-freet0loc(ir)
         do jj=1,imax
            do ii=1,imax
               zunit(ii,jj,ir)=zmatrix(ii,jj)
            end do
         end do
      end do
      if(it.eq.1) then
         write(75,"(4f17.10)") &
           phi,freeloc0(it,ip,1),freeloc0(it,ip,2),freeloc0(it,ip,1)+freeloc0(it,ip,2)
         write(76,"(3f17.10)") phi,free/4.0d0,(freeloc0(it,ip,1)+freeloc0(it,ip,2))/2.0d0


         stspin0=0.0d0
         stspinS0=0.0d0
         stspinL0=0.0d0
         do ir=1,nr
            do ll=1,3
               m=1
               def=0.0d0
               ghi=0.0d0
               do jj=1,imax
                  do ii=1,imax
                     def=def-2.0d0*dreal(conjg(zunit(ii,m,ir))*zspin(ii,jj,ll)*zunit(jj,m,ir))
                     ghi=ghi-dreal(conjg(zunit(ii,m,ir))*zorbit(ii,jj,ll)*zunit(jj,m,ir))
                  end do
               end do
               stspin0(ll,ir)=def+ghi
               stspinS0(ll,ir)=def
               stspinL0(ll,ir)=ghi
            end do ! for spin component
         end do ! for ir=1,2
         do ir=1,nr
            call vector_to_angle(stspin0(1,ir),thetaj0(ir),phij0(ir),absj0(ir))
            call vector_to_angle(stspinS0(1,ir),thetaS0(ir),phiS0(ir),absS0(ir))
            call vector_to_angle(stspinL0(1,ir),thetaL0(ir),phiL0(ir),absL0(ir))
         end do
         write(77,"(3f17.10)") phi,phij0(1),phij0(2)
         write(78,"(5f17.10)") phi,phiS0(1),phiS0(2),phiL0(1),phiL0(2)
      end if ! if it.eq.1
   end do ! for \theta it
end do ! for \phi ip



free0=0.0d0
part0=0.0d0
theta0=0.0d0
phi0=0.0d0
dtheta0=0.5d0*pi/dble(int(0.5d0*pi/dtheta))
open(85,file=""//trim(sys)//"_phi_theta90.txt")
open(86,file=""//trim(sys)//"_phi_zoom_theta90.txt")
open(87,file=""//trim(sys)//"_totangle_theta90.txt")
open(88,file=""//trim(sys)//"_SLangle_theta90.txt")
!do ip=0,40 ! \phi=0
!   phi=pi*dble(ip)/40.0d0
do ip=0,100 ! \phi=pi*dble(ip)/8.0d0
   phi=pi*dble(ip)/100.0d0
   if(ip.eq.3) phi=pi/100.0d0 ! for K3 tilda
   do it=0,1 ! take three point
      theta=0.5d0*pi*dble(it)
      part=0.0d0
      free=0.0d0
      do ir=1,nr
         zmatrix=(0.0d0,0.0d0)
         do jj=1,imax
            do ii=1,imax
               zmatrix(ii,jj)=zH(ii,jj,ir)&
                    +(zspin(ii,jj,1)*dcos(phi)*dsin(theta)&
                    +zspin(ii,jj,2)*dsin(phi)*dsin(theta)&
                    +zspin(ii,jj,3)*dcos(theta))*2.0d0*HexT(ir,itemp)
            end do
         end do
         
         LWORK=2*imax-1
         call zheev('V','L',imax,zmatrix,imax,eigen,WORK,LWORK,RWORK,INFO)
         freeloc(ir)=eigen(1) 
         if(it.eq.0) freet0loc(ir)=freeloc(ir)
         if(it.le.3.and.ip.le.100)  freeloc0(it,ip,ir)=freeloc(ir)-freet0loc(ir) 
         free=free+freeloc(ir)-freet0loc(ir)
         do jj=1,imax
            do ii=1,imax
               zunit(ii,jj,ir)=zmatrix(ii,jj)
            end do
         end do
      end do
      if(it.eq.1) then
         write(85,"(4f17.10)") phi,freeloc0(it,ip,1),freeloc0(it,ip,2),freeloc0(it,ip,1)+freeloc0(it,ip,2)
         write(86,"(3f17.10)") phi,free/4.0d0,(freeloc0(it,ip,1)+freeloc0(it,ip,2))/2.0d0
!         write(6,*) phi,free/4.0d0,(freeloc0(it,ip,1)+freeloc0(it,ip,2))/2.0d0


         stspin0=0.0d0
         stspinS0=0.0d0
         stspinL0=0.0d0
         do ir=1,nr
            do ll=1,3
               m=1
               def=0.0d0
               ghi=0.0d0
               do jj=1,imax
                  do ii=1,imax
                     def=def-2.0d0*dreal(conjg(zunit(ii,m,ir))*zspin(ii,jj,ll)*zunit(jj,m,ir))
                     ghi=ghi-dreal(conjg(zunit(ii,m,ir))*zorbit(ii,jj,ll)*zunit(jj,m,ir))
                  end do
               end do
               stspin0(ll,ir)=def+ghi
               stspinS0(ll,ir)=def
               stspinL0(ll,ir)=ghi
            end do ! for spin component
         end do ! for ir=1,2
         do ir=1,nr
            call vector_to_angle(stspin0(1,ir),thetaj0(ir),phij0(ir),absj0(ir))
            call vector_to_angle(stspinS0(1,ir),thetaS0(ir),phiS0(ir),absS0(ir))
            call vector_to_angle(stspinL0(1,ir),thetaL0(ir),phiL0(ir),absL0(ir))
         end do
         write(87,*) phi,phij0(1),phij0(2)
         write(88,"(5f17.10)") phi,phiS0(1),phiS0(2),phiL0(1),phiL0(2)
      end if ! if it.eq.1

   end do ! for \theta it
end do ! for \phi ip

   
end program phidept


