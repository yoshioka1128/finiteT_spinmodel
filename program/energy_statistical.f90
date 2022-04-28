program energy_statistical
implicit none
integer(4),parameter :: nrmax=30,nJexmax=10,mm=100
real(8),parameter :: pi=dacos(-1.0d0)
real(8),parameter :: kb=1.3806488d-23 ! (J/K)
integer(4) :: nJex,nfe,nr,n4f,ncfp,lm,lp,m,n,i,j,k,l,imax,ii,jj,INFO,it,ip,ir,ir0,LWORK,itemp
real(8) :: S1,L1,J1,lambda,rl(6),s,p,Tc,MTM0,MsT_Kuzmin,Hex0(nrmax),Hex(nrmax),thetamag,thetatot,theta,phi,dcoff
real(8) :: JM(mm,2),SO(mm),Clm(6,-6:6,nrmax),cclm(6,-6:6),rdcoff(0:6)
real(8) :: abc,def,freeloc(nrmax),freet0(nrmax),part(nrmax),temp,dtheta,dtheta0
real(8) :: stspin0(3,nrmax),stspinS0(3,nrmax),stspinL0(3,nrmax),Alm(6,-6:6,nrmax),gibbs
real(8) :: absj0(nrmax),abss0(nrmax),absl0(nrmax),thetaj0(nrmax),phij0(nrmax),WJ(0:nJexmax)
real(8) :: thetaS0(nrmax),phiS0(nrmax),thetaL0(nrmax),phiL0(nrmax),KuTM,MTM
real(8) :: alc,blc,clc,gmm,volume,dmult,rd,K1TM,rK2TM,rK3TM
character(50) :: sys,direction,nJexch,Hch,ratom,tatom,method,sch,wwoK1TM,rdch,Tch
real(8),allocatable :: eigen(:),RWORK(:),Ust(:,:,:,:),Vst(:,:,:)
complex(kind(0d0)),allocatable :: ZUst(:,:,:,:),WORK(:),zH(:,:,:),zH0(:,:,:),zspin(:,:,:),zorbit(:,:,:)

! set paremeter
read(5,*) sys
read(5,*) temp
read(5,*) nJex ! 0:lowestJ, 1:LowestJ+secondJ
read(5,*) method
read(5,*) wwoK1TM
read(5,*) dtheta
write(6,"(2A)") "system=",sys
write(6,"(A,f10.5,A)") "T=",temp," K"
write(6,"(A,I0)") "number of excited states: " ,nJex
write(6,"(2A)") "method=",method
write(6,"(2A)") "Ku^TM=",wwoK1TM
write(6,"(A,f10.5,A,f10.5,A)") "increment of \theta",dtheta," rad.",dtheta*180.0d0/pi," deg."
write(nJexch,"(I2)") nJex
write(Tch,"(f5.1)") temp
Tch=adjustl(Tch)
dtheta0=0.5d0*pi/dble(int(0.5d0*pi/dtheta))
nJexch=adjustl(nJexch)

! input files
open(33,file="~/research/CFPs/Sm_lambda411K/cfp"//trim(sys)//"_"//trim(method)//".txt")
read(33,*) ratom,nr,tatom,nfe
read(33,*) dmult
read(33,*) alc,blc,clc,gmm
read(33,*) lambda
read(33,*) MTM0,Tc,s,p
read(33,*) K1TM,rK2TM,rK3TM
write(6,"(A,A3,I3,A3,I3)") "system: ",ratom,nr,tatom,nfe
write(6,"(A,f15.5)") "multiplicity: ",dmult
write(6,"(A,3f15.5,A)") "(a,b,c)=",alc,blc,clc," ang."
write(6,"(A,f15.5,A)") "lambda=",lambda," K"
write(6,"(A)") "Ms [muB/unitcell], Tc [K], s, p for Kuzmin fitting"
write(6,"(4f15.5)") MTM0, Tc, s, p ! MTM0 [muB per #nfe]
write(6,"(A)") "K1^TM [K/unitcell],K2TM/K1TM,K3TM/K1TM"
write(6,"(3f15.5)") K1TM,rK2TM,rK3TM
write(sch,"(f6.3)") s
Alm=0.0d0
do ir=1,nr
   write(6,"(A,I3)") "atom",ir
   
   read(33,*)
   read(33,*) Hex0(ir),rd
   write(6,"(A,f10.5,A,f10.5)") " Bex(0)=",Hex0(ir),", reduced factor: ",rd
   write(rdch,"(f5.1)") rd
   Hex0(ir)=rd*Hex0(ir)
   read(33,*) rl(2),rl(4),rl(6),ncfp
   write(6,"(A,I3)") " number of LM: ",ncfp
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
volume=alc*blc*clc*dsin(pi*gmm/180.0d0)*1.0d-30
dcoff=kb/volume*1.0d-6 ! convert from [K] to [MJ/m^3]
rdch=adjustl(rdch)

! R atom
if(sys.eq."SmFe12_L3") then
   call rareearth_SmL3(ratom,n4f,nJex,imax,JM,S1,L1,J1)
else 
   call rareearth(ratom,n4f,nJex,imax,JM,S1,L1,J1)
end if
! TM sublattice
abc=MsT_Kuzmin(temp,MTM0,Tc,s,p)/MTM0
do ir=1,nr
   Hex(ir)=-Hex0(ir)*abc
end do
def=abc**3+&
     8.0d0/7.0d0*rK2TM*(abc**3-abc**10)+&
     8.0d0/7.0d0*rK3TM*(abc**3-18.0d0/11.0d0*abc**10+7.0d0/11.0d0*abc**21)
if(wwoK1TM.eq."wK1Fe") then
   KuTM=K1TM*def ! [K/unitcell]
else if(wwoK1TM.eq."woK1TM") then
   KuTM=0.0d0
else
   write(6,*) "read wwoK1TM error",wwoK1TM
   stop
end if
write(6,"(A,f25.15,A)") "K1^TM(T)=",KuTM," K/unitcell"

! analytical calculations
call energy_analytical(temp,sys,method,volume,lambda,Alm,KuTM,Hex,ratom,nr,dmult,rdch)

! statistical calculations
allocate(eigen(imax),RWORK(3*imax-2),Ust(imax,imax,6,-6:6),Vst(imax,imax,-1:1))
allocate(zH0(imax,imax,nr),zH(imax,imax,nr),ZUst(imax,imax,6,-6:6),WORK(2*imax-1),zspin(imax,imax,3),zorbit(imax,imax,3))
LWORK=2*imax-1

! spin orbit coupling
do i=1,imax
   SO(i)=0.5d0*lambda*(JM(i,1)*(JM(i,1)+1.0d0)-L1*(L1+1.0d0)-S1*(S1+1.0d0))
end do

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

! hierachy of energy scale
ir0=1
zH0=0.0d0
do ir=1,nr
   do ii=1,imax
      zH0(ii,ii,ir)=zH0(ii,ii,ir)+SO(ii) ! spin orbit term
   end do
end do
open(270,file=""//trim(sys)//"_"//trim(method)//"_englevel_T"//trim(Tch)//"_diag_so_nJex"//trim(nJexch)//"_rd"//trim(rdch)//".txt")
do ii=1,imax
   write(270,*) ii,real(zH0(ii,ii,ir0))-SO(1),0 ! Hso only
   write(270,*) ii,real(zH0(ii,ii,ir0))-SO(1),1
   write(270,*)
end do
! spin-orbit and exchange term
open(271,file=""//trim(sys)//"_"//trim(method)//"_englevel_T"//trim(Tch)//"_diag_soex_nJex"//trim(nJexch)//"_rd"//trim(rdch)//".txt")
dtheta0=0.5d0*pi/dble(int(0.5d0*pi/dtheta))
call Hexch(0.0d0,0.0d0,imax,zH(:,:,ir0),zH0(:,:,ir0),zspin,Hex(ir0))
call zheev('V','L',imax,zH(:,:,ir0),imax,eigen,WORK,LWORK,RWORK,INFO)
do ii=1,imax
   write(271,*) ii,eigen(ii)-SO(1),0
   write(271,*) ii,eigen(ii)-SO(1),1
   write(271,*)
end do
! spin-orbit, exchange, and CF term
do jj=1,imax
   do ii=1,imax
      abc=0.0d0
      do l=2,6,2
         do m=-l,l
            abc=abc+ZUst(ii,jj,l,m)*Clm(l,m,ir0) ! crystal field term
         end do
      end do
      zH0(ii,jj,ir0)=zH0(ii,jj,ir0)+abc
   end do
end do
open(272,file=""//trim(sys)//"_"//trim(method)//"_englevel_T"//trim(Tch)//"_diag_soexcf_nJex"//trim(nJexch)//"_rd"//trim(rdch)//".txt")
call Hexch(theta,phi,imax,zH(:,:,ir0),zH0(:,:,ir0),zspin,Hex(ir0))
call zheev('V','L',imax,zH(:,:,ir0),imax,eigen,WORK,LWORK,RWORK,INFO)      
do ii=1,imax
   write(272,*) ii,eigen(ii)-SO(1),0 ! Hso+Hex+Hcf
   write(272,*) ii,eigen(ii)-SO(1),1
   write(272,*)
end do

! SO and CF and Hamiltonian
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
! free energy at theta=0 
part=0.0d0
do ir=1,nr
   call Hexch(0.0d0,0.0d0,imax,zH(:,:,ir),zH0(:,:,ir),zspin,Hex(ir))
   call zheev('V','L',imax,zH(:,:,ir),imax,eigen,WORK,LWORK,RWORK,INFO)
   if(temp.eq.0.0d0) then
      freeloc(ir)=eigen(1) 
   else
      do  ii=1,imax
         part(ir)=part(ir)+dexp(-(eigen(ii)-eigen(1))/temp)
      end do
      freeloc(ir)=eigen(1)-temp*dlog(part(ir))
   end if
   freet0(ir)=freeloc(ir)
end do


! angular dpendent free energy
open(95,file=""//trim(sys)//"_"//trim(method)//"_T"//trim(Tch)//"_eng_theta_phi22.5_nJex"//trim(nJexch)//"_rd"//trim(rdch)//".txt")
open(98,file=""//trim(sys)//"_"//trim(method)//"_T"//trim(Tch)//"_SLangle_theta_phi22.5_nJex"//trim(nJexch)//"_rd"//trim(rdch)//".txt")
do it=0,int(0.5d0*pi/dtheta) ! take three point
   theta=dtheta0*dble(it)
   part=0.0d0
   abc=0.0d0
   do ir=1,nr
      call Hexch(theta,pi/8.0d0,imax,zH(:,:,ir),zH0(:,:,ir),zspin,Hex(ir))
      call zheev('V','L',imax,zH(:,:,ir),imax,eigen,WORK,LWORK,RWORK,INFO)
      if(temp.eq.0.0d0) then
         freeloc(ir)=eigen(1) 
      else
         do  ii=1,imax
            part(ir)=part(ir)+dexp(-(eigen(ii)-eigen(1))/temp)
         end do
         freeloc(ir)=eigen(1)-temp*dlog(part(ir))
      end if
      freeloc(ir)=freeloc(ir)-freet0(ir)
      abc=abc+freeloc(ir)
      call momentR_exp(temp,imax,zH(:,:,ir),zspin,zorbit,part(ir),eigen,stspin0(:,ir),stspinS0(:,ir),stspinL0(:,ir))
      call vector_to_angle(stspin0(:,ir),thetaj0(ir),phij0(ir),absj0(ir))
      call vector_to_angle(stspinS0(:,ir),thetaS0(ir),phiS0(ir),absS0(ir))
      call vector_to_angle(stspinL0(:,ir),thetaL0(ir),phiL0(ir),absL0(ir))
   end do ! do ir
   gibbs=abc*dmult+KuTM*dsin(theta)**2
   write(95,"(4f17.10)") theta,freeloc(1),freeloc(2),gibbs*dcoff
   write(98,"(9f17.10)") theta,thetaS0(1),thetaS0(2),thetaL0(1),thetaL0(2),phiS0(1),phiS0(2),phiL0(1),phiL0(2)
end do ! for \theta it

open(105,file=""//trim(sys)//"_"//trim(method)//"_T"//trim(Tch)//"_eng_theta_phi0_nJex"//trim(nJexch)//"_rd"//trim(rdch)//".txt")
open(108,file=""//trim(sys)//"_"//trim(method)//"_T"//trim(Tch)//"_SLangle_theta_phi0_nJex"//trim(nJexch)//"_rd"//trim(rdch)//".txt")
do it=0,int(0.5d0*pi/dtheta) ! take three point
   theta=dtheta0*dble(it)
   part=0.0d0
   abc=0.0d0
   do ir=1,nr
      call Hexch(theta,0.0d0,imax,zH(:,:,ir),zH0(:,:,ir),zspin,Hex(ir))
      call zheev('V','L',imax,zH(:,:,ir),imax,eigen,WORK,LWORK,RWORK,INFO)
      if(temp.eq.0.0d0) then
         freeloc(ir)=eigen(1) 
      else
         do  ii=1,imax
            part(ir)=part(ir)+dexp(-(eigen(ii)-eigen(1))/temp)
         end do
         freeloc(ir)=eigen(1)-temp*dlog(part(ir))
      end if
      freeloc(ir)=freeloc(ir)-freet0(ir) 
      abc=abc+freeloc(ir)
      call momentR_exp(temp,imax,zH(:,:,ir),zspin,zorbit,part(ir),eigen,stspin0(:,ir),stspinS0(:,ir),stspinL0(:,ir))
      call vector_to_angle(stspin0(:,ir),thetaj0(ir),phij0(ir),absj0(ir))
      call vector_to_angle(stspinS0(:,ir),thetaS0(ir),phiS0(ir),absS0(ir))
      call vector_to_angle(stspinL0(:,ir),thetaL0(ir),phiL0(ir),absL0(ir))
   end do
   gibbs=abc*dmult+KuTM*dsin(theta)**2
   write(105,"(4f17.10)") theta,freeloc(1),freeloc(2),gibbs*dcoff
   write(108,"(9f17.10)") theta,thetaS0(1),thetaS0(2),thetaL0(1),thetaL0(2),phiS0(1),phiS0(2),phiL0(1),phiL0(2)
end do ! for \theta it


open(75,file=""//trim(sys)//"_"//trim(method)//"_T"//trim(Tch)//"_eng_phi_theta0_nJex"//trim(nJexch)//"_rd"//trim(rdch)//".txt")
open(78,file=""//trim(sys)//"_"//trim(method)//"_T"//trim(Tch)//"_SLangle_phi_theta0_nJex"//trim(nJexch)//"_rd"//trim(rdch)//".txt")
theta=dtheta0
it=1
do ip=0,100 ! \phi=pi*dble(ip)/8.0d0
   phi=pi*dble(ip)/100.0d0
   if(ip.eq.3) phi=pi/100.0d0 ! for K3 tilda
   part=0.0d0
   abc=0.0d0
   do ir=1,nr
      call Hexch(dtheta0,phi,imax,zH(:,:,ir),zH0(:,:,ir),zspin,Hex(ir))
      call zheev('V','L',imax,zH(:,:,ir),imax,eigen,WORK,LWORK,RWORK,INFO)
      if(temp.eq.0.0d0) then
         freeloc(ir)=eigen(1) 
      else
         do  ii=1,imax
            part(ir)=part(ir)+dexp(-(eigen(ii)-eigen(1))/temp)
         end do
         freeloc(ir)=eigen(1)-temp*dlog(part(ir))
      end if
      freeloc(ir)=freeloc(ir)-freet0(ir) 
      abc=abc+freeloc(ir)
      call momentR_exp(temp,imax,zH(:,:,ir),zspin,zorbit,part(ir),eigen,stspin0(:,ir),stspinS0(:,ir),stspinL0(:,ir))
      call vector_to_angle(stspin0(:,ir),thetaj0(ir),phij0(ir),absj0(ir))
      call vector_to_angle(stspinS0(:,ir),thetaS0(ir),phiS0(ir),absS0(ir))
      call vector_to_angle(stspinL0(:,ir),thetaL0(ir),phiL0(ir),absL0(ir))
   end do ! ir
   gibbs=abc*dmult+KuTM*dsin(theta)**2
   write(75,"(4f17.10)") phi,freeloc(1),freeloc(2),gibbs*dcoff
   write(78,"(9f17.10)") phi,phiS0(1),phiS0(2),phiL0(1),phiL0(2),thetaS0(1),thetaS0(1),thetaL0(2),thetaL0(2)
end do ! for \phi ip

open(85,file=""//trim(sys)//"_"//trim(method)//"_T"//trim(Tch)//"_eng_phi_theta90_nJex"//trim(nJexch)//"_rd"//trim(rdch)//".txt")
open(88,file=""//trim(sys)//"_"//trim(method)//"_T"//trim(Tch)//"_SLangle_phi_theta90_nJex"//trim(nJexch)//"_rd"//trim(rdch)//".txt")
it=1
theta=0.5d0*pi
do ip=0,100 ! \phi=pi*dble(ip)/8.0d0
   phi=pi*dble(ip)/100.0d0
   part=0.0d0
   abc=0.0d0
   do ir=1,nr
      call Hexch(pi/2.0d0,phi,imax,zH(:,:,ir),zH0(:,:,ir),zspin,Hex(ir))
      call zheev('V','L',imax,zH(:,:,ir),imax,eigen,WORK,LWORK,RWORK,INFO)
      if(temp.eq.0.0d0) then
         freeloc(ir)=eigen(1) 
      else
         do  ii=1,imax
            part(ir)=part(ir)+dexp(-(eigen(ii)-eigen(1))/temp)
         end do
         freeloc(ir)=eigen(1)-temp*dlog(part(ir))
      end if
      freeloc(ir)=freeloc(ir)-freet0(ir) 
      abc=abc+freeloc(ir)
      call momentR_exp(temp,imax,zH(:,:,ir),zspin,zorbit,part(ir),eigen,stspin0(:,ir),stspinS0(:,ir),stspinL0(:,ir))
      call vector_to_angle(stspin0(:,ir),thetaj0(ir),phij0(ir),absj0(ir))
      call vector_to_angle(stspinS0(:,ir),thetaS0(ir),phiS0(ir),absS0(ir))
      call vector_to_angle(stspinL0(:,ir),thetaL0(ir),phiL0(ir),absL0(ir))
   end do ! for ir
   gibbs=abc*dmult+KuTM
   write(85,"(4f17.10)") phi,freeloc(1),freeloc(2),gibbs*dcoff
   write(88,"(9f17.10)") phi,phiS0(1),phiS0(2),phiL0(1),phiL0(2),thetaS0(1),thetaS0(2),thetaL0(1),thetaL0(2)
end do ! for \phi ip


! temperature dependence
open(518,file=""//trim(sys)//"_"//trim(method)//"_weight_nJex"//trim(nJexch)//"_rd"//trim(rdch)//".txt")
open(519,file=""//trim(sys)//"_"//trim(method)//"_englevel_J0_nJex"//trim(nJexch)//"_rd"//trim(rdch)//".txt")
open(520,file=""//trim(sys)//"_"//trim(method)//"_englevel_J1_1_nJex"//trim(nJexch)//"_rd"//trim(rdch)//".txt")
open(521,file=""//trim(sys)//"_"//trim(method)//"_englevel_J1_2_nJex"//trim(nJexch)//"_rd"//trim(rdch)//".txt")
open(522,file=""//trim(sys)//"_"//trim(method)//"_englevel_J2_1_nJex"//trim(nJexch)//"_rd"//trim(rdch)//".txt")
open(523,file=""//trim(sys)//"_"//trim(method)//"_englevel_J2_2_nJex"//trim(nJexch)//"_rd"//trim(rdch)//".txt")
open(524,file=""//trim(sys)//"_"//trim(method)//"_englevel_J3_1_nJex"//trim(nJexch)//"_rd"//trim(rdch)//".txt")
open(525,file=""//trim(sys)//"_"//trim(method)//"_englevel_J3_2_nJex"//trim(nJexch)//"_rd"//trim(rdch)//".txt")
do itemp=0,int(Tc)
   temp=dble(itemp)
   abc=MsT_Kuzmin(temp,MTM0,Tc,s,p)/MTM0
   part=0.0d0
   do ir=1,nr
      Hex(ir)=-Hex0(ir)*abc
      call Hexch(0.0d0,0.0d0,imax,zH(:,:,ir),zH0(:,:,ir),zspin,Hex(ir))
      call zheev('V','L',imax,zH(:,:,ir),imax,eigen,WORK,LWORK,RWORK,INFO)
      if(temp.ne.0.0d0) then
         do  ii=1,imax ! state sum
            part(ir)=part(ir)+dexp(-(eigen(ii)-eigen(1))/temp)
         end do
      end if
   end do ! do ir 
   write(519,"(7F15.5)") temp,eigen(1)-SO(1),eigen(2)-SO(1),eigen(3)-SO(1)&
        ,eigen(4)-SO(1),eigen(5)-SO(1),eigen(6)-SO(1)
   write(520,"(7F15.5)") temp,eigen(7)-SO(1),eigen(8)-SO(1),eigen(9)-SO(1)&
        ,eigen(10)-SO(1),eigen(11)-SO(1),eigen(12)-SO(1)
   write(521,"(3F15.5)") temp,eigen(13)-SO(1),eigen(14)-SO(1)
   write(522,"(7F15.5)") temp,eigen(15)-SO(1),eigen(16)-SO(1),eigen(17)-SO(1)&
        ,eigen(18)-SO(1),eigen(19)-SO(1),eigen(20)-SO(1)
   write(523,"(5F15.5)") temp,eigen(21)-SO(1),eigen(22)-SO(1),eigen(23)-SO(1),&
        eigen(24)-SO(1)
   write(524,"(7F15.5)") temp,eigen(25)-SO(1),eigen(26)-SO(1),eigen(27)-SO(1)&
        ,eigen(28)-SO(1),eigen(29)-SO(1),eigen(30)-SO(1)
   write(525,"(7F15.5)") temp,eigen(31)-SO(1),eigen(32)-SO(1),eigen(33)-SO(1),&
        eigen(34)-SO(1),eigen(35)-SO(1),eigen(36)-SO(1)
   ! Weight for J-multiplet
   ir=1
   J1=JM(1,1)
   WJ=0.0d0
   do jj=1,imax
      do i=0,nJex
         if(JM(jj,1).eq.J1+dble(i)) then
            if(temp.eq.0.0d0) then
               WJ(i)=WJ(i)+dreal(conjg(zH(jj,1,ir))*zH(jj,1,ir))
            else
               do m=1,imax
                  WJ(i)=WJ(i)+dreal(conjg(zH(jj,m,ir))*zH(jj,m,ir))*&
                       dexp(-(eigen(m)-eigen(1))/temp)/part(ir)
               end do
            end if
         end if
      end do ! do i
   end do ! do jj
   abc=0.0d0
   do i=0,nJex
      abc=abc+WJ(i)
   end do
   write(518,"(8F15.5)") temp,WJ(0),WJ(1),WJ(2),WJ(3),WJ(4),WJ(5),abc
end do ! do itemp   

end program energy_statistical


