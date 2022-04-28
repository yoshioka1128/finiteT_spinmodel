program magcurve_statistical
implicit none
integer(4),parameter :: nrmax=30
real(8),parameter :: pi=dacos(-1.0d0)
real(8),parameter :: mu0=4.0d0*pi*1.0d-7
real(8),parameter :: kb=1.3806488d-23 ! (J/K)
real(8),parameter :: muB=9.2740100783d-24 ! (Merg/kOe) or (kerg/Oe) or (J/T)
real(8),parameter :: joule=10.0d0 ! (Merg/J)
integer(4) :: nJex,icount,iHmag,nfe,nr,n4f,m,n,i,k,l,imax,ii,jj,it,ip,ll,ir,ncfp,ipmax,itmax,INFO,LWORK
real(8) :: S1,L1,J1,lambda,Hmax,dH,rd,dmult,abc,def
real(8) :: KuFe,Femag0,Rmag0,Tc,s,p,mFe0,Ms,MsT_Kuzmin,Hmag,thetamag,phimag,thetatot,phitot,abstot
real(8) :: aleng,bleng,cleng,gmm,volume,K1Fe,rK2Fe,rK3Fe,temp,free,free0,dtheta,dphi,dphi0
real(8) :: stfct(6,-6:6),JM(100,2),SO(100),cclm(6,-6:6)
real(8) :: Alm(6,-6:6,nrmax),Clm(6,-6:6,nrmax),rl(6),g(100),rdcoff(0:6)
real(8) :: absj0(nrmax),absS0(nrmax),absL0(nrmax),thetaj0(nrmax),phij0(nrmax)
real(8) :: thetaS0(nrmax),phiS0(nrmax),thetaL0(nrmax),phiL0(nrmax)
real(8) :: stspin0(3,nrmax),stspinS0(3,nrmax),stspinL0(3,nrmax)
real(8) :: theta,theta0,phi,phi0,freeloc(nrmax),freeloc0(0:2,0:2,nrmax)
real(8) :: part(nrmax),part0(nrmax),Hex(nrmax),Hex0(nrmax),magtot(3)
real(8),allocatable :: matrix(:,:),eigen(:),eigen0(:),RWORK(:),Ust(:,:,:,:),Vst(:,:,:)
complex(kind(0d0)),allocatable :: ZUst(:,:,:,:),WORK(:),mag(:,:),zU(:,:,:)
complex(kind(0d0)),allocatable :: zspin(:,:,:),zorbit(:,:,:),zH(:,:,:),zH0(:,:,:),zH00(:,:,:)
character(50) :: Tch,sys,direction,nJexch,ratom,tatom,method,rdch,wwoK1Fe

! set parameter
read(5,*) sys
read(5,*) temp ! temperature (K)
read(5,*) nJex ! 0:lowestJ, 1:LowestJ+secondJ
read(5,*) method
read(5,*) wwoK1Fe
read(5,*) direction
read(5,*) Hmax,dH
read(5,*) dtheta
read(5,*) dphi0
write(6,*) "system=",sys
write(6,*) "T=",temp,"(K)"
write(6,*) "number of excited states",nJex
write(6,*) "method=",method
write(6,*) "wwoK1Fe",wwoK1Fe
write(6,*) "direction=",direction
write(6,*) "applied field",Hmax,dH,"(kOe)"
write(6,*) "increment of \theta",dtheta,"(rad)",dtheta*180.0d0/pi,"deg."
write(6,*) "increment of \phi",dphi0,"(rad)",dphi0*180.0d0/pi,"deg."
write(Tch,'(f5.1)') temp
write(nJexch,"(I2)") nJex
Tch=adjustl(Tch)
nJexch=adjustl(nJexch)
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

! input files
open(33,file="~/research/CFPs/Sm_lambda411K/cfp"//trim(sys)//"_"//trim(method)//".txt")
read(33,*) ratom,nr,tatom,nfe
read(33,*) dmult
read(33,*) aleng,bleng,cleng,gmm
read(33,*) lambda
read(33,*) mFe0,Tc,s,p
read(33,*) K1Fe,rK2Fe,rK3Fe
write(6,*) "system=",ratom,nr,tatom,nfe
write(6,*) "(a,b,c)=",aleng,bleng,cleng,"(angstrom)"
write(6,*) "\lambda=",lambda,"(K)"
write(6,*) "Ms, Tc, s, p for Kuzmin fitting"
write(6,*) mFe0, Tc, s, p
write(6,*) "K1Fe,rK2Fe,rK3Fe",K1Fe,rK2Fe,rK3Fe
Alm=0.0d0
do ir=1,nr
   read(33,*)
   read(33,*) Hex0(ir),rd ! Hex
   Hex0(ir)=rd*Hex0(ir)
   read(33,*) rl(2),rl(4),rl(6),ncfp
   do i=1,ncfp
      read(33,*) abc,l,m
      Alm(l,m,ir)=abc*rl(l)
   end do
end do
write(rdch,"(f5.1)") rd
rdch=adjustl(rdch)
if(nr.gt.nrmax) then
   write(6,*) "error: nr > nrmax"
   stop
end if
volume=aleng*bleng*cleng*dsin(gmm/180.0d0*pi)*1.0d-30
write(6,*) "Fe sublattice magnetization @T=0 [T]"
write(6,*) mu0*mFe0*muB/volume

! Fe sublattice
Ms=MsT_Kuzmin(temp,mFe0,Tc,s,p) ! mag. moment per nfe
abc=Ms/mFe0
do ir=1,nr
   Hex(ir)=-Hex0(ir)*abc
end do
def=abc**3+&
     8.0d0/7.0d0*rK2Fe*(abc**3-abc**10)+&
     8.0d0/7.0d0*rK3Fe*(abc**3-18.0d0/11.0d0*abc**10+7.0d0/11.0d0*abc**21)
if(wwoK1Fe.eq."wK1Fe") then
   KuFe=K1Fe*def
else if(wwoK1Fe.eq."woK1Fe") then
   KuFe=0.0d0
else 
   write(6,*) "wwoK1Fe read error"
   stop
end if
write(6,*) "MsFe=",Ms,"[muB/unitcell]"
write(6,*) "KuFe=",KuFe,"[K/unitcell]"

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
call magcurve_analytical(temp,sys,method,wwoK1Fe,volume,lambda,Alm,rdch,KuFe,Hex,Ms,ratom,nr,dmult,direction)

! statistical calculations
allocate(matrix(imax,imax),eigen(imax),eigen0(imax),RWORK(3*imax-2),&
     Ust(imax,imax,6,-6:6),Vst(imax,imax,-1:1))
allocate(ZUst(imax,imax,6,-6:6),WORK(2*imax-1),mag(imax,imax),&
     zU(imax,imax,nr),zH(imax,imax,nr),zH0(imax,imax,nr),zH00(imax,imax,nr),zspin(imax,imax,3),zorbit(imax,imax,3))
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
call moment_op(imax,L1,S1,zspin,zorbit,JM,Vst,Ust)

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
zH00=(0.0d0,0.0d0)
do ir=1,nr
   do jj=1,imax
      zH00(jj,jj,ir)=SO(jj)
      do ii=jj,imax
         do l=2,6,2
            do m=-l,l
               zH00(ii,jj,ir)=zH00(ii,jj,ir)+ZUst(ii,jj,l,m)*Clm(l,m,ir)
            end do
         end do
      end do
   end do
end do

open(70,file=""//trim(sys)//"_"//trim(method)//"_"//trim(wwoK1Fe)//"_&
     magcurve_T"//trim(Tch)//"_"//trim(direction)//"_nJex"//trim(nJexch)//"_rd"//trim(rdch)//".txt")
open(71,file=""//trim(sys)//"_"//trim(method)//"_"//trim(wwoK1Fe)//"_&
     Fe_theta-phi_T"//trim(Tch)//"_"//trim(direction)//"_nJex"//trim(nJexch)//"_rd"//trim(rdch)//".txt")
open(72,file=""//trim(sys)//"_"//trim(method)//"_"//trim(wwoK1Fe)//"_&
     RE1-2_J_theta-phi_T"//trim(Tch)//"_"//trim(direction)//"_nJex"//trim(nJexch)//"_rd"//trim(rdch)//".txt")
open(73,file=""//trim(sys)//"_"//trim(method)//"_"//trim(wwoK1Fe)//"_&
     RE1-2_S_theta-phi_T"//trim(Tch)//"_"//trim(direction)//"_nJex"//trim(nJexch)//"_rd"//trim(rdch)//".txt")
open(74,file=""//trim(sys)//"_"//trim(method)//"_"//trim(wwoK1Fe)//"_&
     RE1-2_L_theta-phi_T"//trim(Tch)//"_"//trim(direction)//"_nJex"//trim(nJexch)//"_rd"//trim(rdch)//".txt")
open(75,file=""//trim(sys)//"_"//trim(method)//"_"//trim(wwoK1Fe)//"_&
     Mtot_theta-phi_T"//trim(Tch)//"_"//trim(direction)//"_nJex"//trim(nJexch)//"_rd"//trim(rdch)//".txt")
write(70,*) "# applied field (T), magnetization (2 f.u.), magnetization (T)"
write(71,*) "# applied field (T), theta (degree), phi (degree), abs"
write(72,*) "# applied field (T), theta (degree), phi (degree), abs for RE 1-2"
write(73,*) "# applied field (T), theta (degree), phi (degree), abs for RE 1-2"
write(74,*) "# applied field (T), theta (degree), phi (degree), abs for RE 1-2"
write(75,*) "# applied field (T), theta (degree), phi (degree), abs [muB/unitcell], abs [T]"
if(nr.gt.2) then
   open(82,file=""//trim(sys)//"_"//trim(method)//"_"//trim(wwoK1Fe)//"_&
        RE3-4_J_theta-phi_T"//trim(Tch)//"_"//trim(direction)//"_nJex"//trim(nJexch)//"_rd"//trim(rdch)//".txt")
   open(83,file=""//trim(sys)//"_"//trim(method)//"_"//trim(wwoK1Fe)//"_&
        RE3-4_S_theta-phi_T"//trim(Tch)//"_"//trim(direction)//"_nJex"//trim(nJexch)//"_rd"//trim(rdch)//".txt")
   open(84,file=""//trim(sys)//"_"//trim(method)//"_"//trim(wwoK1Fe)//"_&
        RE3-4_L_theta-phi_T"//trim(Tch)//"_"//trim(direction)//"_nJex"//trim(nJexch)//"_rd"//trim(rdch)//".txt")
   write(82,*) "# applied field (T), theta (degree), phi (degree), abs for RE 3-4"
   write(83,*) "# applied field (T), theta (degree), phi (degree), abs for RE 3-4"
   write(84,*) "# applied field (T), theta (degree), phi (degree), abs for RE 3-4"
end if
! applied field dependence
do iHmag=0,int(Hmax/dH) ! magnetic field (kOe)
   Hmag=dble(iHmag)*dH*muB/joule/kb ! (K)
   zH0=(0.0d0,0.0d0)
   do ir=1,nr
      do jj=1,imax
         do ii=1,imax
            zH0(ii,jj,ir)=zH00(ii,jj,ir)+&
                 ((zorbit(ii,jj,3)+2.0d0*zspin(ii,jj,3))*dcos(thetamag)&
                 +((zorbit(ii,jj,1)+2.0d0*zspin(ii,jj,1))*dcos(phimag)&
                 +(zorbit(ii,jj,2)+2.0d0*zspin(ii,jj,2))*dsin(phimag))*dsin(thetamag))*Hmag          
         end do
      end do
   end do
   free0=0.0d0
   ! freeneg
   itmax=int(0.5d0*pi/dtheta)
   do it=0,itmax ! take three point
      theta=0.5d0*pi*dble(it)/dble(itmax)
      if(it.ne.0) then
         dphi=dphi0/dsin(theta)
         ipmax=int(0.5d0*pi/dphi)+1
      else
         ipmax=1
      end if
      do ip=0,ipmax ! \phi=0
         phi=0.5d0*pi*dble(ip)/dble(ipmax)
         
         freeloc=0.0d0
         part=0.0d0
         do ir=1,nr ! sum up inequivalent R site
            call Hexch(theta,phi,imax,zH(:,:,ir),zH0(:,:,ir),zspin,Hex(ir))
            call zheev('V','L',imax,zH(:,:,ir),imax,eigen,WORK,LWORK,RWORK,INFO)
            if(temp.eq.0.0d0) then
               freeloc(ir)=eigen(1) 
            else
               do  ii=1,imax
                  part(ir)=part(ir)+dexp(-(eigen(ii)-eigen(1))/temp)
               end do
               freeloc(ir)=eigen(1)-temp*dlog(part(ir))
            end if
         end do ! do ir
         free=0.0d0
         do ir=1,nr
            free=free+freeloc(ir)
         end do
         ! free eng per unitcell including nr*dmult of R
         free=dmult*free-Hmag*Ms*(dsin(thetamag)*dsin(theta)*dcos(phimag-phi)+dcos(thetamag)*dcos(theta))&
              +KuFe*dsin(theta)**2
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
  
! magnetic moment on R ion 
   do ir=1,nr
      call momentR_exp(temp,imax,zU(:,:,ir),zspin,zorbit,part0(ir),eigen0,stspin0(:,ir),stspinS0(:,ir),stspinL0(:,ir))
      call vector_to_angle(stspin0(:,ir),thetaj0(ir),phij0(ir),absj0(ir))
      call vector_to_angle(stspinS0(:,ir),thetaS0(ir),phiS0(ir),absS0(ir))
      call vector_to_angle(stspinL0(:,ir),thetaL0(ir),phiL0(ir),absL0(ir))
   end do

! magnetization of Fe and RE along the applied field direction
   Femag0=Ms*(dsin(thetamag)*dsin(theta0)*dcos(phimag-phi0)+dcos(thetamag)*dcos(theta0))
   Rmag0=0.0d0
   do ir=1,nr
      Rmag0=Rmag0+(dsin(thetamag)*(stspin0(1,ir)*dcos(phimag)+stspin0(2,ir)*dsin(phimag))+stspin0(3,ir)*dcos(thetamag))
   end do
   Rmag0=Rmag0*dmult

! total magnetization
   magtot=0.0d0
   magtot(1)=Ms*dsin(theta0)*dcos(phi0)
   magtot(2)=Ms*dsin(theta0)*dsin(phi0)
   magtot(3)=Ms*dcos(theta0)
   do ir=1,nr
      do ll=1,3
         magtot(ll)=magtot(ll)+dmult*stspin0(ll,ir)
      end do
   end do
   call vector_to_angle(magtot,thetatot,phitot,abstot)
   
   ! output
   write(6,*) 
   write(6,*) "applied field",dble(iHmag)*dH/10.0d0,"(T)",Hmag,"(K)"
   do ir=1,nr
      write(6,"(a20,3f12.5)") "RE theta-phi-abs",thetaj0(ir)/pi*180.0d0,phij0(ir)/pi*180.0d0,absj0(ir)
   end do
   write(6,"(a20,3f12.5)") "Fe theta-phi-abs",theta0/pi*180.0d0,phi0/pi*180.0d0,Ms
   write(6,"(a20,3f12.5)") "total theta-phi-abs",thetatot/pi*180.0d0,phitot/pi*180.0d0,abstot
   write(6,*) "magnetization"," ",direction,Femag0+Rmag0,(Femag0+Rmag0)*muB*mu0/volume
   write(6,*) "energy",free0
   
   write(70,*) dble(iHmag)*dH/10.0d0,Femag0+Rmag0,(Femag0+Rmag0)*muB/(aleng**2*cleng*1.0d-30)*mu0
   write(71,"(4f12.5)") dble(iHmag)*dH/10.0d0,&
        theta0/pi*180.0d0,phi0/pi*180.0d0,Ms
   write(72,"(7f12.5)") dble(iHmag)*dH/10.0d0,&
        thetaj0(1)/pi*180.0d0,phij0(1)/pi*180.0d0,absj0(1),thetaj0(2)/pi*180.0d0,phij0(2)/pi*180.0d0,absj0(2)
   write(73,"(7f12.5)") dble(iHmag)*dH/10.0d0,&
        thetaS0(1)/pi*180.0d0,phiS0(1)/pi*180.0d0,absS0(1),thetaS0(2)/pi*180.0d0,phiS0(2)/pi*180.0d0,absS0(2)
   write(74,"(7f12.5)") dble(iHmag)*dH/10.0d0,&
        thetaL0(1)/pi*180.0d0,phiL0(1)/pi*180.0d0,absL0(1),thetaL0(2)/pi*180.0d0,phiL0(2)/pi*180.0d0,absL0(2)
   if(nr.gt.2) then
      write(82,"(7f12.5)") dble(iHmag)*dH/10.0d0,&
           thetaj0(3)/pi*180.0d0,phij0(3)/pi*180.0d0,absj0(3),thetaj0(4)/pi*180.0d0,phij0(4)/pi*180.0d0,absj0(4)
      write(83,"(7f12.5)") dble(iHmag)*dH/10.0d0,&
           thetaS0(3)/pi*180.0d0,phiS0(3)/pi*180.0d0,absS0(3),thetaS0(4)/pi*180.0d0,phiS0(4)/pi*180.0d0,absS0(4)
      write(84,"(7f12.5)") dble(iHmag)*dH/10.0d0,&
           thetaL0(3)/pi*180.0d0,phiL0(3)/pi*180.0d0,absL0(3),thetaL0(4)/pi*180.0d0,phiL0(4)/pi*180.0d0,absL0(4)
   end if
   write(75,"(8f17.8)") dble(iHmag)*dH/10.0d0,thetatot/pi*180.0d0,phitot/pi*180.0d0,abstot,muB*mu0*abstot/volume
   
end do ! iHmag

end program magcurve_statistical
