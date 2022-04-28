subroutine CF_Ki_surface_Tdept_nobad(ax,cx,Hex,temp0,Klcr,Klcfe,nr,nfe,nunit,irsub,KuFe,HexT)
implicit none
real(8),parameter :: joule=10.0d0 ! (Merg/J)
integer(4) :: nr,nfe,nunit
real(8) :: ax,cx,eng0(4)
real(8) :: J(15),S,Lang,Jzst,Jz,Ost,Ost2,MOst(6,6),abc,stevens(6,15),Blm(6,0:6,16),Blm2(6,0:6),g(15),absj(16),absj0(16)
real(8) :: alpha(15),beta(15),gamma(15),Alm(6,0:6,16),Ecef,Exc,stspin(3),fct(16)
real(8) :: f(6),angle2,pi=dacos(-1.0d0),Eaniso,def,theta,theta0,thetaj(16),phij(16),phi,thetaj0(16),phij0(16),phi0
real(8) :: K11,K22,K1,K2,K3,K4,K5,K11loc(16),K22loc(16),K1loc(16),K2loc(16),K3loc(16),K4loc(16),K5loc(16),Kloc(9,16,0:700),Klocave(9,0:700),HexT(700)
real(8) :: Hex,part(16),temp,free,free0,K2mK3,Dtheta,freetp(0:10000,0:100),freeapp(16),freetploc(0:10000,0:100,16)
real(8),allocatable :: matrix(:,:),eigen(:),spin(:,:),RWORK(:)
complex(kind(0d0)),allocatable :: zmatrix(:,:),ZOst(:,:,:,:),WORK(:),mag(:,:),zmag(:,:),zmoment(:,:,:)
complex(kind(0d0)) :: zabc,zdef
integer(4) :: m,n,i,k,l,imax,ii,jj,INFO,LWORK,it,ip,itnum=400,ll,mm,ir,inddy,idim(15),itemp,ife
character(20) :: rare,Tch,Hexch
real(8) :: aleng,bleng,cleng,volume,angle,K1til(0:2),K2til(0:2),temp0,K1tilloc(0:2,16),K2tilloc(0:2,16),Klcr(9,nr),Klcfe(9,nfe*nunit),KuFe
integer(4) :: ix,iy,iz,ik,irsub(nr*nunit)
real(8),parameter :: kb=1.3806488d-23 ! (J/K)

temp=temp0
! Fe sublattice anisotropy
if(temp.eq.0.0d0.or.temp.eq.4.2d0) then
   KuFe=0.789d0*1.0d6*joule*(ax**2*cx*1.0d-30)/dble(nfe) !(Merg) Hirosawa
else if(temp.eq.200.0d0) then
   KuFe=1.05d0*1.0d6*joule*(ax**2*cx*1.0d-30)/dble(nfe) !(Merg) Hirosawa
else if(temp.eq.300.0d0) then
   KuFe=1.13d0*1.0d6*joule*(ax**2*cx*1.0d-30)/dble(nfe) ! 1.0 (MJ/m^3)
else if(temp.eq.325.0d0) then
   KuFe=1.125d0*1.0d6*joule*(ax**2*cx*1.0d-30)/dble(nfe) ! 1.0 (MJ/m^3)
else if(temp.eq.350.0d0) then
   KuFe=1.0875*1.0d6*joule*(ax**2*cx*1.0d-30)/dble(nfe) ! 1.0 (MJ/m^3)
else if(temp.eq.375.0d0) then
   KuFe=1.0375*1.0d6*joule*(ax**2*cx*1.0d-30)/dble(nfe) ! 1.0 (MJ/m^3)
else if(temp.eq.400.0d0) then
   KuFe=0.97d0*1.0d6*joule*(ax**2*cx*1.0d-30)/dble(nfe)
else if(temp.eq.425.0d0) then
   KuFe=0.9d0*1.0d6*joule*(ax**2*cx*1.0d-30)/dble(nfe)
else if(temp.eq.450.0d0) then
   KuFe=0.8d0*1.0d6*joule*(ax**2*cx*1.0d-30)/dble(nfe)
else if(temp.eq.475.0d0) then
   KuFe=0.65d0*1.0d6*joule*(ax**2*cx*1.0d-30)/dble(nfe)
else if(temp.eq.500.0d0) then
   KuFe=0.5d0*1.0d6*joule*(ax**2*cx*1.0d-30)/dble(nfe)
else 
   write(6,*) "Fe Ku error"
   stop
end if

do ife=1,nfe*nunit
   Klcfe(1,ife)=KuFe
   do ik=2,9
      Klcfe(ik,ife)=0.0d0
   end do
end do

rare="Nd"
write(Hexch,'(I5)') int(Hex)
write(Tch,'(I5)') int(temp0)

Dtheta=0.001d0

!do inddy=1,5
if(rare.eq."Nd") inddy=1
if(rare.eq."Dy") inddy=2
if(rare.eq."Pr") inddy=3
if(rare.eq."Tb") inddy=4
if(rare.eq."Ho") inddy=5

! stevens parameter
stevens=0.0d0
!if(rare.eq."Nd") then
! inddy=1
   open(34,file="~/CFPs/cfp"//trim(adjustl(rare))//"opencore_k100.txt")
   open(35,file="Ki_surfaceNd11.txt")
! f1,f2,g1,g2 according to Yamada et al.
! n=1
   open(66,file="freeeng_f2_phi0_"//trim(adjustl(rare))//"_Hex"//trim(adjustl(Hexch))//"_T"//trim(adjustl(Tch))//".txt")
! n=2
   open(77,file="freeeng_g2_phi0_"//trim(adjustl(rare))//"_Hex"//trim(adjustl(Hexch))//"_T"//trim(adjustl(Tch))//".txt")
! n=3
   open(666,file="freeeng_f1_phi0_"//trim(adjustl(rare))//"_Hex"//trim(adjustl(Hexch))//"_T"//trim(adjustl(Tch))//".txt")
! n=4
   open(777,file="freeeng_g1_phi0_"//trim(adjustl(rare))//"_Hex"//trim(adjustl(Hexch))//"_T"//trim(adjustl(Tch))//".txt")

   open(88,file="freeeng_f2_phi0.25pi_"//trim(adjustl(rare))//"_Hex"//trim(adjustl(Hexch))//"_T"//trim(adjustl(Tch))//".txt")
   open(99,file="freeeng_g2_phi0.25pi_"//trim(adjustl(rare))//"_Hex"//trim(adjustl(Hexch))//"_T"//trim(adjustl(Tch))//".txt")
   open(888,file="freeeng_f1_phi0.25pi_"//trim(adjustl(rare))//"_Hex"//trim(adjustl(Hexch))//"_T"//trim(adjustl(Tch))//".txt")
   open(999,file="freeeng_g1_phi0.25pi_"//trim(adjustl(rare))//"_Hex"//trim(adjustl(Hexch))//"_T"//trim(adjustl(Tch))//".txt")

   open(166,file="freeeng_f2_theta0.5pi_"//trim(adjustl(rare))//"_Hex"//trim(adjustl(Hexch))//"_T"//trim(adjustl(Tch))//".txt")
   open(177,file="freeeng_g2_theta0.5pi_"//trim(adjustl(rare))//"_Hex"//trim(adjustl(Hexch))//"_T"//trim(adjustl(Tch))//".txt")
   open(1666,file="freeeng_f1_theta0.5pi_"//trim(adjustl(rare))//"_Hex"//trim(adjustl(Hexch))//"_T"//trim(adjustl(Tch))//".txt")
   open(1777,file="freeeng_g1_theta0.5pi_"//trim(adjustl(rare))//"_Hex"//trim(adjustl(Hexch))//"_T"//trim(adjustl(Tch))//".txt")

   open(101,file="K1_surface_nobad_"//trim(adjustl(rare))//".txt")
   open(102,file="K2_surface_nobad_"//trim(adjustl(rare))//".txt")
   open(106,file="Kave_surface_nobad_"//trim(adjustl(rare))//".txt")
if(inddy.eq.1) then ! Nd
   J(inddy)=9.0d0/2.0d0
   S=3.0d0/2.0d0
   Lang=J(inddy)+S
   alpha(inddy)=-7.0d0/(3.0d0**2*11.0d0**2)
   beta(inddy)=-2.0d0**3*17.0d0/(3.0d0**3*11.0d0**3*13.0d0)
   gamma(inddy)=-5.0d0*17.0d0*19.0d0/(3.0d0**3*7.0d0*11.0d0**3*13.0d0**2)
   idim(inddy)=10
   aleng=8.81d0 ! angstrom
   bleng=aleng
   cleng=12.21d0 ! angstrom
   angle=pi/2.0d0
else if(inddy.eq.2) then ! Dy
   J(inddy)=15.0d0/2.0d0
   S=5.0d0/2.0d0
   Lang=J(inddy)-S
   alpha(inddy)=-2.0d0/(3.0d0**2*5.0d0*7.0d0)
   beta(inddy)=-2.0d0**3/(3.0d0**3*5.0d0*7.0d0*11.0d0*13.0d0)
   gamma(inddy)=2.0d0**2/(3.0d0**3*7.0d0*11.0d0**2*13.0**2)
   idim(inddy)=16
   aleng=8.76d0 ! angstrom
   bleng=aleng
   cleng=11.99d0 ! angstrom
   angle=pi/2.0d0
else if(inddy.eq.3) then ! Pr
   J(inddy)=4.0d0
   S=1.0d0
   Lang=J(inddy)+S
   alpha(inddy)=-2.0d0**2*13.0d0/(3.0d0**2*5.0d0**2*11.0d0)
   beta(inddy)=-2.0d0**2/(3.0d0**2*5.0d0*11.0d0**2)
   gamma(inddy)=2.0d0**4*17.0d0/(3.0d0**4*5.0d0*7.0d0*11.0d0**2*13.0d0)
   idim(inddy)=9
   aleng=8.81d0 ! angstrom
   bleng=aleng
   cleng=12.27d0 ! angstrom
   angle=pi/2.0d0
else if(inddy.eq.4) then ! Tb
   J(inddy)=6.0d0
   S=3.0d0
   Lang=J(inddy)-S
   alpha(inddy)=-1.0/(3.0d0**2*11.0d0)
   beta(inddy)=2.0d0/(3.0d0**3*5.0d0*11.0d0**2)
   gamma(inddy)=-1.0d0/(3.0d0**4*7.0d0*11.0d0**2*13.0d0)
   idim(inddy)=13
   aleng=8.81d0 ! angstrom
   bleng=aleng
   cleng=12.27d0 ! angstrom
   angle=pi/2.0d0
else if(inddy.eq.5) then ! Ho
   J(inddy)=8.0d0
   S=2.0d0
   Lang=J(inddy)-S
   alpha(inddy)=-1.0d0/(2.0d0*3.0d0**2*5.0d0**2)
   beta(inddy)=-1.0d0/(2.0d0*3.0d0*5.0d0*7.0d0*11.0d0*13.0d0)
   gamma(inddy)=-5.0d0/(3.0d0**3*7.0d0*11.0d0**2*13.0d0**2)
   idim(inddy)=17
   aleng=8.75d0 ! angstrom
   bleng=aleng
   cleng=11.99d0 ! angstrom
   angle=pi/2.0d0
end if
   g(inddy)=1.0d0+(J(inddy)*(J(inddy)+1.0d0)+S*(S+1.0d0)-Lang*(Lang+1.0d0))/(2.0d0*J(inddy)*(J(inddy)+1.0d0))

   imax=idim(inddy)
allocate(matrix(imax,imax),eigen(imax),spin(3,imax),RWORK(3*imax-2))
allocate(zmatrix(imax,imax),ZOst(imax,imax,6,0:6),WORK(2*imax-1),mag(imax,imax),zmag(imax,imax),zmoment(imax,imax,3))

stevens(2,inddy)=alpha(inddy)
stevens(4,inddy)=beta(inddy)
stevens(6,inddy)=gamma(inddy)

f(2)=J(inddy)*(J(inddy)-0.5d0)
f(4)=J(inddy)*(J(inddy)-0.5d0)*(J(inddy)-1.0d0)*(J(inddy)-1.5d0)
f(6)=J(inddy)*(J(inddy)-0.5d0)*(J(inddy)-1.0d0)*(J(inddy)-1.5d0)*(J(inddy)-2.0d0)*(J(inddy)-2.5d0)

Blm=0.0d0
do n=1,2
   read(34,*) rare
   do l=2,6,2
      do m=0,l,2
         read(34,*) Alm(l,m,n)
         if(n.eq.1) then
            Blm(l,m,1)=Alm(l,m,n)*stevens(l,inddy)
            Blm(l,m,2)=Alm(l,m,n)*stevens(l,inddy)
            Blm(l,m,5)=Alm(l,m,n)*stevens(l,inddy)
            Blm(l,m,6)=Alm(l,m,n)*stevens(l,inddy)
         else 
            Blm(l,m,3)=Alm(l,m,n)*stevens(l,inddy)
            Blm(l,m,4)=Alm(l,m,n)*stevens(l,inddy)
            Blm(l,m,7)=Alm(l,m,n)*stevens(l,inddy)
            Blm(l,m,8)=Alm(l,m,n)*stevens(l,inddy)
         end if
      end do
   end do
end do
do l=2,6,2
   do m=2,6,4
      Blm(l,m,1)=-Blm(l,m,5)
      Blm(l,m,2)=-Blm(l,m,5)
      Blm(l,m,3)=-Blm(l,m,7)
      Blm(l,m,4)=-Blm(l,m,7)
   end do
end do

call Stevenseq(imax,J(inddy),ZOst)

fct=1.0d0
do ir=1,4
   fct(ir)=0.5d0
end do


do itemp=1,700
temp=dble(itemp)

! ground state energy
Ecef=0.0d0
do i=1,8
   zmatrix=0.0d0
   do jj=1,imax
      do ii=jj,imax
         do l=2,6,2
            do m=0,l,2
               zmatrix(ii,jj)=zmatrix(ii,jj)+ZOst(ii,jj,l,m)*Blm(l,m,i)
            end do
         end do
      end do
   end do

   LWORK=2*imax-1
   call zheev('V','L',imax,zmatrix,imax,eigen,WORK,LWORK,RWORK,INFO)
   Ecef=Ecef+eigen(imax)-eigen(1)

!   do ii=1,imax
!      write(6,*) eigen (ii)
!   end do
end do ! for rare earth site

call moment(imax,J(inddy),zmoment,imax)

!do jj=1,imax
!   do ii=jj,imax
!      write(6,*) ii,jj
!      zabc=zmag(ii,jj)-(zmoment(ii,jj,1)*dsin(theta)/dsqrt(2.0d0)+zmoment(ii,jj,2)*dsin(theta)/dsqrt(2.0d0)+zmoment(ii,jj,3)*dcos(theta))
!      if(abs(zabc).gt.1.0d-10) write(6,*) abs(zabc)
!   end do
!end do
!stop ! change do it

free0=0.0d0
theta0=0.0d0
phi0=0.0d0
thetaj0=0.0d0
phij0=0.0d0
absj0=0.0d0
freetp=0.0d0
do ip=0,8 ! \phi=0
   phi=pi*dble(ip)/8.0d0
!   write(6,*) "phi=",phi/pi,"pi"
   do it=0,int(pi/Dtheta/2.0d0) ! take three point
      if((it.ge.3.and.temp0.ne.temp).or.(ip.ge.3.and.temp0.ne.temp)) cycle
!      write(6,*) it,ip
      theta=Dtheta*dble(it)
      !   theta=pi/6.0d0
      !   write(6,*) "theta=",theta/pi*180.0d0
      call exch(imax,J(inddy),zmag,theta,phi,imax)

      free=0.0d0
      part=0.0d0
      do ir=1,8 ! sum up inequivalent R site
         zmatrix=0.0d0
         do jj=1,imax
            do ii=jj,imax
               do l=2,6,2
                  do m=0,l,2
                     zmatrix(ii,jj)=zmatrix(ii,jj)+ZOst(ii,jj,l,m)*Blm(l,m,ir)
                  end do
               end do
            end do
         end do
         do jj=1,imax
            do ii=1,imax
               zmatrix(ii,jj)=zmatrix(ii,jj)+2.0d0*(g(inddy)-1.0d0)*fct(ir)*HexT(itemp)*zmag(ii,jj)
            end do
         end do
         
         LWORK=2*imax-1
         call zheev('V','L',imax,zmatrix,imax,eigen,WORK,LWORK,RWORK,INFO)
         def=0.0d0
         
         if(temp.eq.0.0d0) then
            free=free+eigen(1) 
!            write(6,*) eigen(1)
            freetploc(it,ip,ir)=eigen(1)
         else
            do  ii=1,imax
               part(ir)=part(ir)+dexp(-eigen(ii)/temp)
               !         write(6,*) "part",part(ir)
            end do
            free=free-temp*dlog(part(ir))
            freetploc(it,ip,ir)=-temp*dlog(part(ir))
         end if
         
         ! mag moment
         spin=0.0d0
         stspin=0.0d0
         if(temp.eq.0.0d0) then
            do ll=1,3
               do jj=1,imax
                  do ii=1,imax
                     spin(ll,1)=spin(ll,1)+dreal(conjg(zmatrix(ii,1))*zmoment(ii,jj,ll)*zmatrix(jj,1))
                  end do
               end do
               stspin(ll)=spin(ll,1)
            end do ! for spin component
         else
            do ll=1,3
               do m=1,imax
                  
                  do jj=1,imax
                     do ii=1,imax
                        spin(ll,m)=spin(ll,m)+dreal(conjg(zmatrix(ii,m))*zmoment(ii,jj,ll)*zmatrix(jj,m))
                     end do
                  end do
                  stspin(ll)=stspin(ll)+spin(ll,m)*dexp(-eigen(m)/temp)/part(ir)
               end do ! for eng level
            end do ! for spin component
         end if
         
         absj(ir)=0.0d0
         do ll=1,3
            absj(ir)=absj(ir)+stspin(ll)**2
            !      write(6,*) stspin(ll)
         end do
         absj(ir)=dsqrt(absj(ir))
         thetaj(ir)=dacos(stspin(3)/absj(ir))
         phij(ir)=dacos(stspin(1)/dsqrt(stspin(1)**2+stspin(2)**2))
         !   write(6,*) stspin(1),stspin(2),stspin(3)
         if(stspin(2).ge.0.0d0) then
            phij(ir)=phij(ir)
         else
            phij(ir)=2.0d0*pi-phij(ir)
         end if
      end do ! for rare earth site ir
      
      freetp(it,ip)=free
      
      if(free.lt.free0) then
         !      write(6,*) theta0/pi*180.0d0,free,free0
         free0=free
         theta0=theta
         phi0=phi
         thetaj0=thetaj
         phij0=phij
         absj0=absj
      end if
   end do ! for \theta it
end do ! for \phi ip

volume=2.0d0*kb/(aleng*bleng*cleng*dsin(angle)*1.0d-30)*1.0d-6 ! temperature to J/m^3
do ir=1,8
   do i=0,2 ! derivative 
      K1tilloc(i,ir)=0.5d0*(2.0d0*(freetploc(1,i,ir)-freetploc(0,i,ir))/Dtheta**2)
      K2tilloc(i,ir)=K1tilloc(i,ir)/3.0d0+&
           (2.0d0*(freetploc(2,i,ir)-4.0d0*freetploc(1,i,ir)+3.0d0*freetploc(0,i,ir))/Dtheta**4)/24.0d0
   end do
end do

Klocave=0.0d0
do ir=1,8
   Kloc(1,ir,itemp)=K1tilloc(0,ir)
   Kloc(6,ir,itemp)=K1tilloc(2,ir)-Kloc(1,ir,itemp)
   Kloc(7,ir,itemp)=(dsqrt(2.0d0)+1.0d0)*(2.0d0*K2tilloc(1,ir)-K2tilloc(0,ir)-K2tilloc(2,ir))
   Kloc(2,ir,itemp)=K2tilloc(1,ir)-Kloc(7,ir,itemp)/dsqrt(2.0d0)
   Kloc(3,ir,itemp)=K2tilloc(0,ir)-Kloc(2,ir,itemp)

   Klocave(1,itemp)=Klocave(1,itemp)+Kloc(1,ir,itemp)
   Klocave(2,itemp)=Klocave(2,itemp)+Kloc(2,ir,itemp)
   Klocave(3,itemp)=Klocave(3,itemp)+Kloc(3,ir,itemp)
   Klocave(6,itemp)=Klocave(6,itemp)+Kloc(6,ir,itemp)
   Klocave(7,itemp)=Klocave(7,itemp)+Kloc(7,ir,itemp)
end do
Klocave=Klocave/8.0d0
write(106,*) temp,Klocave(1,itemp)

!write(6,*) irsub(1),irsub(3),irsub(5),irsub(7)
!write(6,*) 3,4,1,2
!stop

if(temp.eq.temp0) then
!   Ku=(K1loc(1)+K1loc(2)+K1loc(3)+K1loc(4))/4.0d0
   Klcr=0.0d0
   do ir=1,8
      do ik=1,9
         Klcr(ik,ir)=Kloc(ik,ir,itemp)
      end do
   end do


   phi=0.0d0
   do it=1,int(pi/Dtheta/2.0d0)
      theta=dble(it)*Dtheta
      do ir=1,4
         freeapp(ir)=(Kloc(1,ir,itemp)+Kloc(6,ir,itemp)*dsin(2.0d0*phi))*dsin(theta)**2+(Kloc(2,ir,itemp)+Kloc(7,ir,itemp)*dsin(2.0d0*phi)+Kloc(3,ir,itemp)*dcos(4.0d0*phi))*dsin(theta)**4
      end do

      write(66,*) 2.0d0*theta/pi,(freetploc(it,0,1)-freetploc(0,0,1)),freeapp(1) ! Kelvin
      write(77,*) 2.0d0*theta/pi,(freetploc(it,0,2)-freetploc(0,0,2)),freeapp(2)
      write(666,*) 2.0d0*theta/pi,(freetploc(it,0,3)-freetploc(0,0,3)),freeapp(3)
      write(777,*) 2.0d0*theta/pi,(freetploc(it,0,4)-freetploc(0,0,4)),freeapp(4)
!      if(it.eq.int(pi/Dtheta/2.0d0)) write(6,*) "freeapp(3)",theta,phi,freeapp(3)
   end do

   phi=pi/4.0d0
   do it=1,int(pi/Dtheta/2.0d0)
      theta=dble(it)*Dtheta
      do ir=1,4
         freeapp(ir)=(Kloc(1,ir,itemp)+Kloc(6,ir,itemp)*dsin(2.0d0*phi))*dsin(theta)**2+(Kloc(2,ir,itemp)+Kloc(7,ir,itemp)*dsin(2.0d0*phi)+Kloc(3,ir,itemp)*dcos(4.0d0*phi))*dsin(theta)**4
      end do
      write(88,*) 2.0d0*theta/pi,(freetploc(it,2,1)-freetploc(0,2,1)),freeapp(1) ! Kelvin
      write(99,*) 2.0d0*theta/pi,(freetploc(it,2,2)-freetploc(0,2,2)),freeapp(2)
      write(888,*) 2.0d0*theta/pi,(freetploc(it,2,3)-freetploc(0,2,3)),freeapp(3)
      write(999,*) 2.0d0*theta/pi,(freetploc(it,2,4)-freetploc(0,2,4)),freeapp(4)
   end do
   do i=1,4
      eng0(i)=freetploc(0,0,i)
   end do

end if ! if temp.eq.temp0
!if(temp.eq.temp0) write(6,*) "freeloc",
end do



do ir=1,8
   do itemp=0,700
      temp=dble(itemp)
      write(101,*) temp,Kloc(1,ir,itemp),ir
      write(102,*) temp,Kloc(2,ir,itemp),ir
      if(ir.eq.10.and.temp.eq.300.0d0) then
         write(35,*) "T, K1, K2, K3, K6, K7"
         write(35,*) temp,Kloc(1,ir,itemp),Kloc(2,ir,itemp),Kloc(3,ir,itemp),Kloc(6,ir,itemp),Kloc(7,ir,itemp)
      end if
!   write(103,*) temp,Kloc(3,1),Kloc(3,2),(Kloc(3,1)+Kloc(3,2)+Kloc(3,3)+Kloc(3,4))*volume
!   write(104,*) temp,Kloc(6,1),Kloc(6,2),(Kloc(6,1)+Kloc(6,2)+Kloc(6,3)+Kloc(6,4))*volume
!   write(105,*) temp,Kloc(7,1),Kloc(7,2),(Kloc(7,1)+Kloc(7,2)+Kloc(7,3)+Kloc(7,4))*volume
   end do
   write(101,*) 
   write(102,*) 
end do

end subroutine CF_Ki_surface_Tdept_nobad

