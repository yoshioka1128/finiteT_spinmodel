subroutine fundip_spin2(dip,TQ,f,ik,Happl,dK1,dK2,dK3,nunit,nunitx,nunity,nunitz,alphar,gammar,gammars,numint0,&
     nr,nrtot,JRR, magR, iIntrr, numintrr, vr,vrs,irunit,irsub,diprr,ratio,magRT,magRST,&
     nfe,nfetot,JRFe,JRFeS,magFe,iIntrfe,numintrfe,vfe,ifeunit,ifesub,diprfe)
implicit none 
integer(4) :: nr,nfe,numint0,nunit,nunitx,nunity,nunitz,idim,idim2,ik,iIntrr(nr*nunit,numint0),iIntrfe(nr*nunit,numint0)
integer(4) :: ir,ir0,ir00,ix,iy,iz,i,j,ixd,iyd,izd,i2,nrtot,nfetot
integer(4) :: ix2,iy2,iz2,ir2,ife2,ife00
integer(4) :: numintrr(nr*nunit),numintrfe(nr*nunit)
real(8) :: vr(nr*nunit,3),vrs(nr*nunit,3),vfe(nfe*nunit,3),f(nr*nunit,3,4),HE(3),HA(3),Hk(3),Happl(3),T(3)
real(8) :: vH,TQ,abc,def,ghi,gammar,gammars,JRFe,JRFeS,JRR,magR(nr*nunit),magFe(nfe*nunit),magRT,magRST
real(8) :: alphar(nr*nunit),rfept(nfe*nunit,3),rrpt(nr*nunit,3),dis,ratio
real(8) :: fein(nfe*nunit),rin(nr*nunit),dfept(nfe*nunit,3),drpt(nr*nunit,3),delta(3,3)
real(8) :: drfept(-nunitx:nunitx,-nunity:nunity,-nunitz:nunitz,nr,nfe,3),drrpt(-nunitx:nunitx,-nunity:nunity,-nunitz:nunitz,nr,nr,3)
real(8) :: diprr(-nunitx:nunitx,-nunity:nunity,-nunitz:nunitz,nr,nr,3,3),diprfe(-nunitx:nunitx,-nunity:nunity,-nunitz:nunitz,nr,nfe,3,3)
integer(4) :: irunit(nr*nunit,3),ifeunit(nfe*nunit,3),irsub(nr*nunit),ifesub(nfe*nunit)
character(10) :: dip
real(8) :: dK1(7,nr*nunit),dK2(7,nr*nunit),dK3(7,nr*nunit)

delta=0.0d0
do idim=1,3
   delta(idim,idim)=1.0d0
end do

TQ=0.0d0
do ir=1,nrtot
   ix=irunit(ir,1)
   iy=irunit(ir,2)
   iz=irunit(ir,3)
   i=irsub(ir)
   
!   if(magR(ir).eq.0.0d0) cycle
! exchange interaction
   HA=0.0d0
   do idim=1,3
      abc=0.0d0
      do j=1,numintrr(ir)
         abc=abc+vr(iIntrr(ir,j),idim)*abs(magR(iIntrr(ir,j)))
      end do
      def=0.0d0
      do j=1,numintrfe(ir)
         def=def+vfe(iIntrfe(ir,j),idim)*abs(magFe(iIntrfe(ir,j)))
      end do
!      HA(idim)=2.0d0*JRR*abc+JRFe*def ! factor 2*JRR is originate from double counting of interacting pair
      HA(idim)=2.0d0*JRR*abc+JRFeS*def ! factor 2*JRR is originate from double counting of interacting pair
!      write(6,*)  idim,HA(idim)
   end do

   if(dip.eq."on") then
      do idim=1,3
         abc=0.0d0
         do idim2=1,3

            do ir2=1,nrtot
               ix2=irunit(ir2,1)
               iy2=irunit(ir2,2)
               iz2=irunit(ir2,3)
               i2=irsub(ir2)
               abc=abc+diprr(ix-ix2,iy-iy2,iz-iz2,i,i2,idim,idim2)*vr(ir2,idim2)*magR(ir2)
            end do
            do ife2=1,nfetot
               ix2=ifeunit(ife2,1)
               iy2=ifeunit(ife2,2)
               iz2=ifeunit(ife2,3)
               i2=ifesub(ife2)
               abc=abc+diprfe(ix-ix2,iy-iy2,iz-iz2,i,i2,idim,idim2)*vfe(ife2,idim2)*magFe(ife2)
!               write(6,*) i2,ix2,iy2,iz2,abc
            end do

         end do ! idim2
         HA(idim)=HA(idim)-abc*1.0d24 ! angstrome^3 to cm^3
      end do ! idim
   end if ! dipole interaction      

! local anisotropy
   izd=3
   ixd=1
   iyd=2
   Hk=0.0d0
   do ixd=1,2
      if(ixd.eq.1) iyd=2
      if(ixd.eq.2) iyd=1
      Hk(ixd)=-(&
           2.0d0*dK1(1,ir)*vrs(ir,ixd)&
           +2.0d0*dK1(2,ir)*vrs(ir,iyd)&
           +2.0d0*dK1(3,ir)*vrs(ir,ixd)& ! new
           +4.0d0*dK2(1,ir)*(1.0d0-vrs(ir,izd)**2)*vrs(ir,ixd)&
           +dK2(2,ir)*(2.0d0*vrs(ir,iyd)*(1.0d0-vrs(ir,izd)**2)+4.0d0*vrs(ir,iyd)*vrs(ir,ixd)**2)&
           +4.0d0*dK2(3,ir)*vrs(ir,ixd)& ! new
           +dK2(5,ir)*(4.0d0*(1.0d0-vrs(ir,izd)**2)*vrs(ir,ixd)-16.0d0*vrs(ir,iyd)**2*vrs(ir,ixd))&
           +6.0d0*dK3(1,ir)*(1.0d0-vrs(ir,izd)**2)**2*vrs(ir,ixd)&
           +dK3(2,ir)*(2.0d0*vrs(ir,iyd)*(1.0d0-vrs(ir,izd)**2)**2+8.0d0*vrs(ir,iyd)*vrs(ir,ixd)**2*(1.0d0-vrs(ir,izd)**2))&
           +dK3(5,ir)*(6.0d0*(1.0d0-vrs(ir,izd)**2)*vrs(ir,ixd)-16.0d0*vrs(ir,iyd)**2*vrs(ir,ixd)*(1.0d0-vrs(ir,izd)**2)-16.0d0*vrs(ir,iyd)**2*vrs(ir,ixd)**3)&
           +2.0d0*dK3(6,ir)*(3.0d0*vrs(ir,iyd)*(1.0d0-vrs(ir,izd)**2)**2+12.0d0*vrs(ir,iyd)*vrs(ir,ixd)**2*(1.0d0-vrs(ir,izd)**2)-48.0d0*vrs(ir,iyd)**3*vrs(ir,ixd)**2))
   end do
   
   ! effective field kOe
   do idim=1,3
      HE(idim)=HA(idim)+Hk(idim)/abs(magRST)+Happl(idim)*magRT/magRST
   end do

!   if(ik.eq.2) then
!      write(6,"(6I5,A,3f15.5)") nr,nfe,i,ix,iy,iz,"HA",HA(1),HA(2),HA(3)
!      write(6,"(6I5,A,3f15.5)") nr,nfe,i,ix,iy,iz,"Hk",Hk(1),Hk(2),Hk(3)
!   end if

   ! Torque
   T(1)=vrs(ir,2)*HE(3)-vrs(ir,3)*HE(2)
   T(2)=vrs(ir,3)*HE(1)-vrs(ir,1)*HE(3)
   T(3)=vrs(ir,1)*HE(2)-vrs(ir,2)*HE(1)
   vH=0.0d0 ! inner product
   do idim=1,3
      vH=vH+vrs(ir,idim)*HE(idim)
   end do

   abc=0.0d0
   do idim=1,3 ! make rhs of eq
!      f(ir,idim,ik)=-ratio*gammar*(T(idim)+alphar(ir)*(vH*vr(ir,idim)-HE(idim)))/(1.0d0+alphar(ir)**2)
      f(ir,idim,ik)=-gammars*(-T(idim)+alphar(ir)*(vH*vrs(ir,idim)-HE(idim)))/(1.0d0+alphar(ir)**2)
      abc=abc+T(idim)**2
   end do
   TQ=TQ+dsqrt(abc)
!end do ! i
!end do ! ix
!end do ! iy
!end do ! iz
end do
!write(6,*) "TQ=",TQ


end subroutine fundip_spin2


subroutine fundip_spin(dip,TQ,f,ik,Happl,dK1,dK2,dK3,nunit,nunitx,nunity,nunitz,alphar,gammar,gammars,numint0,&
     nr,nrtot,JRR, magR, iIntrr, numintrr, vr,vrs,irunit,irsub,diprr,ratio,magRT,magRST,&
     nfe,nfetot,JRFe,JRFeS,magFe,iIntrfe,numintrfe,vfe,ifeunit,ifesub,diprfe)
implicit none 
integer(4) :: nr,nfe,numint0,nunit,nunitx,nunity,nunitz,idim,idim2,ik,iIntrr(nr*nunit,numint0),iIntrfe(nr*nunit,numint0)
integer(4) :: ir,ir0,ir00,ix,iy,iz,i,j,ixd,iyd,izd,i2,nrtot,nfetot
integer(4) :: ix2,iy2,iz2,ir2,ife2,ife00
integer(4) :: numintrr(nr*nunit),numintrfe(nr*nunit)
real(8) :: vr(nr*nunit,3),vrs(nr*nunit,3),vfe(nfe*nunit,3),f(nr*nunit,3,4),HE(3),HA(3),Hk(3),Happl(3),T(3)
real(8) :: vH,TQ,abc,def,ghi,gammar,gammars,JRFe,JRFeS,JRR,magR(nr*nunit),magFe(nfe*nunit),magRT,magRST
real(8) :: alphar(nr*nunit),rfept(nfe*nunit,3),rrpt(nr*nunit,3),dis,ratio
real(8) :: fein(nfe*nunit),rin(nr*nunit),dfept(nfe*nunit,3),drpt(nr*nunit,3),delta(3,3)
real(8) :: drfept(-nunitx:nunitx,-nunity:nunity,-nunitz:nunitz,nr,nfe,3),drrpt(-nunitx:nunitx,-nunity:nunity,-nunitz:nunitz,nr,nr,3)
real(8) :: diprr(-nunitx:nunitx,-nunity:nunity,-nunitz:nunitz,nr,nr,3,3),diprfe(-nunitx:nunitx,-nunity:nunity,-nunitz:nunitz,nr,nfe,3,3)
integer(4) :: irunit(nr*nunit,3),ifeunit(nfe*nunit,3),irsub(nr*nunit),ifesub(nfe*nunit)
character(10) :: dip
real(8) :: dK1(7,nr*nunit),dK2(7,nr*nunit),dK3(7,nr*nunit)

delta=0.0d0
do idim=1,3
   delta(idim,idim)=1.0d0
end do

TQ=0.0d0
do ir=1,nrtot
   ix=irunit(ir,1)
   iy=irunit(ir,2)
   iz=irunit(ir,3)
   i=irsub(ir)
   
!   if(magR(ir).eq.0.0d0) cycle
! exchange interaction
   HA=0.0d0
   do idim=1,3
      abc=0.0d0
      do j=1,numintrr(ir)
         abc=abc+vr(iIntrr(ir,j),idim)*abs(magR(iIntrr(ir,j)))
      end do
      def=0.0d0
      do j=1,numintrfe(ir)
         def=def+vfe(iIntrfe(ir,j),idim)*abs(magFe(iIntrfe(ir,j)))
      end do
!      HA(idim)=2.0d0*JRR*abc+JRFe*def ! factor 2*JRR is originate from double counting of interacting pair
      HA(idim)=2.0d0*JRR*abc+JRFeS*def ! factor 2*JRR is originate from double counting of interacting pair
!      write(6,*)  idim,HA(idim)
   end do

   if(dip.eq."on") then
      do idim=1,3
         abc=0.0d0
         do idim2=1,3

            do ir2=1,nrtot
               ix2=irunit(ir2,1)
               iy2=irunit(ir2,2)
               iz2=irunit(ir2,3)
               i2=irsub(ir2)
               abc=abc+diprr(ix-ix2,iy-iy2,iz-iz2,i,i2,idim,idim2)*vr(ir2,idim2)*magR(ir2)
            end do
            do ife2=1,nfetot
               ix2=ifeunit(ife2,1)
               iy2=ifeunit(ife2,2)
               iz2=ifeunit(ife2,3)
               i2=ifesub(ife2)
               abc=abc+diprfe(ix-ix2,iy-iy2,iz-iz2,i,i2,idim,idim2)*vfe(ife2,idim2)*magFe(ife2)
!               write(6,*) i2,ix2,iy2,iz2,abc
            end do

         end do ! idim2
         HA(idim)=HA(idim)-abc*1.0d24 ! angstrome^3 to cm^3
      end do ! idim
   end if ! dipole interaction      

! local anisotropy
   izd=3
   ixd=1
   iyd=2
   Hk=0.0d0
   do ixd=1,2
      if(ixd.eq.1) iyd=2
      if(ixd.eq.2) iyd=1
      Hk(ixd)=-(&
           2.0d0*dK1(1,ir)*vrs(ir,ixd)&
           +2.0d0*dK1(2,ir)*vrs(ir,iyd)&
           +2.0d0*dK1(3,ir)*vrs(ir,ixd)& ! new
           +4.0d0*dK2(1,ir)*(1.0d0-vrs(ir,izd)**2)*vrs(ir,ixd)&
           +dK2(2,ir)*(2.0d0*vrs(ir,iyd)*(1.0d0-vrs(ir,izd)**2)+4.0d0*vrs(ir,iyd)*vrs(ir,ixd)**2)&
           +4.0d0*dK2(3,ir)*vrs(ir,ixd)& ! new
           +dK2(5,ir)*(4.0d0*(1.0d0-vrs(ir,izd)**2)*vrs(ir,ixd)-16.0d0*vrs(ir,iyd)**2*vrs(ir,ixd))&
           +6.0d0*dK3(1,ir)*(1.0d0-vrs(ir,izd)**2)**2*vrs(ir,ixd)&
           +dK3(2,ir)*(2.0d0*vrs(ir,iyd)*(1.0d0-vrs(ir,izd)**2)**2+8.0d0*vrs(ir,iyd)*vrs(ir,ixd)**2*(1.0d0-vrs(ir,izd)**2))&
           +dK3(5,ir)*(6.0d0*(1.0d0-vrs(ir,izd)**2)*vrs(ir,ixd)-16.0d0*vrs(ir,iyd)**2*vrs(ir,ixd)*(1.0d0-vrs(ir,izd)**2)-16.0d0*vrs(ir,iyd)**2*vrs(ir,ixd)**3)&
           +2.0d0*dK3(6,ir)*(3.0d0*vrs(ir,iyd)*(1.0d0-vrs(ir,izd)**2)**2+12.0d0*vrs(ir,iyd)*vrs(ir,ixd)**2*(1.0d0-vrs(ir,izd)**2)-48.0d0*vrs(ir,iyd)**3*vrs(ir,ixd)**2))
   end do
   
   ! effective field kOe
   do idim=1,3
      HE(idim)=HA(idim)+Hk(idim)/abs(magRST)+Happl(idim)*magRT/magRST
   end do

!   if(ik.eq.2) then
!      write(6,"(6I5,A,3f15.5)") nr,nfe,i,ix,iy,iz,"HA",HA(1),HA(2),HA(3)
!      write(6,"(6I5,A,3f15.5)") nr,nfe,i,ix,iy,iz,"Hk",Hk(1),Hk(2),Hk(3)
!   end if

   ! Torque
   T(1)=vrs(ir,2)*HE(3)-vrs(ir,3)*HE(2)
   T(2)=vrs(ir,3)*HE(1)-vrs(ir,1)*HE(3)
   T(3)=vrs(ir,1)*HE(2)-vrs(ir,2)*HE(1)
   vH=0.0d0 ! inner product
   do idim=1,3
      vH=vH+vrs(ir,idim)*HE(idim)
   end do

   abc=0.0d0
   do idim=1,3 ! make rhs of eq
!      f(ir,idim,ik)=-ratio*gammar*(T(idim)+alphar(ir)*(vH*vr(ir,idim)-HE(idim)))/(1.0d0+alphar(ir)**2)
      f(ir,idim,ik)=-gammars*(T(idim)+alphar(ir)*(vH*vrs(ir,idim)-HE(idim)))/(1.0d0+alphar(ir)**2)
      abc=abc+T(idim)**2
   end do
   TQ=TQ+dsqrt(abc)
!end do ! i
!end do ! ix
!end do ! iy
!end do ! iz
end do
!write(6,*) "TQ=",TQ


end subroutine fundip_spin
