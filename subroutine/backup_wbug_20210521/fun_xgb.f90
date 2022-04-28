subroutine fun(TQ,f,ik,HE,HA,Hk,Happl,Klc,nr,nfe,nunit,nunitx,nunity,nunitz,alpha,alphagb,gamma,&
     JRR,magR,iIntrr,numintrr,vr,irpt,&
     JRFe,magFe,iIntrfe,numintrfe,vfe,ifept,&
     igb)

implicit none 
integer(4) :: nr,nfe,nunit,nunitx,nunity,nunitz,idim,ik,iIntrr(nr*nunit,20),iIntrfe(nr*nunit,20),ir,ix,iy,iz,iunit,i,j,ixd,iyd,izd
integer(4) :: irpt(nr,nunitx,nunity,nunitz),numintrr(nr),numintrfe(nr),ifept(nfe,nunitx,nunity,nunitz)
real(8) :: vr(nr*nunit,3),vfe(nfe*nunit,3),f(nr*nunit,3,4),HE(nr*nunit,3),HA(nr*nunit,3),Hk(nr*nunit,3),Happl(3),T(3)
real(8) :: vH,alpha,alphagb,TQ,abc,def,gamma,JRFe,JRR,Klc(9,nr*nunit),magR(nr*nunit),magFe(nfe*nunit)
integer(4) :: igb(2),ife

TQ=0.0d0
Hk=0.0d0
HA=0.0d0
do iz=1,nunitz
do iy=1,nunity
do ix=1,nunitx
do i=1,nr
   ir=irpt(i,ix,iy,iz)
   ife=ifept(i,ix,iy,iz)
   if(dabs(vr(ir,1)**2+vr(ir,2)**2+vr(ir,3)**2).eq.0.0d0) cycle
      do idim=1,3
         abc=0.0d0
         do j=1,numintrr(i)
            abc=abc+vr(iIntrr(ir,j),idim)
         end do
         def=0.0d0
         do j=1,numintrfe(i)
            def=def+vfe(iIntrfe(ir,j),idim)
         end do
         HA(ir,idim)=JRR*abc*magR(ir)+JRFe*def*magFe(ife)
      end do

      izd=3
      ixd=1
      iyd=2
      ir=irpt(i,ix,iy,iz)
      do ixd=1,2
         if(ixd.eq.1) iyd=2
         if(ixd.eq.2) iyd=1
         Hk(ir,ixd)=-(&
         2.0d0*Klc(1,ir)*vr(ir,ixd)&
        +4.0d0*Klc(2,ir)*(1.0d0-vr(ir,izd)**2)*vr(ir,ixd)&
        +6.0d0*Klc(4,ir)*(1.0d0-vr(ir,izd)**2)**2*vr(ir,ixd)&
        +2.0d0*Klc(6,ir)*vr(ir,iyd)&
        +Klc(7,ir)*(2.0d0*vr(ir,iyd)*(1.0d0-vr(ir,izd)**2)+4.0d0*vr(ir,iyd)*vr(ir,ixd)**2)&
        +Klc(8,ir)*(2.0d0*vr(ir,iyd)*(1.0d0-vr(ir,izd)**2)**2+8.0d0*vr(ir,iyd)*vr(ir,ixd)**2*(1.0d0-vr(ir,izd)**2))&
        +Klc(3,ir)*(4.0d0*(1.0d0-vr(ir,izd)**2)*vr(ir,ixd)-16.0d0*vr(ir,iyd)**2*vr(ir,ixd))&
        +Klc(5,ir)*(6.0d0*(1.0d0-vr(ir,izd)**2)*vr(ir,ixd)-16.0d0*vr(ir,iyd)**2*vr(ir,ixd)*(1.0d0-vr(ir,izd)**2)-16.0d0*vr(ir,iyd)**2*vr(ir,ixd)**3)&
        +2.0d0*Klc(9,ir)*(3.0d0*vr(ir,iyd)*(1.0d0-vr(ir,izd)**2)**2+12.0d0*vr(ir,iyd)*vr(ir,ixd)**2*(1.0d0-vr(ir,izd)**2)-48.0d0*vr(ir,iyd)**3*vr(ir,ixd)**2))
      end do

   ! effective field
      ir=irpt(i,ix,iy,iz)
      if(ix.ge.igb(1).and.ix.le.igb(2)) then
         do idim=1,3 ! for gb phase
            HE(ir,idim)=HA(ir,idim)+Happl(idim) ! no anisotropy field
         end do
      else
         do idim=1,3
            HE(ir,idim)=HA(ir,idim)+Hk(ir,idim)/magR(ir)+Happl(idim)
         end do
      end if

   ! Torque
      ir=irpt(i,ix,iy,iz)
      T(1)=vr(ir,2)*HE(ir,3)-vr(ir,3)*HE(ir,2)
      T(2)=vr(ir,3)*HE(ir,1)-vr(ir,1)*HE(ir,3)
      T(3)=vr(ir,1)*HE(ir,2)-vr(ir,2)*HE(ir,1)
      vH=0.0d0 ! inner product
      do idim=1,3
         vH=vH+vr(ir,idim)*HE(ir,idim)
      end do

      abc=0.0d0
      do idim=1,3 ! make rhs of eq
         if(ix.ge.igb(1).and.ix.le.igb(2)) then
            f(ir,idim,ik)=-gamma*(T(idim)+alphagb*(vH*vr(ir,idim)-HE(ir,idim)))/(1.0d0+alpha**2)
            abc=abc+T(idim)**2
         else
            f(ir,idim,ik)=-gamma*(T(idim)+alpha*(vH*vr(ir,idim)-HE(ir,idim)))/(1.0d0+alpha**2)
            abc=abc+T(idim)**2
         end if
      end do
      TQ=TQ+dsqrt(abc)

end do ! for i
end do ! for ix
end do ! for iy
end do ! for iz

end subroutine fun
