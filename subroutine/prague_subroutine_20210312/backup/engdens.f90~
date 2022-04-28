subroutine engdens(engr,engrunit,engrtot,engr0,numr,Happl,Klc,nunit,nunitx,nunity,nunitz,irunit,iunit,&
     nr, JRR, magR, iIntrr, numintrr, vr,nrtot,&
     nfe,JRFe,magFe,iIntrfe,numintrfe,vfe,nfetot,nunit0)
implicit none 
integer(4) :: nr,nfe,nunit,nunitx,nunity,nunitz,idim,ik,iIntrr(nr*nunit,20),iIntrfe(nr*nunit,20),ir,ix,iy,iz,iunit(nunitx,nunity,nunitz),i,j,ixd,iyd,izd
integer(4) :: irpt(nr,nunitx,nunity,nunitz),numintrr(nr),numintrfe(nr),irunit(nr*nunit,3),nunit0
real(8) :: vr(nr*nunit,3),vfe(nfe*nunit,3),Happl(3),T(3)
real(8) :: vH,abc,def,JRFe,JRR,Klc(9,nr*nunit),magR(nr*nunit),magFe(nfe*nunit),engr(nr*nunit,4),engrunit(nunit,4),engr0(nr*nunit,4),engrtot(4)
integer(4) :: ife,numr,inum,nrtot,nfetot,ife0,ir0
real(8),parameter :: muB=9.274d-24 ! (Merg/kOe) or (kerg/Oe)
real(8),parameter :: mev=1.60218d-21 ! (Merg/mev) 
real(8),parameter :: kb=8.6173303d-2*mev ! (Merg/K)

engr=0.0d0
do ir=1,nrtot
   
   if(magR(ir).eq.0.0d0) cycle
      ! exchange energy
      do idim=1,3
         abc=0.0d0
         do j=1,numintrr(ir)
            abc=abc+vr(iIntrr(ir,j),idim)*magR(iIntrr(ir,j))
        end do
         def=0.0d0
         do j=1,numintrfe(ir)
            def=def+vfe(iIntrfe(ir,j),idim)*magFe(iIntrfe(ir,j))
         end do
         engr(ir,1)=engr(ir,1)-(2.0d0*JRR*abc+JRFe*def)*vr(ir,idim)
      end do
      engr(ir,1)=engr(ir,1)*magR(ir)

      ! anisotropy energy
      engr(ir,2)=(&
           (Klc(1,ir)+Klc(2,ir)*(vr(ir,1)**2+vr(ir,2)**2)+Klc(4,ir)*(vr(ir,1)**2+vr(ir,2)**2)**2)*(vr(ir,1)**2+vr(ir,2)**2)&
           +(Klc(6,ir)+Klc(7,ir)*(vr(ir,1)**2+vr(ir,2)**2)+Klc(8,ir)*(vr(ir,1)**2+vr(ir,2)**2)**2)*2.0d0*vr(ir,1)*vr(ir,2)&
           +(Klc(3,ir)+Klc(5,ir)*(vr(ir,1)**2+vr(ir,2)**2))*((vr(ir,1)**2+vr(ir,2)**2)**2-8.0d0*vr(ir,1)**2*vr(ir,2)**2)&
           +Klc(9,ir)*2.0d0*(3.0d0*vr(ir,1)*vr(ir,2)*(vr(ir,1)**2+vr(ir,2)**2)**2-16.0d0*vr(ir,1)**3*vr(ir,2)**3))-engr0(ir,2)

      ! Zeeman energy
      do idim=1,3
         engr(ir,3)=engr(ir,3)-vr(ir,idim)*Happl(idim)*magR(ir)
      end do
      
      engr(ir,1)=engr(ir,1)-engr0(ir,1) ! exchange
      engr(ir,3)=engr(ir,3)-engr0(ir,3) ! Zeeman
      engr(ir,4)=engr(ir,1)+engr(ir,2)+engr(ir,3) ! total
end do ! for i=1,nr

! local energy per unit cell
engrunit=0.0d0
do ir=1,nrtot
   do inum=1,4
      ix=irunit(ir,1)
      iy=irunit(ir,2)
      iz=irunit(ir,3)
      engrunit(iunit(ix,iy,iz),inum)=engrunit(iunit(ix,iy,iz),inum)+engr(ir,inum) 
   end do
end do
  
! total energy per unit cell
engrtot=0.0d0
do inum=1,4
   do i=1,nunit
      engrtot(inum)=engrtot(inum)+engrunit(i,inum)
   end do
end do
engrtot=engrtot/dble(nunit0)
!write(6,*) "total eng",engrtot(4)

end subroutine engdens
