subroutine engdens(engr,engrave,engr0,icountr,Happl,Klc,nunit,nunitx,nunity,nunitz,irpt,&
     nr, JRR, magR, iIntrr, numintrr, vr,&
     nfe,JRFe,magFe,iIntrfe,numintrfe,vfe)
implicit none 
integer(4) :: nr,nfe,nunit,nunitx,nunity,nunitz,idim,ik,iIntrr(nr*nunit,20),iIntrfe(nr*nunit,20),ir,ix,iy,iz,iunit,i,j,ixd,iyd,izd
integer(4) :: irpt(nr,nunitx,nunity,nunitz),numintrr(nr),numintrfe(nr)
real(8) :: vr(nr*nunit,3),vfe(nfe*nunit,3),Happl(3),T(3)
real(8) :: vH,abc,def,JRFe,JRR,Klc(9,nr*nunit),magR(nr*nunit),magFe(nfe*nunit),engr(nr*nunit,4),engrave(nr),engr0(nr,4)
integer(4) :: ife,icountr
real(8),parameter :: muB=9.274d-24 ! (Merg/kOe) or (kerg/Oe)
real(8),parameter :: mev=1.60218d-21 ! (Merg/mev) 
real(8),parameter :: kb=8.6173303d-2*mev ! (Merg/K)

engr=0.0d0
do iz=1,nunitz
write(6,*) iz
do iy=1,nunity
do ix=1,nunitx
do i=1,nr
   ir=irpt(i,ix,iy,iz)
   if(dabs(vr(ir,1)**2+vr(ir,2)**2+vr(ir,3)**2).eq.0.0d0) cycle

      ! exchange energy
      do idim=1,3
         abc=0.0d0
         do j=1,numintrr(i)
            abc=abc+vr(iIntrr(ir,j),idim)*magR(iIntrr(ir,j))
        end do
         def=0.0d0
         do j=1,numintrfe(i)
            def=def+vfe(iIntrfe(ir,j),idim)*magFe(iIntrfe(ir,j))
         end do
         engr(ir,1)=engr(ir,1)-(JRR*abc+JRFe*def)*vr(ir,idim)*magR(ir)
      end do

      ! anisotropy energy
      engr(ir,2)=(&
           (Klc(1,ir)+Klc(2,ir)*(vr(ir,1)**2+vr(ir,2)**2)+Klc(4,ir)*(vr(ir,1)**2+vr(ir,2)**2)**2)*(vr(ir,1)**2+vr(ir,2)**2)&
           +(Klc(6,ir)+Klc(7,ir)*(vr(ir,1)**2+vr(ir,2)**2)+Klc(8,ir)*(vr(ir,1)**2+vr(ir,2)**2)**2)*2.0d0*vr(ir,1)*vr(ir,2)&
           +(Klc(3,ir)+Klc(5,ir)*(vr(ir,1)**2+vr(ir,2)**2))*((vr(ir,1)**2+vr(ir,2)**2)**2-8.0d0*vr(ir,1)**2*vr(ir,2)**2)&
           +Klc(9,ir)*2.0d0*(3.0d0*vr(ir,1)*vr(ir,2)*(vr(ir,1)**2+vr(ir,2)**2)**2-16.0d0*vr(ir,1)**3*vr(ir,2)**3))-engr0(i,2)

      ! Zeeman energy
      do idim=1,3
         engr(ir,3)=engr(ir,2)-vr(ir,idim)*Happl(idim)*magR(ir)
      end do
      if(i.eq.17) write(6,*) engr(ir,1),engr(ir,1)/kb
      
      engr(ir,1)=engr(ir,1)-engr0(i,1) ! exchange
      engr(ir,3)=engr(ir,3)-engr0(i,3) ! Zeeman
      engr(ir,4)=engr(ir,1)+engr(ir,2)+engr(ir,3) ! total
end do ! for i=1,nr
end do ! for ix
end do ! for iy
end do ! for iz
engrave=0.0d0
do i=1,4
   do ir=1,nr*nunit
      engrave(i)=engrave(i)+engr(ir,i)
   end do
end do
engrave=engrave/dble(icountr)

end subroutine engdens
