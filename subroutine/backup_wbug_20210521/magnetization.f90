subroutine magnetization2(numr,numfe,nr,nfe,nfu,nunit,vr,vfe,magR,magFe,magtot,magunit,nunitx,nunity,nunitz,iunit,irunit,ifeunit,&
     nrtot,nfetot,nunit0)
implicit none
integer(4) :: numr,numfe,ir,ife,nr,nfe,nunit,nunit0,i,j,nunitx,nunity,nunitz,iunit(nunitx,nunity,nunitz)
integer(4) :: idim,irunit(nr*nunit,3),ifeunit(nfe*nunit,3),nrtot,nfetot,ir0,ife0
integer(4) :: ix,iy,iz,nux,nuy,nuz,nfu
real(8) :: vr(nr*nunit,3),vfe(nfe*nunit,3),magR(nunit*nr),magFe(nunit*nfe)
real(8) :: magtot(3),magunit(nunit,3)
character(10) :: pbcx,pbcy,pbcz

! magnetization calc.
magunit=0.0d0
do idim=1,3
   do ir=1,nrtot
      ix=irunit(ir,1)
      iy=irunit(ir,2)
      iz=irunit(ir,3)
      magunit(iunit(ix,iy,iz),idim)=magunit(iunit(ix,iy,iz),idim)+vr(ir,idim)*abs(magR(ir))
   end do
   do ife=1,nfetot
      ix=ifeunit(ife,1)
      iy=ifeunit(ife,2)
      iz=ifeunit(ife,3)
      magunit(iunit(ix,iy,iz),idim)=magunit(iunit(ix,iy,iz),idim)+vfe(ife,idim)*abs(magFe(ife))
   end do
end do

magunit=magunit/dble(nfu) ! per f.u.
!magunit=magunit/dble(nfu) ! per f.u.
magtot=0.0d0
do idim=1,3
   do i=1,nunit
      magtot(idim)=magtot(idim)+magunit(i,idim)
   end do
end do

magtot=magtot/dble(nunit0) 
!write(6,*) "magtot",magtot

end subroutine magnetization2




subroutine magnetization(numr,numfe,nr,nfe,nunit,vr,vfe,magR,magFe,magtot,magunit,nunitx,nunity,nunitz,iunit,irunit,ifeunit,&
     nrtot,nfetot,nunit0)
implicit none
integer(4) :: numr,numfe,ir,ife,nr,nfe,nunit,nunit0,i,j,nunitx,nunity,nunitz,iunit(nunitx,nunity,nunitz)
integer(4) :: idim,irunit(nr*nunit,3),ifeunit(nfe*nunit,3),nrtot,nfetot,ir0,ife0
integer(4) :: ix,iy,iz,nux,nuy,nuz
real(8) :: vr(nr*nunit,3),vfe(nfe*nunit,3),magR(nunit*nr),magFe(nunit*nfe)
real(8) :: magtot(3),magunit(nunit,3)
character(10) :: pbcx,pbcy,pbcz

! magnetization calc.
magunit=0.0d0
do idim=1,3
   do ir=1,nrtot
      ix=irunit(ir,1)
      iy=irunit(ir,2)
      iz=irunit(ir,3)
      magunit(iunit(ix,iy,iz),idim)=magunit(iunit(ix,iy,iz),idim)+vr(ir,idim)*magR(ir)
   end do
   do ife=1,nfetot
      ix=ifeunit(ife,1)
      iy=ifeunit(ife,2)
      iz=ifeunit(ife,3)
      magunit(iunit(ix,iy,iz),idim)=magunit(iunit(ix,iy,iz),idim)+vfe(ife,idim)*magFe(ife)
   end do
end do

magunit=magunit/4.0d0 ! for f.u.
magtot=0.0d0
do idim=1,3
   do i=1,nunit
      magtot(idim)=magtot(idim)+magunit(i,idim)
   end do
end do

magtot=magtot/dble(nunit0) 
!write(6,*) "magtot",magtot

end subroutine magnetization
