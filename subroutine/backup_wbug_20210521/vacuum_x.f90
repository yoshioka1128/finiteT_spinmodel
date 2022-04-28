subroutine vacuum_x(nunitx,nunity,nunitz,nunit,nr,nfe,vr,vfe,irpt,ifept,magR,magFe,nunitxbd,nunitybd,nunitzbd)
integer(4) :: nunitx,nunity,nunitz,nr,nfe,ix,iy,i,nunit,irpt(nr,nunitx,nunity,nunitz),ifept(nfe,nunitx,nunity,nunitz),ixv1,ixv2
integer(4) :: nunitxbd(2),nunitybd(2),nunitzbd(2)
real(8) :: vr(nr*nunit,3),vfe(nfe*nunit,3),magR(nr*nunit),magFe(nfe*nunit)

! vacuum cell
do iz=1,nunitz
   do iy=1,nunity
      do ix=1,nunitx
         if(ix.lt.nunitxbd(1).or.ix.gt.nunitxbd(2).or.iy.lt.nunitybd(1).or.iy.gt.nunitybd(2).or.iz.lt.nunitzbd(1).or.iz.gt.nunitzbd(2)) then
            do i=1,nr
               ixv1=irpt(i,1,iy,iz)
               ixv2=irpt(i,nunitx,iy,iz)
               magR(ixv1)=0.0d0
               magR(ixv2)=0.0d0
            end do
            do i=1,nfe
               ixv1=ifept(i,1,iy,iz)
               ixv2=ifept(i,nunitx,iy,iz)
               magFe(ixv1)=0.0d0
               magFe(ixv2)=0.0d0
            end do
         end if
      end do
   end do
end do

end subroutine vacuum_x
