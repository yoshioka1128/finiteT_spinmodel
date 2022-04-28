subroutine vacuum(nunitx,nunity,nunitz,nunit,nr,nfe,irpt,ifept,magR,magFe,nunitxbd,nunitybd,nunitzbd)
integer(4) :: nunitx,nunity,nunitz,nr,nfe,ix,iy,i,nunit,irpt(nr,nunitx,nunity,nunitz),ifept(nfe,nunitx,nunity,nunitz),ir,ife
integer(4) :: nunitxbd(2),nunitybd(2),nunitzbd(2),icount
real(8) :: magR(nr*nunit),magFe(nfe*nunit)

! vacuum cell
icount=0
do iz=1,nunitz
   do iy=1,nunity
      do ix=1,nunitx
         if(ix.lt.nunitxbd(1).or.ix.gt.nunitxbd(2).or.iy.lt.nunitybd(1).or.iy.gt.nunitybd(2).or.iz.lt.nunitzbd(1).or.iz.gt.nunitzbd(2)) then
            do i=1,nr
               ir=irpt(i,ix,iy,iz)
               magR(ir)=0.0d0
            end do
            do i=1,nfe
               ife=ifept(i,ix,iy,iz)
               magFe(ife)=0.0d0
            end do
         else
            icount=icount+1
         end if
      end do
   end do
end do

end subroutine vacuum
