!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine position(irpt,ifept,nunit,nunitx,nunity,nunitz,ipx,imx,ipy,imy,ipz,imz,nr,nfe,rfe,rr,ax,cx)
implicit none
integer(4) :: nunit,nr,nfe,nunitx,nunity,nunitz,irpt(nr,nunitx,nunity,nunitz),ifept(nfe,nunitx,nunity,nunitz)
integer(4) :: i,ir,iunit,ipx(nunitx),imx(nunitx),ipy(nunity),imy(nunity),ipz(nunitz),imz(nunitz),ix,iy,iz
integer(4) :: ipt1,ipt2,idim
real(8) :: rfe(nfe,3),rr(nr,3),ax,cx,abc,rfept(nfe*nunit,3),rrpt(nr*nunit,3)

!write(6,*) nunitx,nunity,nunitz
do iz=1,nunitz
   do iy=1,nunity
      do ix=1,nunitx
         do i=1,nr
            irpt(i,ix,iy,iz)=i+(iz-1)*nr+(ix-1)*nr*nunitz+(iy-1)*nr*nunitz*nunitx
!            write(6,*) irpt(i,ix,iy,iz),ix,iy,iz
         end do
         do i=1,nfe
            ifept(i,ix,iy,iz)=i+(iz-1)*nfe+(ix-1)*nfe*nunitz+(iy-1)*nfe*nunitz*nunitx
!            write(6,*) ifept(i,ix,iy,iz),ix,iy,iz 
         end do
      end do
   end do
end do

! shift unitcell
do ix=1,nunitx
   ipx(ix)=ix+1
   imx(ix)=ix-1
end do
do iy=1,nunity
   ipy(iy)=iy+1
   imy(iy)=iy-1
end do
do iz=1,nunitz
   ipz(iz)=iz+1
   imz(iz)=iz-1
end do
! pbc and vacuum for x-y and z direction, respectively
ipx(nunitx)=1
imx(1)=nunitx
ipy(nunity)=1
imy(1)=nunity
ipz(nunitz)=1
imz(1)=nunitz

!write(6,*) "x direction"
!do ix=1,nunitx
!   write(6,*) ix,ipx(ix),imx(ix)
!end do

!write(6,*) "y direction"
!do iy=1,nunity
!   write(6,*) iy,ipy(iy),imx(iy)
!end do

!write(6,*) "z direction"
!do iz=1,nunitz
!   write(6,*) iz,ipz(iz),imz(iz)
!end do
! RE
 ! RE1
rr(1,1)=0.268d0*ax
rr(1,2)=0.268d0*ax
rr(1,3)=1.000d0*cx
 ! RE2
rr(2,1)=0.732d0*ax
rr(2,2)=0.732d0*ax
rr(2,3)=1.000d0*cx
 ! RE3 temp check
rr(3,1)=0.140d0*ax
rr(3,2)=0.860d0*ax
rr(3,3)=1.000d0*cx
 ! RE4 temp check
rr(4,1)=0.860d0*ax
rr(4,2)=0.140d0*ax
rr(4,3)=1.000d0*cx
 ! RE5
rr(5,1)=0.768d0*ax
rr(5,2)=0.232d0*ax
rr(5,3)=0.500d0*cx
 ! RE6
rr(6,1)=0.232d0*ax
rr(6,2)=0.732d0*ax
rr(6,3)=0.500d0*cx
 ! RE7 temp check
rr(7,1)=0.64d0*ax
rr(7,2)=0.640d0*ax
rr(7,3)=0.500d0*cx
 ! RE8 temp check
rr(8,1)=0.360d0*ax
rr(8,2)=0.360d0*ax
rr(8,3)=0.500d0*cx

! Fe
 ! gsite, Fe1
rfe(1,1)=0.223d0*ax
rfe(1,2)=0.567d0*ax
rfe(1,3)=0.127d0*cx
 ! gsite, Fe2
rfe(2,1)=0.777d0*ax
rfe(2,2)=0.433d0*ax
rfe(2,3)=0.873d0*cx
 ! gsite, Fe3
rfe(3,1)=0.933d0*ax
rfe(3,2)=0.723d0*ax
rfe(3,3)=0.627d0*cx
 ! gsite, Fe4
rfe(4,1)=0.067d0*ax
rfe(4,2)=0.277d0*ax
rfe(4,3)=0.373d0*cx
 ! gsite, Fe5
rfe(5,1)=0.277d0*ax
rfe(5,2)=0.067d0*ax
rfe(5,3)=0.627d0*cx

 ! gsite, Fe6
rfe(6,1)=0.723d0*ax
rfe(6,2)=0.933d0*ax
rfe(6,3)=0.373d0*cx
 ! gsite, Fe7
rfe(7,1)=0.777d0*ax
rfe(7,2)=0.433d0*ax
rfe(7,3)=0.127d0*cx
 ! gsite, Fe8
rfe(8,1)=0.223d0*ax
rfe(8,2)=0.567d0*ax
rfe(8,3)=0.873d0*cx
 ! gsite, Fe9
rfe(9,1)=0.433d0*ax
rfe(9,2)=0.777d0*ax
rfe(9,3)=0.127d0*cx
 ! gsite, Fe10
rfe(10,1)=0.567d0*ax
rfe(10,2)=0.223d0*ax
rfe(10,3)=0.873d0*cx

 ! gsite, Fe11
rfe(11,1)=0.067d0*ax
rfe(11,2)=0.277d0*ax
rfe(11,3)=0.627d0*cx
 ! gsite, Fe12
rfe(12,1)=0.933d0*ax
rfe(12,2)=0.723d0*ax
rfe(12,3)=0.373d0*cx
 ! gsite, Fe13
rfe(13,1)=0.723d0*ax
rfe(13,2)=0.933d0*ax
rfe(13,3)=0.627d0*cx
 ! gsite, Fe14
rfe(14,1)=0.277d0*ax
rfe(14,2)=0.067d0*ax
rfe(14,3)=0.373d0*cx
 ! gsite, Fe15
rfe(15,1)=0.567d0*ax
rfe(15,2)=0.223d0*ax
rfe(15,3)=0.127d0*cx

 ! gsite, Fe16
rfe(16,1)=0.433d0*ax
rfe(16,2)=0.777d0*ax
rfe(16,3)=0.873d0*cx
 ! gsite, Fe17
rfe(17,1)=0.037d0*ax
rfe(17,2)=0.360d0*ax
rfe(17,3)=0.176d0*cx
 ! gsite, Fe18
rfe(18,1)=0.963d0*ax
rfe(18,2)=0.640d0*ax
rfe(18,3)=0.824d0*cx
 ! gsite, Fe19
rfe(19,1)=0.140d0*ax
rfe(19,2)=0.537d0*ax
rfe(19,3)=0.676d0*cx
 ! gsite, Fe20
rfe(20,1)=0.860d0*ax
rfe(20,2)=0.463d0*ax
rfe(20,3)=0.324d0*cx

 ! gsite, Fe21
rfe(21,1)=0.463d0*ax
rfe(21,2)=0.860d0*ax
rfe(21,3)=0.676d0*cx
 ! gsite, Fe22
rfe(22,1)=0.537d0*ax
rfe(22,2)=0.140d0*ax
rfe(22,3)=0.324d0*cx
 ! gsite, Fe23
rfe(23,1)=0.963d0*ax
rfe(23,2)=0.640d0*ax
rfe(23,3)=0.176d0*cx
 ! gsite, Fe24
rfe(24,1)=0.037d0*ax
rfe(24,2)=0.360d0*ax
rfe(24,3)=0.824d0*cx
 ! gsite, Fe25
rfe(25,1)=0.640d0*ax
rfe(25,2)=0.963d0*ax
rfe(25,3)=0.176d0*cx

 ! gsite, Fe26
rfe(26,1)=0.360d0*ax
rfe(26,2)=0.037d0*ax
rfe(26,3)=0.824d0*cx
 ! gsite, Fe27
rfe(27,1)=0.860d0*ax
rfe(27,2)=0.463d0*ax
rfe(27,3)=0.676d0*cx
 ! gsite, Fe28
rfe(28,1)=0.140d0*ax
rfe(28,2)=0.537d0*ax
rfe(28,3)=0.324d0*cx
 ! gsite, Fe29
rfe(29,1)=0.537d0*ax
rfe(29,2)=0.140d0*ax
rfe(29,3)=0.676d0*cx
 ! gsite, Fe30
rfe(30,1)=0.463d0*ax
rfe(30,2)=0.860d0*ax
rfe(30,3)=0.324d0*cx

 ! gsite, Fe31
rfe(31,1)=0.360d0*ax
rfe(31,2)=0.037d0*ax
rfe(31,3)=0.176d0*cx
 ! gsite, Fe32
rfe(32,1)=0.640d0*ax
rfe(32,2)=0.963d0*ax
rfe(32,3)=0.824d0*cx
 ! gsite, Fe33
rfe(33,1)=0.098d0*ax
rfe(33,2)=0.098d0*ax
rfe(33,3)=0.204d0*cx
 ! gsite, Fe34
rfe(34,1)=0.902d0*ax
rfe(34,2)=0.902d0*ax
rfe(34,3)=0.796d0*cx
 ! gsite, Fe35
rfe(35,1)=0.402d0*ax
rfe(35,2)=0.598d0*ax
rfe(35,3)=0.704d0*cx

 ! gsite, Fe36
rfe(36,1)=0.598d0*ax
rfe(36,2)=0.402d0*ax
rfe(36,3)=0.296d0*cx
 ! gsite, Fe37
rfe(37,1)=0.902d0*ax
rfe(37,2)=0.902d0*ax
rfe(37,3)=0.204d0*cx
 ! gsite, Fe38
rfe(38,1)=0.098d0*ax
rfe(38,2)=0.098d0*ax
rfe(38,3)=0.796d0*cx
 ! gsite, Fe39
rfe(39,1)=0.598d0*ax
rfe(39,2)=0.402d0*ax
rfe(39,3)=0.704d0*cx
 ! gsite, Fe40
rfe(40,1)=0.402d0*ax
rfe(40,2)=0.598d0*ax
rfe(40,3)=0.296d0*cx

 ! gsite, Fe41
rfe(41,1)=0.317d0*ax
rfe(41,2)=0.317d0*ax
rfe(41,3)=0.246d0*cx
 ! gsite, Fe42
rfe(42,1)=0.683d0*ax
rfe(42,2)=0.683d0*ax
rfe(42,3)=0.754d0*cx
 ! gsite, Fe43
rfe(43,1)=0.183d0*ax
rfe(43,2)=0.817d0*ax
rfe(43,3)=0.746d0*cx
 ! gsite, Fe44
rfe(44,1)=0.817d0*ax
rfe(44,2)=0.183d0*ax
rfe(44,3)=0.254d0*cx
 ! gsite, Fe45
rfe(45,1)=0.683d0*ax
rfe(45,2)=0.683d0*ax
rfe(45,3)=0.246d0*cx

 ! gsite, Fe46
rfe(46,1)=0.317d0*ax
rfe(46,2)=0.317d0*ax
rfe(46,3)=0.754d0*cx
 ! gsite, Fe47
rfe(47,1)=0.817d0*ax
rfe(47,2)=0.183d0*ax
rfe(47,3)=0.746d0*cx
 ! gsite, Fe48
rfe(48,1)=0.183d0*ax
rfe(48,2)=0.817d0*ax
rfe(48,3)=0.254d0*cx
 ! gsite, Fe49
rfe(49,1)=0.500d0*ax
rfe(49,2)=0.500d0*ax
rfe(49,3)=0.114d0*cx
 ! gsite, Fe50
rfe(50,1)=0.500d0*ax
rfe(50,2)=0.500d0*ax
rfe(50,3)=0.886d0*cx

 ! gsite, Fe51
rfe(51,1)=0.000d0*ax
rfe(51,2)=0.000d0*ax
rfe(51,3)=0.614d0*cx
 ! gsite, Fe52
rfe(52,1)=0.000d0*ax
rfe(52,2)=0.000d0*ax
rfe(52,3)=0.386d0*cx
 ! gsite, Fe53 shift
rfe(53,1)=0.000d0*ax
rfe(53,2)=0.500d0*ax
rfe(53,3)=1.000d0*cx
 ! gsite, Fe54
rfe(54,1)=0.000d0*ax
rfe(54,2)=0.500d0*ax
rfe(54,3)=0.500d0*cx
 ! gsite, Fe55
rfe(55,1)=0.500d0*ax
rfe(55,2)=0.000d0*ax
rfe(55,3)=0.500d0*cx

 ! gsite, Fe56 shift
rfe(56,1)=0.500d0*ax
rfe(56,2)=0.000d0*ax
rfe(56,3)=1.000d0*cx


do ipt1=1,nfe
   do iz=-1,1
   do iy=-1,1
   do ix=-1,1
      do ipt2=1,nfe
         abc=0.0d0
         abc=abc+(rr(ipt2,1)+dble(ix)*ax-rr(ipt1,1))**2
         abc=abc+(rr(ipt2,2)+dble(iy)*ax-rr(ipt1,2))**2
         abc=abc+(rr(ipt2,3)+dble(iz)*ax-rr(ipt1,3))**2
         abc=dsqrt(abc)
!         if(abc.lt.4.0d0) write(6,*) ipt1,ipt2,abc
      end do
   end do
   end do
   end do
end do

end subroutine position
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine shift(ax,cx,range,nunit,nunitx,nunity,nunitz,ipx,imx,ipy,imy,ipz,imz,iIntr,numintr,irpt,nr,rr,ifept,nfe,rfe)
implicit none
integer(4) :: nunit,nr,nfe,iIntr(nr*nunit,20),nunitx,nunity,nunitz,irpt(nr,nunitx,nunity,nunitz),ifept(nfe,nunitx,nunity,nunitz)
integer(4) :: i,ir,iunit,ipx(nunitx),imx(nunitx),ipy(nunity),imy(nunity),ipz(nunitz),imz(nunitz)
integer(4) :: ix,iy,iz,ipt1,ipt2,icount,icount0,idx,idy,idz,numintr(nr)
real(8) :: ax,cx,abc,rr(nr,3),rfe(nfe,3),range

do ipt1=1,nr
!   write(6,*) "ipt1=",ipt1
   icount=0
   do idz=-1,1
   do idy=-1,1
   do idx=-1,1
!      write(6,*) idx,idy,idz
      do ipt2=1,nfe
         abc=0.0d0
         abc=abc+(rfe(ipt2,1)+dble(idx)*ax-rr(ipt1,1))**2
         abc=abc+(rfe(ipt2,2)+dble(idy)*ax-rr(ipt1,2))**2
         abc=abc+(rfe(ipt2,3)+dble(idz)*cx-rr(ipt1,3))**2
         abc=dsqrt(abc)
         if(abc.lt.range.and.abc.ne.0.0d0) then
            icount=icount+1
!            write(6,*) icount,ipt2,abc
!            write(6,*) ipt1,ipt2,ix,iy,idz

            do ix=1,nunitx
            do iy=1,nunity
            do iz=1,nunitz
!            write(6,*) ipt1,ix,iy,iz
!            write(6,*) ipt2,ix+idx,iy+idy,iz+idz

               ir=irpt(ipt1,ix,iy,iz)
            if(idx.eq.-1.and.idy.eq.-1.and.idz.eq.-1) iIntr(ir,icount)=ifept(ipt2,imx(ix),imy(iy),imz(iz))
            if(idx.eq.0.and.idy.eq.-1.and.idz.eq.-1)  iIntr(ir,icount)=ifept(ipt2,ix,imy(iy),imz(iz))
            if(idx.eq.1.and.idy.eq.-1.and.idz.eq.-1)  iIntr(ir,icount)=ifept(ipt2,ipx(ix),imy(iy),imz(iz))

            if(idx.eq.-1.and.idy.eq.0.and.idz.eq.-1) iIntr(ir,icount)=ifept(ipt2,imx(ix),iy,imz(iz))
            if(idx.eq.0.and.idy.eq.0.and.idz.eq.-1)  iIntr(ir,icount)=ifept(ipt2,ix,iy,imz(iz))
            if(idx.eq.1.and.idy.eq.0.and.idz.eq.-1)  iIntr(ir,icount)=ifept(ipt2,ipx(ix),iy,imz(iz))

            if(idx.eq.-1.and.idy.eq.1.and.idz.eq.-1) iIntr(ir,icount)=ifept(ipt2,imx(ix),ipy(iy),imz(iz))
            if(idx.eq.0.and.idy.eq.1.and.idz.eq.-1)  iIntr(ir,icount)=ifept(ipt2,ix,ipy(iy),imz(iz))
            if(idx.eq.1.and.idy.eq.1.and.idz.eq.-1)  iIntr(ir,icount)=ifept(ipt2,ipx(ix),ipy(iy),imz(iz))

            if(idx.eq.-1.and.idy.eq.-1.and.idz.eq.0) iIntr(ir,icount)=ifept(ipt2,imx(ix),imy(iy),iz)
            if(idx.eq.0.and.idy.eq.-1.and.idz.eq.0)  iIntr(ir,icount)=ifept(ipt2,ix,imy(iy),iz)
            if(idx.eq.1.and.idy.eq.-1.and.idz.eq.0)  iIntr(ir,icount)=ifept(ipt2,ipx(ix),imy(iy),iz)

            if(idx.eq.-1.and.idy.eq.0.and.idz.eq.0) iIntr(ir,icount)=ifept(ipt2,imx(ix),iy,iz)
            if(idx.eq.0.and.idy.eq.0.and.idz.eq.0)  iIntr(ir,icount)=ifept(ipt2,ix,iy,iz)
            if(idx.eq.1.and.idy.eq.0.and.idz.eq.0)  iIntr(ir,icount)=ifept(ipt2,ipx(ix),iy,iz)

            if(idx.eq.-1.and.idy.eq.1.and.idz.eq.0) iIntr(ir,icount)=ifept(ipt2,imx(ix),ipy(iy),iz)
            if(idx.eq.0.and.idy.eq.1.and.idz.eq.0)  iIntr(ir,icount)=ifept(ipt2,ix,ipy(iy),iz)
            if(idx.eq.1.and.idy.eq.1.and.idz.eq.0)  iIntr(ir,icount)=ifept(ipt2,ipx(ix),ipy(iy),iz)

            if(idx.eq.-1.and.idy.eq.-1.and.idz.eq.1) iIntr(ir,icount)=ifept(ipt2,imx(ix),imy(iy),ipz(iz))
            if(idx.eq.0.and.idy.eq.-1.and.idz.eq.1)  iIntr(ir,icount)=ifept(ipt2,ix,imy(iy),ipz(iz))
            if(idx.eq.1.and.idy.eq.-1.and.idz.eq.1)  iIntr(ir,icount)=ifept(ipt2,ipx(ix),imy(iy),ipz(iz))

            if(idx.eq.-1.and.idy.eq.0.and.idz.eq.1) iIntr(ir,icount)=ifept(ipt2,imx(ix),iy,ipz(iz))
            if(idx.eq.0.and.idy.eq.0.and.idz.eq.1)  iIntr(ir,icount)=ifept(ipt2,ix,iy,ipz(iz))
            if(idx.eq.1.and.idy.eq.0.and.idz.eq.1)  iIntr(ir,icount)=ifept(ipt2,ipx(ix),iy,ipz(iz))

            if(idx.eq.-1.and.idy.eq.1.and.idz.eq.1) iIntr(ir,icount)=ifept(ipt2,imx(ix),ipy(iy),ipz(iz))
            if(idx.eq.0.and.idy.eq.1.and.idz.eq.1)  iIntr(ir,icount)=ifept(ipt2,ix,ipy(iy),ipz(iz))
            if(idx.eq.1.and.idy.eq.1.and.idz.eq.1)  iIntr(ir,icount)=ifept(ipt2,ipx(ix),ipy(iy),ipz(iz))
         end do
         end do
         end do
         end if

      end do
   end do
   end do
   end do
   numintr(ipt1)=icount
end do

end subroutine shift
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
