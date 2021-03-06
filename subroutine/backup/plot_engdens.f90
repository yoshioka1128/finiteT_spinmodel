subroutine plot_engdens(itime,nunit,nr,nfe,engr,engfe,nunitx,nunity,nunitz,rr,rfe,irpt,ifept,&
     magR,magFe,ax,cx,nplotx,nploty,nplotz)
implicit none
integer(4) :: nunit,nr,nfe,itime,iunit,i,nunitx,nunity,nunitz,ix,iy,iz,ir,irpt(nr,nunitx,nunity,nunitz),ifept(nfe,nunitx,nunity,nunitz)
real(8) :: magR(nr*nunit),magFe(nfe*nunit)
integer(4) :: ipx(nunitx),ipy(nunity),ipz(nunitz),imx(nunitx),imy(nunity),imz(nunitz),nplotx(2),nploty(2),nplotz(2)
real(8) :: engr(nr*nunit,4),engfe(nfe*nunit,4),rr(nr,3),rfe(nfe,3),ax,cx
character(5) :: time

write(time,'(I5)') itime
open(108,file="movie_source/Esum_TDf"//trim(adjustl(time))//".txt") 
open(109,file="movie_source/Esum_TDg"//trim(adjustl(time))//".txt") 
open(110,file="movie_source/Esum_TDFe"//trim(adjustl(time))//".txt") 

open(111,file="movie_source/Eex_TDf"//trim(adjustl(time))//".txt") 
open(112,file="movie_source/Eex_TDg"//trim(adjustl(time))//".txt") 
open(113,file="movie_source/Eex_TDFe"//trim(adjustl(time))//".txt") 

open(114,file="movie_source/EA_TDf"//trim(adjustl(time))//".txt") 
open(115,file="movie_source/EA_TDg"//trim(adjustl(time))//".txt") 
open(116,file="movie_source/EA_TDFe"//trim(adjustl(time))//".txt") 

open(117,file="movie_source/EH_TDf"//trim(adjustl(time))//".txt") 
open(118,file="movie_source/EH_TDg"//trim(adjustl(time))//".txt") 
open(119,file="movie_source/EH_TDFe"//trim(adjustl(time))//".txt") 

do iz=nplotz(1),nplotz(2)
do iy=nploty(1),nploty(2)
do ix=nplotx(1),nplotx(2)
!do ix=nunitx-nplotx,nunitx-1
   do ir=1,2
      ! inside
      ! fsite
      if(magR(irpt(ir,ix,iy,iz)).ne.0.0d0) then
      write(108,"(4E20.10e2)") rr(ir,1)+dble(ix-nplotx(1))*ax,rr(ir,2)+dble(iy-nploty(1))*ax,rr(ir,3)+dble(iz-nplotz(1))*cx,&
           engr(irpt(ir,ix,iy,iz),4)
      write(111,"(4E20.10e2)") rr(ir,1)+dble(ix-nplotx(1))*ax,rr(ir,2)+dble(iy-nploty(1))*ax,rr(ir,3)+dble(iz-nplotz(1))*cx,&
           engr(irpt(ir,ix,iy,iz),1)
      write(114,"(4E20.10e2)") rr(ir,1)+dble(ix-nplotx(1))*ax,rr(ir,2)+dble(iy-nploty(1))*ax,rr(ir,3)+dble(iz-nplotz(1))*cx,&
           engr(irpt(ir,ix,iy,iz),2)
      write(117,"(4E20.10e2)") rr(ir,1)+dble(ix-nplotx(1))*ax,rr(ir,2)+dble(iy-nploty(1))*ax,rr(ir,3)+dble(iz-nplotz(1))*cx,&
           engr(irpt(ir,ix,iy,iz),3)
      end if
      if(magR(irpt(ir+4,ix,iy,iz)).ne.0.0d0) then
      write(108,"(4E20.10e2)") rr(ir+4,1)+dble(ix-nplotx(1))*ax,rr(ir+4,2)+dble(iy-nploty(1))*ax,rr(ir+4,3)+dble(iz-nplotz(1))*cx,&
           engr(irpt(ir+4,ix,iy,iz),4)
      write(111,"(4E20.10e2)") rr(ir+4,1)+dble(ix-nplotx(1))*ax,rr(ir+4,2)+dble(iy-nploty(1))*ax,rr(ir+4,3)+dble(iz-nplotz(1))*cx,&
           engr(irpt(ir+4,ix,iy,iz),1)
      write(114,"(4E20.10e2)") rr(ir+4,1)+dble(ix-nplotx(1))*ax,rr(ir+4,2)+dble(iy-nploty(1))*ax,rr(ir+4,3)+dble(iz-nplotz(1))*cx,&
           engr(irpt(ir+4,ix,iy,iz),2)
      write(117,"(4E20.10e2)") rr(ir+4,1)+dble(ix-nplotx(1))*ax,rr(ir+4,2)+dble(iy-nploty(1))*ax,rr(ir+4,3)+dble(iz-nplotz(1))*cx,&
           engr(irpt(ir+4,ix,iy,iz),3)
      end if
      ! gsite
      if(magR(irpt(ir+2,ix,iy,iz)).ne.0.0d0) then
      write(109,"(4E20.10e2)") rr(ir+2,1)+dble(ix-nplotx(1))*ax,rr(ir+2,2)+dble(iy-nploty(1))*ax,rr(ir+2,3)+dble(iz-nplotz(1))*cx,&
           engr(irpt(ir+2,ix,iy,iz),4)
      write(112,"(4E20.10e2)") rr(ir+2,1)+dble(ix-nplotx(1))*ax,rr(ir+2,2)+dble(iy-nploty(1))*ax,rr(ir+2,3)+dble(iz-nplotz(1))*cx,&
           engr(irpt(ir+2,ix,iy,iz),1)
      write(115,"(4E20.10e2)") rr(ir+2,1)+dble(ix-nplotx(1))*ax,rr(ir+2,2)+dble(iy-nploty(1))*ax,rr(ir+2,3)+dble(iz-nplotz(1))*cx,&
           engr(irpt(ir+2,ix,iy,iz),2)
      write(118,"(4E20.10e2)") rr(ir+2,1)+dble(ix-nplotx(1))*ax,rr(ir+2,2)+dble(iy-nploty(1))*ax,rr(ir+2,3)+dble(iz-nplotz(1))*cx,&
           engr(irpt(ir+2,ix,iy,iz),3)
      end if
      if(magR(irpt(ir+6,ix,iy,iz)).ne.0.0d0) then
      write(109,"(4E20.10e2)") rr(ir+6,1)+dble(ix-nplotx(1))*ax,rr(ir+6,2)+dble(iy-nploty(1))*ax,rr(ir+6,3)+dble(iz-nplotz(1))*cx,&
           engr(irpt(ir+6,ix,iy,iz),4)
      write(112,"(4E20.10e2)") rr(ir+6,1)+dble(ix-nplotx(1))*ax,rr(ir+6,2)+dble(iy-nploty(1))*ax,rr(ir+6,3)+dble(iz-nplotz(1))*cx,&
           engr(irpt(ir+6,ix,iy,iz),1)
      write(115,"(4E20.10e2)") rr(ir+6,1)+dble(ix-nplotx(1))*ax,rr(ir+6,2)+dble(iy-nploty(1))*ax,rr(ir+6,3)+dble(iz-nplotz(1))*cx,&
           engr(irpt(ir+6,ix,iy,iz),2)
      write(118,"(4E20.10e2)") rr(ir+6,1)+dble(ix-nplotx(1))*ax,rr(ir+6,2)+dble(iy-nploty(1))*ax,rr(ir+6,3)+dble(iz-nplotz(1))*cx,&
           engr(irpt(ir+6,ix,iy,iz),3)
      end if
   end do
end do
end do
end do

do iz=nplotz(1),nplotz(2)
do iy=nploty(1),nploty(2)
do ix=nplotx(1),nplotx(2)
!do ix=nunitx-nplotx,nunitx-1
   ! inside
   do ir=1,nfe
      if(magFe(ifept(ir,ix,iy,iz)).ne.0.0d0) then
      write(110,"(4E20.10e2)") rfe(ir,1)+dble(ix-nplotx(1))*ax,rfe(ir,2)+dble(iy-nploty(1))*ax,rfe(ir,3)+dble(iz-nplotz(1))*cx,&
           engfe(ifept(ir,ix,iy,iz),4)
      write(113,"(4E20.10e2)") rfe(ir,1)+dble(ix-nplotx(1))*ax,rfe(ir,2)+dble(iy-nploty(1))*ax,rfe(ir,3)+dble(iz-nplotz(1))*cx,&
           engfe(ifept(ir,ix,iy,iz),1)
      write(116,"(4E20.10e2)") rfe(ir,1)+dble(ix-nplotx(1))*ax,rfe(ir,2)+dble(iy-nploty(1))*ax,rfe(ir,3)+dble(iz-nplotz(1))*cx,&
           engfe(ifept(ir,ix,iy,iz),2)
      write(119,"(4E20.10e2)") rfe(ir,1)+dble(ix-nplotx(1))*ax,rfe(ir,2)+dble(iy-nploty(1))*ax,rfe(ir,3)+dble(iz-nplotz(1))*cx,&
           engfe(ifept(ir,ix,iy,iz),3)
      end if
   end do
end do
end do
end do

!! bottom
!do iy=nploty(1),nploty(2)
!do ix=nplotx(1),nplotx(2)
!iz=nplotz(1)
!!do ix=nunitx-nplotx,nunitx-1
!   do ir=1,2
!      ! fsite
!      if(magR(irpt(ir,ix,iy,imz(iz))).ne.0.0d0) then
!      write(108,"(4E20.10e2)") rr(ir,1)+dble(ix-nplotx(2))*ax,rr(ir,2)+dble(iy-nploty(2))*ax,rr(ir,3)+dble(iz-nplotz(2)-1)*cx,&
!           engr(irpt(ir,ix,iy,imz(iz)),4)
!      write(111,"(4E20.10e2)") rr(ir,1)+dble(ix-nplotx(2))*ax,rr(ir,2)+dble(iy-nploty(2))*ax,rr(ir,3)+dble(iz-nplotz(2)-1)*cx,&
!           engr(irpt(ir,ix,iy,imz(iz)),1)
!      write(114,"(4E20.10e2)") rr(ir,1)+dble(ix-nplotx(2))*ax,rr(ir,2)+dble(iy-nploty(2))*ax,rr(ir,3)+dble(iz-nplotz(2)-1)*cx,&
!           engr(irpt(ir,ix,iy,imz(iz)),2)
!      write(117,"(4E20.10e2)") rr(ir,1)+dble(ix-nplotx(2))*ax,rr(ir,2)+dble(iy-nploty(2))*ax,rr(ir,3)+dble(iz-nplotz(2)-1)*cx,&
!           engr(irpt(ir,ix,iy,imz(iz)),3)
!
!      end if
!      ! gsite
!      if(magR(irpt(ir+2,ix,iy,imz(iz))).ne.0.0d0) then
!      write(109,"(4E20.10e2)") rr(ir+2,1)+dble(ix-nplotx(2))*ax,rr(ir+2,2)+dble(iy-nploty(2))*ax,rr(ir+2,3)+dble(iz-nplotz(2)-1)*cx,&
!           engr(irpt(ir+2,ix,iy,imz(iz)),4)
!      write(112,"(4E20.10e2)") rr(ir+2,1)+dble(ix-nplotx(2))*ax,rr(ir+2,2)+dble(iy-nploty(2))*ax,rr(ir+2,3)+dble(iz-nplotz(2)-1)*cx,&
!           engr(irpt(ir+2,ix,iy,imz(iz)),1)
!      write(115,"(4E20.10e2)") rr(ir+2,1)+dble(ix-nplotx(2))*ax,rr(ir+2,2)+dble(iy-nploty(2))*ax,rr(ir+2,3)+dble(iz-nplotz(2)-1)*cx,&
!           engr(irpt(ir+2,ix,iy,imz(iz)),2)
!      write(118,"(4E20.10e2)") rr(ir+2,1)+dble(ix-nplotx(2))*ax,rr(ir+2,2)+dble(iy-nploty(2))*ax,rr(ir+2,3)+dble(iz-nplotz(2)-1)*cx,&
!           engr(irpt(ir+2,ix,iy,imz(iz)),3)
!      end if
!   end do
!end do
!end do


!! surface
!! +y
!do ix=nplotx(1),nplotx(2)
!do iz=nplotz(1),nplotz(2)
!iy=nploty(2)
!   do ir=51,52 !+y
!      if(magFe(ifept(ir,ix,ipy(iy),iz)).ne.0.0d0) then
!      write(110,"(4E20.10e2)") rfe(ir,1)+dble(ix-nplotx(2))*ax,rfe(ir,2)+dble(iy+1-nploty(2))*ax,rfe(ir,3)+dble(iz-nplotz(2))*cx,&
!           engfe(ifept(ir,ix,ipy(iy),iz),4)
!      write(113,"(4E20.10e2)") rfe(ir,1)+dble(ix-nplotx(2))*ax,rfe(ir,2)+dble(iy+1-nploty(2))*ax,rfe(ir,3)+dble(iz-nplotz(2))*cx,&
!           engfe(ifept(ir,ix,ipy(iy),iz),1)
!      write(116,"(4E20.10e2)") rfe(ir,1)+dble(ix-nplotx(2))*ax,rfe(ir,2)+dble(iy+1-nploty(2))*ax,rfe(ir,3)+dble(iz-nplotz(2))*cx,&
!           engfe(ifept(ir,ix,ipy(iy),iz),2)
!      write(119,"(4E20.10e2)") rfe(ir,1)+dble(ix-nplotx(2))*ax,rfe(ir,2)+dble(iy+1-nploty(2))*ax,rfe(ir,3)+dble(iz-nplotz(2))*cx,&
!           engfe(ifept(ir,ix,ipy(iy),iz),3)
!      end if
!   end do
!   do ir=55,56 !+y
!      if(magFe(ifept(ir,ix,ipy(iy),iz)).ne.0.0d0) then
!      write(110,"(4E20.10e2)") rfe(ir,1)+dble(ix-nplotx(2))*ax,rfe(ir,2)+dble(iy+1-nploty(2))*ax,rfe(ir,3)+dble(iz-nplotz(2))*cx,&
!           engfe(ifept(ir,ix,ipy(iy),iz),4)
!      write(113,"(4E20.10e2)") rfe(ir,1)+dble(ix-nplotx(2))*ax,rfe(ir,2)+dble(iy+1-nploty(2))*ax,rfe(ir,3)+dble(iz-nplotz(2))*cx,&
!           engfe(ifept(ir,ix,ipy(iy),iz),1)
!      write(116,"(4E20.10e2)") rfe(ir,1)+dble(ix-nplotx(2))*ax,rfe(ir,2)+dble(iy+1-nploty(2))*ax,rfe(ir,3)+dble(iz-nplotz(2))*cx,&
!           engfe(ifept(ir,ix,ipy(iy),iz),2)
!      write(119,"(4E20.10e2)") rfe(ir,1)+dble(ix-nplotx(2))*ax,rfe(ir,2)+dble(iy+1-nploty(2))*ax,rfe(ir,3)+dble(iz-nplotz(2))*cx,&
!           engfe(ifept(ir,ix,ipy(iy),iz),3)
!      end if
!   end do
!end do
!end do
!
!!+x
!do iy=nploty(1),nploty(2)
!do iz=nplotz(1),nplotz(2)
!ix=nplotx(2)
!   do ir=51,54
!      if(magFe(ifept(ir,ipx(ix),iy,iz)).ne.0.0d0) then
!      write(110,"(4E20.10e2)") rfe(ir,1)+dble(ix+1-nplotx(2))*ax,rfe(ir,2)+dble(iy-nploty(2))*ax,rfe(ir,3)+dble(iz-nplotz(2))*cx,&
!           engfe(ifept(ir,ipx(ix),iy,iz),4)
!      write(113,"(4E20.10e2)") rfe(ir,1)+dble(ix+1-nplotx(2))*ax,rfe(ir,2)+dble(iy-nploty(2))*ax,rfe(ir,3)+dble(iz-nplotz(2))*cx,&
!           engfe(ifept(ir,ipx(ix),iy,iz),1)
!      write(116,"(4E20.10e2)") rfe(ir,1)+dble(ix+1-nplotx(2))*ax,rfe(ir,2)+dble(iy-nploty(2))*ax,rfe(ir,3)+dble(iz-nplotz(2))*cx,&
!           engfe(ifept(ir,ipx(ix),iy,iz),2)
!      write(119,"(4E20.10e2)") rfe(ir,1)+dble(ix+1-nplotx(2))*ax,rfe(ir,2)+dble(iy-nploty(2))*ax,rfe(ir,3)+dble(iz-nplotz(2))*cx,&
!           engfe(ifept(ir,ipx(ix),iy,iz),3)
!      end if
!   end do
!end do
!end do
!
!!-z
!do ix=nplotx(1),nplotx(2)
!do iy=nploty(1),nploty(2)
!iz=nplotz(1)
!   ir=53 !-z
!   if(magFe(ifept(ir,ix,iy,imz(iz))).ne.0.0d0) then
!   write(110,"(4E20.10e2)") rfe(ir,1)+dble(ix-nplotx(2))*ax,rfe(ir,2)+dble(iy-nploty(2))*ax,rfe(ir,3)+dble(iz-nplotz(2)-1)*cx,&
!        engfe(ifept(ir,ix,iy,imz(iz)),4)
!   write(113,"(4E20.10e2)") rfe(ir,1)+dble(ix-nplotx(2))*ax,rfe(ir,2)+dble(iy-nploty(2))*ax,rfe(ir,3)+dble(iz-nplotz(2)-1)*cx,&
!        engfe(ifept(ir,ix,iy,imz(iz)),1)
!   write(116,"(4E20.10e2)") rfe(ir,1)+dble(ix-nplotx(2))*ax,rfe(ir,2)+dble(iy-nploty(2))*ax,rfe(ir,3)+dble(iz-nplotz(2)-1)*cx,&
!        engfe(ifept(ir,ix,iy,imz(iz)),2)
!   write(119,"(4E20.10e2)") rfe(ir,1)+dble(ix-nplotx(2))*ax,rfe(ir,2)+dble(iy-nploty(2))*ax,rfe(ir,3)+dble(iz-nplotz(2)-1)*cx,&
!        engfe(ifept(ir,ix,iy,imz(iz)),3)
!   end if
!   ir=56 !-z
!   if(magFe(ifept(ir,ix,iy,imz(iz))).ne.0.0d0) then
!   write(110,"(4E20.10e2)") rfe(ir,1)+dble(ix-nplotx(2))*ax,rfe(ir,2)+dble(iy-nploty(2))*ax,rfe(ir,3)+dble(iz-nplotz(2)-1)*cx,&
!        engfe(ifept(ir,ix,iy,imz(iz)),4)
!   write(113,"(4E20.10e2)") rfe(ir,1)+dble(ix-nplotx(2))*ax,rfe(ir,2)+dble(iy-nploty(2))*ax,rfe(ir,3)+dble(iz-nplotz(2)-1)*cx,&
!        engfe(ifept(ir,ix,iy,imz(iz)),1)
!   write(116,"(4E20.10e2)") rfe(ir,1)+dble(ix-nplotx(2))*ax,rfe(ir,2)+dble(iy-nploty(2))*ax,rfe(ir,3)+dble(iz-nplotz(2)-1)*cx,&
!        engfe(ifept(ir,ix,iy,imz(iz)),2)
!   write(119,"(4E20.10e2)") rfe(ir,1)+dble(ix-nplotx(2))*ax,rfe(ir,2)+dble(iy-nploty(2))*ax,rfe(ir,3)+dble(iz-nplotz(2)-1)*cx,&
!        engfe(ifept(ir,ix,iy,imz(iz)),3)
!   end if
!end do
!end do
!
!! +y-z
!do ix=nplotx(1),nplotx(2)
!iy=nploty(2)
!iz=nplotz(1)
!   ir=56
!   if(magFe(ifept(ir,ix,ipy(iy),imz(iz))).ne.0.0d0) then
!   write(110,"(4E20.10e2)") rfe(ir,1)+dble(ix-nplotx(2))*ax,rfe(ir,2)+dble(iy+1-nploty(2))*ax,rfe(ir,3)+dble(iz-1-nplotz(2))*cx,&
!        engfe(ifept(ir,ix,ipy(iy),imz(iz)),4)
!   write(113,"(4E20.10e2)") rfe(ir,1)+dble(ix-nplotx(2))*ax,rfe(ir,2)+dble(iy+1-nploty(2))*ax,rfe(ir,3)+dble(iz-1-nplotz(2))*cx,&
!        engfe(ifept(ir,ix,ipy(iy),imz(iz)),1)
!   write(116,"(4E20.10e2)") rfe(ir,1)+dble(ix-nplotx(2))*ax,rfe(ir,2)+dble(iy+1-nploty(2))*ax,rfe(ir,3)+dble(iz-1-nplotz(2))*cx,&
!        engfe(ifept(ir,ix,ipy(iy),imz(iz)),2)
!   write(119,"(4E20.10e2)") rfe(ir,1)+dble(ix-nplotx(2))*ax,rfe(ir,2)+dble(iy+1-nploty(2))*ax,rfe(ir,3)+dble(iz-1-nplotz(2))*cx,&
!        engfe(ifept(ir,ix,ipy(iy),imz(iz)),3)
!   end if
!end  do
!
!! +x-z
!do iy=nploty(1),nploty(2)
!ix=nplotx(2)
!iz=nplotz(1)
!   ir=53
!   if(magFe(ifept(ir,ipx(ix),iy,imz(iz))).ne.0.0d0) then
!   write(110,"(4E20.10e2)") rfe(ir,1)+dble(ix+1-nplotx(2))*ax,rfe(ir,2)+dble(iy-nploty(2))*ax,rfe(ir,3)+dble(iz-1-nplotz(2))*cx,&
!        engfe(ifept(ir,ipx(ix),iy,imz(iz)),4)
!   write(113,"(4E20.10e2)") rfe(ir,1)+dble(ix+1-nplotx(2))*ax,rfe(ir,2)+dble(iy-nploty(2))*ax,rfe(ir,3)+dble(iz-1-nplotz(2))*cx,&
!        engfe(ifept(ir,ipx(ix),iy,imz(iz)),1)
!   write(116,"(4E20.10e2)") rfe(ir,1)+dble(ix+1-nplotx(2))*ax,rfe(ir,2)+dble(iy-nploty(2))*ax,rfe(ir,3)+dble(iz-1-nplotz(2))*cx,&
!        engfe(ifept(ir,ipx(ix),iy,imz(iz)),2)
!   write(119,"(4E20.10e2)") rfe(ir,1)+dble(ix+1-nplotx(2))*ax,rfe(ir,2)+dble(iy-nploty(2))*ax,rfe(ir,3)+dble(iz-1-nplotz(2))*cx,&
!        engfe(ifept(ir,ipx(ix),iy,imz(iz)),3)
!   end if
!end  do
!
!! +x+y
!do iz=nplotz(1),nplotz(2)
!ix=nplotx(2)
!iy=nploty(2)
!   do ir=51,52
!      if(magFe(ifept(ir,ipx(ix),ipy(iy),iz)).ne.0.0d0) then
!         write(110,"(4E20.10e2)") rfe(ir,1)+dble(ix+1-nplotx(2))*ax,rfe(ir,2)+dble(iy+1-nploty(2))*ax,rfe(ir,3)+dble(iz-nplotz(2))*cx,&
!              engfe(ifept(ir,ipx(ix),ipy(iy),iz),4)
!         write(113,"(4E20.10e2)") rfe(ir,1)+dble(ix+1-nplotx(2))*ax,rfe(ir,2)+dble(iy+1-nploty(2))*ax,rfe(ir,3)+dble(iz-nplotz(2))*cx,&
!              engfe(ifept(ir,ipx(ix),ipy(iy),iz),1)
!         write(116,"(4E20.10e2)") rfe(ir,1)+dble(ix+1-nplotx(2))*ax,rfe(ir,2)+dble(iy+1-nploty(2))*ax,rfe(ir,3)+dble(iz-nplotz(2))*cx,&
!              engfe(ifept(ir,ipx(ix),ipy(iy),iz),2)
!         write(119,"(4E20.10e2)") rfe(ir,1)+dble(ix+1-nplotx(2))*ax,rfe(ir,2)+dble(iy+1-nploty(2))*ax,rfe(ir,3)+dble(iz-nplotz(2))*cx,&
!              engfe(ifept(ir,ipx(ix),ipy(iy),iz),3)
!      end if
!   end do
!end  do


end subroutine plot_engdens
