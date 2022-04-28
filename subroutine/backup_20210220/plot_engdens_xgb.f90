subroutine plot_engdens(itime,nunit,nr,nfe,amp,engr,engfe,ipx,imx,ipy,imy,ipz,imz,nunitx,nunity,nunitz,rr,rfe,irpt,ifept,&
     irpt0,ifept0,ax,cx,nplotx,nploty,nplotz,igb)
implicit none
integer(4) :: nunit,nr,nfe,itime,iunit,i,nunitx,nunity,nunitz,ix,iy,iz,ir,irpt(nr,nunitx,nunity,nunitz),ifept(nfe,nunitx,nunity,nunitz)
integer(4) :: irpt0(nr,nunitx,nunity,nunitz),ifept0(nfe,nunitx,nunity,nunitz)
integer(4) :: ipx(nunitx),ipy(nunity),ipz(nunitz),imx(nunitx),imy(nunity),imz(nunitz),nplotx(2),nploty(2),nplotz(2),igb(2)
real(8) :: amp(2),engr(nr*nunit,4),engfe(nfe*nunit,4),rr(nr,3),rfe(nfe,3),ax,cx,amp0
character(5) :: time

write(time,'(I5)') itime
open(108,file="../movie_source/engtotTDf"//trim(adjustl(time))//".txt") 
open(109,file="../movie_source/engtotTDg"//trim(adjustl(time))//".txt") 
open(110,file="../movie_source/engtotTDFe"//trim(adjustl(time))//".txt") 

open(111,file="../movie_source/engHATDf"//trim(adjustl(time))//".txt") 
open(112,file="../movie_source/engHATDg"//trim(adjustl(time))//".txt") 
open(113,file="../movie_source/engHATDFe"//trim(adjustl(time))//".txt") 

open(114,file="../movie_source/engHkTDf"//trim(adjustl(time))//".txt") 
open(115,file="../movie_source/engHkTDg"//trim(adjustl(time))//".txt") 
open(116,file="../movie_source/engHkTDFe"//trim(adjustl(time))//".txt") 

open(117,file="../movie_source/engHextTDf"//trim(adjustl(time))//".txt") 
open(118,file="../movie_source/engHextTDg"//trim(adjustl(time))//".txt") 
open(119,file="../movie_source/engHextTDFe"//trim(adjustl(time))//".txt") 

do iz=nplotz(1),nplotz(2)
do iy=nploty(1),nploty(2)
do ix=nplotx(1),nplotx(2)
!do ix=nunitx-nplotx,nunitx-1
   if(ix.ge.igb(1).and.ix.le.igb(2)) then
      amp0=amp(2)
   else
      amp0=amp(1)
   end if
   do ir=1,2
      ! inside
      ! fsite
      if(irpt0(ir,ix,iy,iz).eq.1) then
      write(108,"(3f10.5,f15.5)") rr(ir,1)+dble(ix-nplotx(2))*ax,rr(ir,2)+dble(iy-nploty(2))*ax,rr(ir,3)+dble(iz-nplotz(2))*cx,&
           amp0*engr(irpt(ir,ix,iy,iz),4)
      write(111,"(3f10.5,f15.5)") rr(ir,1)+dble(ix-nplotx(2))*ax,rr(ir,2)+dble(iy-nploty(2))*ax,rr(ir,3)+dble(iz-nplotz(2))*cx,&
           amp0*engr(irpt(ir,ix,iy,iz),1)
      write(114,"(3f10.5,f15.5)") rr(ir,1)+dble(ix-nplotx(2))*ax,rr(ir,2)+dble(iy-nploty(2))*ax,rr(ir,3)+dble(iz-nplotz(2))*cx,&
           amp0*engr(irpt(ir,ix,iy,iz),2)
      write(117,"(3f10.5,f15.5)") rr(ir,1)+dble(ix-nplotx(2))*ax,rr(ir,2)+dble(iy-nploty(2))*ax,rr(ir,3)+dble(iz-nplotz(2))*cx,&
           amp0*engr(irpt(ir,ix,iy,iz),3)
      end if
      if(irpt0(ir+4,ix,iy,iz).eq.1) then
      write(108,"(3f10.5,f15.5)") rr(ir+4,1)+dble(ix-nplotx(2))*ax,rr(ir+4,2)+dble(iy-nploty(2))*ax,rr(ir+4,3)+dble(iz-nplotz(2))*cx,&
           amp0*engr(irpt(ir+4,ix,iy,iz),4)
      write(111,"(3f10.5,f15.5)") rr(ir+4,1)+dble(ix-nplotx(2))*ax,rr(ir+4,2)+dble(iy-nploty(2))*ax,rr(ir+4,3)+dble(iz-nplotz(2))*cx,&
           amp0*engr(irpt(ir+4,ix,iy,iz),1)
      write(114,"(3f10.5,f15.5)") rr(ir+4,1)+dble(ix-nplotx(2))*ax,rr(ir+4,2)+dble(iy-nploty(2))*ax,rr(ir+4,3)+dble(iz-nplotz(2))*cx,&
           amp0*engr(irpt(ir+4,ix,iy,iz),2)
      write(117,"(3f10.5,f15.5)") rr(ir+4,1)+dble(ix-nplotx(2))*ax,rr(ir+4,2)+dble(iy-nploty(2))*ax,rr(ir+4,3)+dble(iz-nplotz(2))*cx,&
           amp0*engr(irpt(ir+4,ix,iy,iz),3)
      end if
      ! gsite
      if(irpt0(ir+2,ix,iy,iz).eq.1) then
      write(109,"(3f10.5,f15.5)") rr(ir+2,1)+dble(ix-nplotx(2))*ax,rr(ir+2,2)+dble(iy-nploty(2))*ax,rr(ir+2,3)+dble(iz-nplotz(2))*cx,&
           amp0*engr(irpt(ir+2,ix,iy,iz),4)
      write(112,"(3f10.5,f15.5)") rr(ir+2,1)+dble(ix-nplotx(2))*ax,rr(ir+2,2)+dble(iy-nploty(2))*ax,rr(ir+2,3)+dble(iz-nplotz(2))*cx,&
           amp0*engr(irpt(ir+2,ix,iy,iz),1)
      write(115,"(3f10.5,f15.5)") rr(ir+2,1)+dble(ix-nplotx(2))*ax,rr(ir+2,2)+dble(iy-nploty(2))*ax,rr(ir+2,3)+dble(iz-nplotz(2))*cx,&
           amp0*engr(irpt(ir+2,ix,iy,iz),2)
      write(118,"(3f10.5,f15.5)") rr(ir+2,1)+dble(ix-nplotx(2))*ax,rr(ir+2,2)+dble(iy-nploty(2))*ax,rr(ir+2,3)+dble(iz-nplotz(2))*cx,&
           amp0*engr(irpt(ir+2,ix,iy,iz),3)
      end if
      if(irpt0(ir+6,ix,iy,iz).eq.1) then
      write(109,"(3f10.5,f15.5)") rr(ir+6,1)+dble(ix-nplotx(2))*ax,rr(ir+6,2)+dble(iy-nploty(2))*ax,rr(ir+6,3)+dble(iz-nplotz(2))*cx,&
           amp0*engr(irpt(ir+6,ix,iy,iz),4)
      write(112,"(3f10.5,f15.5)") rr(ir+6,1)+dble(ix-nplotx(2))*ax,rr(ir+6,2)+dble(iy-nploty(2))*ax,rr(ir+6,3)+dble(iz-nplotz(2))*cx,&
           amp0*engr(irpt(ir+6,ix,iy,iz),1)
      write(115,"(3f10.5,f15.5)") rr(ir+6,1)+dble(ix-nplotx(2))*ax,rr(ir+6,2)+dble(iy-nploty(2))*ax,rr(ir+6,3)+dble(iz-nplotz(2))*cx,&
           amp0*engr(irpt(ir+6,ix,iy,iz),2)
      write(118,"(3f10.5,f15.5)") rr(ir+6,1)+dble(ix-nplotx(2))*ax,rr(ir+6,2)+dble(iy-nploty(2))*ax,rr(ir+6,3)+dble(iz-nplotz(2))*cx,&
           amp0*engr(irpt(ir+6,ix,iy,iz),3)
      end if
   end do
end do
end do
end do
! bottom
do iy=nploty(1),nploty(2)
do ix=nplotx(1),nplotx(2)
iz=nplotz(1)
!do ix=nunitx-nplotx,nunitx-1
   if(ix.ge.igb(1).and.ix.le.igb(2)) then
      amp0=amp(2)
   else
      amp0=amp(1)
   end if

   do ir=1,2
      ! fsite
      if(irpt0(ir,ix,iy,imz(iz)).eq.1) then
      write(108,"(3f10.5,f15.5)") rr(ir,1)+dble(ix-nplotx(2))*ax,rr(ir,2)+dble(iy-nploty(2))*ax,rr(ir,3)+dble(iz-nplotz(2)-1)*cx,&
           amp0*engr(irpt(ir,ix,iy,imz(iz)),4)
      write(111,"(3f10.5,f15.5)") rr(ir,1)+dble(ix-nplotx(2))*ax,rr(ir,2)+dble(iy-nploty(2))*ax,rr(ir,3)+dble(iz-nplotz(2)-1)*cx,&
           amp0*engr(irpt(ir,ix,iy,imz(iz)),1)
      write(114,"(3f10.5,f15.5)") rr(ir,1)+dble(ix-nplotx(2))*ax,rr(ir,2)+dble(iy-nploty(2))*ax,rr(ir,3)+dble(iz-nplotz(2)-1)*cx,&
           amp0*engr(irpt(ir,ix,iy,imz(iz)),2)
      write(117,"(3f10.5,f15.5)") rr(ir,1)+dble(ix-nplotx(2))*ax,rr(ir,2)+dble(iy-nploty(2))*ax,rr(ir,3)+dble(iz-nplotz(2)-1)*cx,&
           amp0*engr(irpt(ir,ix,iy,imz(iz)),3)

      end if
      ! gsite
      if(irpt0(ir+2,ix,iy,imz(iz)).eq.1) then
      write(109,"(3f10.5,f15.5)") rr(ir+2,1)+dble(ix-nplotx(2))*ax,rr(ir+2,2)+dble(iy-nploty(2))*ax,rr(ir+2,3)+dble(iz-nplotz(2)-1)*cx,&
           amp0*engr(irpt(ir+2,ix,iy,imz(iz)),4)
      write(112,"(3f10.5,f15.5)") rr(ir+2,1)+dble(ix-nplotx(2))*ax,rr(ir+2,2)+dble(iy-nploty(2))*ax,rr(ir+2,3)+dble(iz-nplotz(2)-1)*cx,&
           amp0*engr(irpt(ir+2,ix,iy,imz(iz)),1)
      write(115,"(3f10.5,f15.5)") rr(ir+2,1)+dble(ix-nplotx(2))*ax,rr(ir+2,2)+dble(iy-nploty(2))*ax,rr(ir+2,3)+dble(iz-nplotz(2)-1)*cx,&
           amp0*engr(irpt(ir+2,ix,iy,imz(iz)),2)
      write(118,"(3f10.5,f15.5)") rr(ir+2,1)+dble(ix-nplotx(2))*ax,rr(ir+2,2)+dble(iy-nploty(2))*ax,rr(ir+2,3)+dble(iz-nplotz(2)-1)*cx,&
           amp0*engr(irpt(ir+2,ix,iy,imz(iz)),3)
      end if
   end do
end do
end do

do iz=nplotz(1),nplotz(2)
do iy=nploty(1),nploty(2)
do ix=nplotx(1),nplotx(2)
!do ix=nunitx-nplotx,nunitx-1
   if(ix.ge.igb(1).and.ix.le.igb(2)) then
      amp0=amp(2)
   else
      amp0=amp(1)
   end if

   ! inside
   do ir=1,nfe
      if(ifept0(ir,ix,iy,iz).eq.1) then
      write(110,"(3f10.5,f15.5)") rfe(ir,1)+dble(ix-nplotx(2))*ax,rfe(ir,2)+dble(iy-nploty(2))*ax,rfe(ir,3)+dble(iz-nplotz(2))*cx,&
           amp0*engfe(ifept(ir,ix,iy,iz),4)
      write(113,"(3f10.5,f15.5)") rfe(ir,1)+dble(ix-nplotx(2))*ax,rfe(ir,2)+dble(iy-nploty(2))*ax,rfe(ir,3)+dble(iz-nplotz(2))*cx,&
           amp0*engfe(ifept(ir,ix,iy,iz),1)
      write(116,"(3f10.5,f15.5)") rfe(ir,1)+dble(ix-nplotx(2))*ax,rfe(ir,2)+dble(iy-nploty(2))*ax,rfe(ir,3)+dble(iz-nplotz(2))*cx,&
           amp0*engfe(ifept(ir,ix,iy,iz),2)
      write(119,"(3f10.5,f15.5)") rfe(ir,1)+dble(ix-nplotx(2))*ax,rfe(ir,2)+dble(iy-nploty(2))*ax,rfe(ir,3)+dble(iz-nplotz(2))*cx,&
           amp0*engfe(ifept(ir,ix,iy,iz),3)
      end if
   end do
end do
end do
end do

! surface
! +y
do ix=nplotx(1),nplotx(2)
do iz=nplotz(1),nplotz(2)
iy=nploty(2)
   if(ix.ge.igb(1).and.ix.le.igb(2)) then
      amp0=amp(2)
   else
      amp0=amp(1)
   end if

   do ir=51,52 !+y
      if(ifept0(ir,ix,ipy(iy),iz).eq.1) then
      write(110,"(3f10.5,f15.5)") rfe(ir,1)+dble(ix-nplotx(2))*ax,rfe(ir,2)+dble(iy+1-nploty(2))*ax,rfe(ir,3)+dble(iz-nplotz(2))*cx,&
           amp0*engfe(ifept(ir,ix,ipy(iy),iz),4)
      write(113,"(3f10.5,f15.5)") rfe(ir,1)+dble(ix-nplotx(2))*ax,rfe(ir,2)+dble(iy+1-nploty(2))*ax,rfe(ir,3)+dble(iz-nplotz(2))*cx,&
           amp0*engfe(ifept(ir,ix,ipy(iy),iz),1)
      write(116,"(3f10.5,f15.5)") rfe(ir,1)+dble(ix-nplotx(2))*ax,rfe(ir,2)+dble(iy+1-nploty(2))*ax,rfe(ir,3)+dble(iz-nplotz(2))*cx,&
           amp0*engfe(ifept(ir,ix,ipy(iy),iz),2)
      write(119,"(3f10.5,f15.5)") rfe(ir,1)+dble(ix-nplotx(2))*ax,rfe(ir,2)+dble(iy+1-nploty(2))*ax,rfe(ir,3)+dble(iz-nplotz(2))*cx,&
           amp0*engfe(ifept(ir,ix,ipy(iy),iz),3)
      end if
   end do
   do ir=55,56 !+y
      if(ifept0(ir,ix,ipy(iy),iz).eq.1) then
      write(110,"(3f10.5,f15.5)") rfe(ir,1)+dble(ix-nplotx(2))*ax,rfe(ir,2)+dble(iy+1-nploty(2))*ax,rfe(ir,3)+dble(iz-nplotz(2))*cx,&
           amp0*engfe(ifept(ir,ix,ipy(iy),iz),4)
      write(113,"(3f10.5,f15.5)") rfe(ir,1)+dble(ix-nplotx(2))*ax,rfe(ir,2)+dble(iy+1-nploty(2))*ax,rfe(ir,3)+dble(iz-nplotz(2))*cx,&
           amp0*engfe(ifept(ir,ix,ipy(iy),iz),1)
      write(116,"(3f10.5,f15.5)") rfe(ir,1)+dble(ix-nplotx(2))*ax,rfe(ir,2)+dble(iy+1-nploty(2))*ax,rfe(ir,3)+dble(iz-nplotz(2))*cx,&
           amp0*engfe(ifept(ir,ix,ipy(iy),iz),2)
      write(119,"(3f10.5,f15.5)") rfe(ir,1)+dble(ix-nplotx(2))*ax,rfe(ir,2)+dble(iy+1-nploty(2))*ax,rfe(ir,3)+dble(iz-nplotz(2))*cx,&
           amp0*engfe(ifept(ir,ix,ipy(iy),iz),3)
      end if
   end do
end do
end do

!+x
do iy=nploty(1),nploty(2)
do iz=nplotz(1),nplotz(2)
ix=nplotx(2)
   if(ix.ge.igb(1).and.ix.le.igb(2)) then
      amp0=amp(2)
   else
      amp0=amp(1)
   end if

   do ir=51,54
      if(ifept0(ir,ipx(ix),iy,iz).eq.1) then
      write(110,"(3f10.5,f15.5)") rfe(ir,1)+dble(ix+1-nplotx(2))*ax,rfe(ir,2)+dble(iy-nploty(2))*ax,rfe(ir,3)+dble(iz-nplotz(2))*cx,&
           amp0*engfe(ifept(ir,ipx(ix),iy,iz),4)
      write(113,"(3f10.5,f15.5)") rfe(ir,1)+dble(ix+1-nplotx(2))*ax,rfe(ir,2)+dble(iy-nploty(2))*ax,rfe(ir,3)+dble(iz-nplotz(2))*cx,&
           amp0*engfe(ifept(ir,ipx(ix),iy,iz),1)
      write(116,"(3f10.5,f15.5)") rfe(ir,1)+dble(ix+1-nplotx(2))*ax,rfe(ir,2)+dble(iy-nploty(2))*ax,rfe(ir,3)+dble(iz-nplotz(2))*cx,&
           amp0*engfe(ifept(ir,ipx(ix),iy,iz),2)
      write(119,"(3f10.5,f15.5)") rfe(ir,1)+dble(ix+1-nplotx(2))*ax,rfe(ir,2)+dble(iy-nploty(2))*ax,rfe(ir,3)+dble(iz-nplotz(2))*cx,&
           amp0*engfe(ifept(ir,ipx(ix),iy,iz),3)
      end if
   end do
end do
end do

!-z
do ix=nplotx(1),nplotx(2)
do iy=nploty(1),nploty(2)
iz=nplotz(1)
   if(ix.ge.igb(1).and.ix.le.igb(2)) then
      amp0=amp(2)
   else
      amp0=amp(1)
   end if

   ir=53 !-z
   if(ifept0(ir,ix,iy,imz(iz)).eq.1) then
   write(110,"(3f10.5,f15.5)") rfe(ir,1)+dble(ix-nplotx(2))*ax,rfe(ir,2)+dble(iy-nploty(2))*ax,rfe(ir,3)+dble(iz-nplotz(2)-1)*cx,&
        amp0*engfe(ifept(ir,ix,iy,imz(iz)),4)
   write(113,"(3f10.5,f15.5)") rfe(ir,1)+dble(ix-nplotx(2))*ax,rfe(ir,2)+dble(iy-nploty(2))*ax,rfe(ir,3)+dble(iz-nplotz(2)-1)*cx,&
        amp0*engfe(ifept(ir,ix,iy,imz(iz)),1)
   write(116,"(3f10.5,f15.5)") rfe(ir,1)+dble(ix-nplotx(2))*ax,rfe(ir,2)+dble(iy-nploty(2))*ax,rfe(ir,3)+dble(iz-nplotz(2)-1)*cx,&
        amp0*engfe(ifept(ir,ix,iy,imz(iz)),2)
   write(119,"(3f10.5,f15.5)") rfe(ir,1)+dble(ix-nplotx(2))*ax,rfe(ir,2)+dble(iy-nploty(2))*ax,rfe(ir,3)+dble(iz-nplotz(2)-1)*cx,&
        amp0*engfe(ifept(ir,ix,iy,imz(iz)),3)
   end if
   ir=56 !-z
   if(ifept0(ir,ix,iy,imz(iz)).eq.1) then
   write(110,"(3f10.5,f15.5)") rfe(ir,1)+dble(ix-nplotx(2))*ax,rfe(ir,2)+dble(iy-nploty(2))*ax,rfe(ir,3)+dble(iz-nplotz(2)-1)*cx,&
        amp0*engfe(ifept(ir,ix,iy,imz(iz)),4)
   write(113,"(3f10.5,f15.5)") rfe(ir,1)+dble(ix-nplotx(2))*ax,rfe(ir,2)+dble(iy-nploty(2))*ax,rfe(ir,3)+dble(iz-nplotz(2)-1)*cx,&
        amp0*engfe(ifept(ir,ix,iy,imz(iz)),1)
   write(116,"(3f10.5,f15.5)") rfe(ir,1)+dble(ix-nplotx(2))*ax,rfe(ir,2)+dble(iy-nploty(2))*ax,rfe(ir,3)+dble(iz-nplotz(2)-1)*cx,&
        amp0*engfe(ifept(ir,ix,iy,imz(iz)),2)
   write(119,"(3f10.5,f15.5)") rfe(ir,1)+dble(ix-nplotx(2))*ax,rfe(ir,2)+dble(iy-nploty(2))*ax,rfe(ir,3)+dble(iz-nplotz(2)-1)*cx,&
        amp0*engfe(ifept(ir,ix,iy,imz(iz)),3)
   end if
end do
end do

! +y-z
do ix=nplotx(1),nplotx(2)
iy=nploty(2)
iz=nplotz(1)
   if(ix.ge.igb(1).and.ix.le.igb(2)) then
      amp0=amp(2)
   else
      amp0=amp(1)
   end if

   ir=56
   if(ifept0(ir,ix,ipy(iy),imz(iz)).eq.1) then
   write(110,"(3f10.5,f15.5)") rfe(ir,1)+dble(ix-nplotx(2))*ax,rfe(ir,2)+dble(iy+1-nploty(2))*ax,rfe(ir,3)+dble(iz-1-nplotz(2))*cx,&
        amp0*engfe(ifept(ir,ix,ipy(iy),imz(iz)),4)
   write(113,"(3f10.5,f15.5)") rfe(ir,1)+dble(ix-nplotx(2))*ax,rfe(ir,2)+dble(iy+1-nploty(2))*ax,rfe(ir,3)+dble(iz-1-nplotz(2))*cx,&
        amp0*engfe(ifept(ir,ix,ipy(iy),imz(iz)),1)
   write(116,"(3f10.5,f15.5)") rfe(ir,1)+dble(ix-nplotx(2))*ax,rfe(ir,2)+dble(iy+1-nploty(2))*ax,rfe(ir,3)+dble(iz-1-nplotz(2))*cx,&
        amp0*engfe(ifept(ir,ix,ipy(iy),imz(iz)),2)
   write(119,"(3f10.5,f15.5)") rfe(ir,1)+dble(ix-nplotx(2))*ax,rfe(ir,2)+dble(iy+1-nploty(2))*ax,rfe(ir,3)+dble(iz-1-nplotz(2))*cx,&
        amp0*engfe(ifept(ir,ix,ipy(iy),imz(iz)),3)
   end if
end  do

! +x-z
do iy=nploty(1),nploty(2)
ix=nplotx(2)
iz=nplotz(1)
   if(ix.ge.igb(1).and.ix.le.igb(2)) then
      amp0=amp(2)
   else
      amp0=amp(1)
   end if

   ir=53
   if(ifept0(ir,ipx(ix),iy,imz(iz)).eq.1) then
   write(110,"(3f10.5,f15.5)") rfe(ir,1)+dble(ix+1-nplotx(2))*ax,rfe(ir,2)+dble(iy-nploty(2))*ax,rfe(ir,3)+dble(iz-1-nplotz(2))*cx,&
        amp0*engfe(ifept(ir,ipx(ix),iy,imz(iz)),4)
   write(113,"(3f10.5,f15.5)") rfe(ir,1)+dble(ix+1-nplotx(2))*ax,rfe(ir,2)+dble(iy-nploty(2))*ax,rfe(ir,3)+dble(iz-1-nplotz(2))*cx,&
        amp0*engfe(ifept(ir,ipx(ix),iy,imz(iz)),1)
   write(116,"(3f10.5,f15.5)") rfe(ir,1)+dble(ix+1-nplotx(2))*ax,rfe(ir,2)+dble(iy-nploty(2))*ax,rfe(ir,3)+dble(iz-1-nplotz(2))*cx,&
        amp0*engfe(ifept(ir,ipx(ix),iy,imz(iz)),2)
   write(119,"(3f10.5,f15.5)") rfe(ir,1)+dble(ix+1-nplotx(2))*ax,rfe(ir,2)+dble(iy-nploty(2))*ax,rfe(ir,3)+dble(iz-1-nplotz(2))*cx,&
        amp0*engfe(ifept(ir,ipx(ix),iy,imz(iz)),3)
   end if
end  do

! +x+y
do iz=nplotz(1),nplotz(2)
ix=nplotx(2)
iy=nploty(2)
   if(ix.ge.igb(1).and.ix.le.igb(2)) then
      amp0=amp(2)
   else
      amp0=amp(1)
   end if

   do ir=51,52
      if(ifept0(ir,ipx(ix),ipy(iy),iz).eq.1) then
         write(110,"(3f10.5,f15.5)") rfe(ir,1)+dble(ix+1-nplotx(2))*ax,rfe(ir,2)+dble(iy+1-nploty(2))*ax,rfe(ir,3)+dble(iz-nplotz(2))*cx,&
              amp0*engfe(ifept(ir,ipx(ix),ipy(iy),iz),4)
         write(113,"(3f10.5,f15.5)") rfe(ir,1)+dble(ix+1-nplotx(2))*ax,rfe(ir,2)+dble(iy+1-nploty(2))*ax,rfe(ir,3)+dble(iz-nplotz(2))*cx,&
              amp0*engfe(ifept(ir,ipx(ix),ipy(iy),iz),1)
         write(116,"(3f10.5,f15.5)") rfe(ir,1)+dble(ix+1-nplotx(2))*ax,rfe(ir,2)+dble(iy+1-nploty(2))*ax,rfe(ir,3)+dble(iz-nplotz(2))*cx,&
              amp0*engfe(ifept(ir,ipx(ix),ipy(iy),iz),2)
         write(119,"(3f10.5,f15.5)") rfe(ir,1)+dble(ix+1-nplotx(2))*ax,rfe(ir,2)+dble(iy+1-nploty(2))*ax,rfe(ir,3)+dble(iz-nplotz(2))*cx,&
              amp0*engfe(ifept(ir,ipx(ix),ipy(iy),iz),3)
      end if
   end do
end  do


end subroutine plot_engdens
