subroutine plot_configxyz(itime,nunit,nr,nfe,amp,vr,vfe,ipx,imx,ipy,imy,ipz,imz,nunitx,nunity,nunitz,rr,rfe,irpt,ifept,ax,cx,nplotx,nploty,nplotz,igb)

implicit none
integer(4) :: nunit,nr,nfe,itime,iunit,i,nunitx,nunity,nunitz,ix,iy,iz,ir,irpt(nr,nunitx,nunity,nunitz),ifept(nfe,nunitx,nunity,nunitz)
integer(4) :: ipx(nunitx),ipy(nunity),ipz(nunitz),imx(nunitx),imy(nunity),imz(nunitz),nplotx(2),nploty(2),nplotz(2),igb(2)
real(8) :: amp(2),vr(nr*nunit,3),vfe(nfe*nunit,3),rr(nr,3),rfe(nfe,3),ax,cx,amp0
character(5) :: time

write(time,'(I5)') itime
open(8,file="movie_source/TDf"//trim(adjustl(time))//".txt") 
open(9,file="movie_source/TDg"//trim(adjustl(time))//".txt") 
open(10,file="movie_source/TDFe"//trim(adjustl(time))//".txt") 

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
      write(8,"(6f10.5)") rr(ir,1)+dble(ix-nplotx(2))*ax,rr(ir,2)+dble(iy-nploty(2))*ax,rr(ir,3)+dble(iz-nplotz(2))*cx,&
           amp0*vr(irpt(ir,ix,iy,iz),1),amp0*vr(irpt(ir,ix,iy,iz),2),amp0*vr(irpt(ir,ix,iy,iz),3) ! fsite, Nd1
      write(8,"(6f10.5)") rr(ir+4,1)+dble(ix-nplotx(2))*ax,rr(ir+4,2)+dble(iy-nploty(2))*ax,rr(ir+4,3)+dble(iz-nplotz(2))*cx,&
           amp0*vr(irpt(ir+4,ix,iy,iz),1),amp0*vr(irpt(ir+4,ix,iy,iz),2),amp0*vr(irpt(ir+4,ix,iy,iz),3) ! fsite, Nd1
      ! gsite
      write(9,"(6f10.5)") rr(ir+2,1)+dble(ix-nplotx(2))*ax,rr(ir+2,2)+dble(iy-nploty(2))*ax,rr(ir+2,3)+dble(iz-nplotz(2))*cx,&
           amp0*vr(irpt(ir+2,ix,iy,iz),1),amp0*vr(irpt(ir+2,ix,iy,iz),2),amp0*vr(irpt(ir+2,ix,iy,iz),3) ! gsite, Nd1
      write(9,"(6f10.5)") rr(ir+6,1)+dble(ix-nplotx(2))*ax,rr(ir+6,2)+dble(iy-nploty(2))*ax,rr(ir+6,3)+dble(iz-nplotz(2))*cx,&
           amp0*vr(irpt(ir+6,ix,iy,iz),1),amp0*vr(irpt(ir+6,ix,iy,iz),2),amp0*vr(irpt(ir+6,ix,iy,iz),3) ! gsite, Nd1
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
      write(8,"(6f10.5)") rr(ir,1)+dble(ix-nplotx(2))*ax,rr(ir,2)+dble(iy-nploty(2))*ax,rr(ir,3)+dble(iz-nplotz(2)-1)*cx,&
           amp0*vr(irpt(ir,ix,iy,imz(iz)),1),amp0*vr(irpt(ir,ix,iy,imz(iz)),2),amp0*vr(irpt(ir,ix,iy,imz(iz)),3) ! fsite, Nd1
      ! gsite
      write(9,"(6f10.5)") rr(ir+2,1)+dble(ix-nplotx(2))*ax,rr(ir+2,2)+dble(iy-nploty(2))*ax,rr(ir+2,3)+dble(iz-nplotz(2)-1)*cx,&
           amp0*vr(irpt(ir+2,ix,iy,imz(iz)),1),amp0*vr(irpt(ir+2,ix,iy,imz(iz)),2),amp0*vr(irpt(ir+2,ix,iy,imz(iz)),3) ! gsite, Nd1
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
      write(10,"(6f10.5)") rfe(ir,1)+dble(ix-nplotx(2))*ax,rfe(ir,2)+dble(iy-nploty(2))*ax,rfe(ir,3)+dble(iz-nplotz(2))*cx,&
           amp0*vfe(ifept(ir,ix,iy,iz),1),amp0*vfe(ifept(ir,ix,iy,iz),2),amp0*vfe(ifept(ir,ix,iy,iz),3) ! Fesite
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

   ir=51 !+y
   write(10,"(6f10.5)") rfe(ir,1)+dble(ix-nplotx(2))*ax,rfe(ir,2)+dble(iy+1-nploty(2))*ax,rfe(ir,3)+dble(iz-nplotz(2))*cx,&
        amp0*vfe(ifept(ir,ix,ipy(iy),iz),1),amp0*vfe(ifept(ir,ix,ipy(iy),iz),2),amp0*vfe(ifept(ir,ix,ipy(iy),iz),3) ! Fesite
   ir=52 !+y
   write(10,"(6f10.5)") rfe(ir,1)+dble(ix-nplotx(2))*ax,rfe(ir,2)+dble(iy+1-nploty(2))*ax,rfe(ir,3)+dble(iz-nplotz(2))*cx,&
        amp0*vfe(ifept(ir,ix,ipy(iy),iz),1),amp0*vfe(ifept(ir,ix,ipy(iy),iz),2),amp0*vfe(ifept(ir,ix,ipy(iy),iz),3) ! Fesite
   ir=55 !+y
   write(10,"(6f10.5)") rfe(ir,1)+dble(ix-nplotx(2))*ax,rfe(ir,2)+dble(iy+1-nploty(2))*ax,rfe(ir,3)+dble(iz-nplotz(2))*cx,&
        amp0*vfe(ifept(ir,ix,ipy(iy),iz),1),amp0*vfe(ifept(ir,ix,ipy(iy),iz),2),amp0*vfe(ifept(ir,ix,ipy(iy),iz),3) ! Fesite
   ir=56 !+y
   write(10,"(6f10.5)") rfe(ir,1)+dble(ix-nplotx(2))*ax,rfe(ir,2)+dble(iy+1-nploty(2))*ax,rfe(ir,3)+dble(iz-nplotz(2))*cx,&
        amp0*vfe(ifept(ir,ix,ipy(iy),iz),1),amp0*vfe(ifept(ir,ix,ipy(iy),iz),2),amp0*vfe(ifept(ir,ix,ipy(iy),iz),3) ! Fesite
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

   ir=51 
   write(10,"(6f10.5)") rfe(ir,1)+dble(ix+1-nplotx(2))*ax,rfe(ir,2)+dble(iy-nploty(2))*ax,rfe(ir,3)+dble(iz-nplotz(2))*cx,&
        amp0*vfe(ifept(ir,ipx(ix),iy,iz),1),amp0*vfe(ifept(ir,ipx(ix),iy,iz),2),amp0*vfe(ifept(ir,ipx(ix),iy,iz),3) ! Fesite
   ir=52 
   write(10,"(6f10.5)") rfe(ir,1)+dble(ix+1-nplotx(2))*ax,rfe(ir,2)+dble(iy-nploty(2))*ax,rfe(ir,3)+dble(iz-nplotz(2))*cx,&
        amp0*vfe(ifept(ir,ipx(ix),iy,iz),1),amp0*vfe(ifept(ir,ipx(ix),iy,iz),2),amp0*vfe(ifept(ir,ipx(ix),iy,iz),3) ! Fesite
   ir=53 
   write(10,"(6f10.5)") rfe(ir,1)+dble(ix+1-nplotx(2))*ax,rfe(ir,2)+dble(iy-nploty(2))*ax,rfe(ir,3)+dble(iz-nplotz(2))*cx,&
        amp0*vfe(ifept(ir,ipx(ix),iy,iz),1),amp0*vfe(ifept(ir,ipx(ix),iy,iz),2),amp0*vfe(ifept(ir,ipx(ix),iy,iz),3) ! Fesite
   ir=54 
   write(10,"(6f10.5)") rfe(ir,1)+dble(ix+1-nplotx(2))*ax,rfe(ir,2)+dble(iy-nploty(2))*ax,rfe(ir,3)+dble(iz-nplotz(2))*cx,&
        amp0*vfe(ifept(ir,ipx(ix),iy,iz),1),amp0*vfe(ifept(ir,ipx(ix),iy,iz),2),amp0*vfe(ifept(ir,ipx(ix),iy,iz),3) ! Fesite
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
   write(10,"(6f10.5)") rfe(ir,1)+dble(ix-nplotx(2))*ax,rfe(ir,2)+dble(iy-nploty(2))*ax,rfe(ir,3)+dble(iz-nplotz(2)-1)*cx,&
        amp0*vfe(ifept(ir,ix,iy,imz(iz)),1),amp0*vfe(ifept(ir,ix,iy,imz(iz)),2),amp0*vfe(ifept(ir,ix,iy,imz(iz)),3) ! Fesite
   ir=56 !-z
   write(10,"(6f10.5)") rfe(ir,1)+dble(ix-nplotx(2))*ax,rfe(ir,2)+dble(iy-nploty(2))*ax,rfe(ir,3)+dble(iz-nplotz(2)-1)*cx,&
        amp0*vfe(ifept(ir,ix,iy,imz(iz)),1),amp0*vfe(ifept(ir,ix,iy,imz(iz)),2),amp0*vfe(ifept(ir,ix,iy,imz(iz)),3) ! Fesite
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
   write(10,"(6f10.5)") rfe(ir,1)+dble(ix-nplotx(2))*ax,rfe(ir,2)+dble(iy+1-nploty(2))*ax,rfe(ir,3)+dble(iz-1-nplotz(2))*cx,&
        amp0*vfe(ifept(ir,ix,ipy(iy),imz(iz)),1),amp0*vfe(ifept(ir,ix,ipy(iy),imz(iz)),2),amp0*vfe(ifept(ir,ix,ipy(iy),imz(iz)),3) ! Fesite
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
   write(10,"(6f10.5)") rfe(ir,1)+dble(ix+1-nplotx(2))*ax,rfe(ir,2)+dble(iy-nploty(2))*ax,rfe(ir,3)+dble(iz-1-nplotz(2))*cx,&
        amp0*vfe(ifept(ir,ipx(ix),iy,imz(iz)),1),amp0*vfe(ifept(ir,ipx(ix),iy,imz(iz)),2),amp0*vfe(ifept(ir,ipx(ix),iy,imz(iz)),3) ! Fesite
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

   ir=51
   write(10,"(6f10.5)") rfe(ir,1)+dble(ix+1-nplotx(2))*ax,rfe(ir,2)+dble(iy+1-nploty(2))*ax,rfe(ir,3)+dble(iz-nplotz(2))*cx,&
        amp0*vfe(ifept(ir,ipx(ix),ipy(iy),iz),1),amp0*vfe(ifept(ir,ipx(ix),ipy(iy),iz),2),amp0*vfe(ifept(ir,ipx(ix),ipy(iy),iz),3) ! Fesite
   ir=52
   write(10,"(6f10.5)") rfe(ir,1)+dble(ix+1-nplotx(2))*ax,rfe(ir,2)+dble(iy+1-nploty(2))*ax,rfe(ir,3)+dble(iz-nplotz(2))*cx,&
        amp0*vfe(ifept(ir,ipx(ix),ipy(iy),iz),1),amp0*vfe(ifept(ir,ipx(ix),ipy(iy),iz),2),amp0*vfe(ifept(ir,ipx(ix),ipy(iy),iz),3) ! Fesite

end  do


end subroutine plot_configxyz
