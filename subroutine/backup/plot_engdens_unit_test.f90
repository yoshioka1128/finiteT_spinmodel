subroutine plot_engdens_unit(itime,nunit,nunitxbd,nunitybd,nunitzbd,engrunit,engfeunit,magunit,nunitx,nunity,nunitz,iunit,ax,cx)
implicit none
integer(4) :: nunit,itime,nunitx,nunity,nunitz,ix,iy,iz,ir,iunit(nunitx,nunity,nunitz)
integer(4) :: nunitxbd(2),nunitybd(2),nunitzbd(2),i
real(8) :: engrunit(nunit,4),engfeunit(nunit,4),ax,cx,magunit(nunit,3)
character(5) :: time

write(time,'(I5)') itime
open(108,file="movie_source/Esum_unit_TD"//trim(adjustl(time))//".txt") 
open(111,file="movie_source/Eex_unit_TD"//trim(adjustl(time))//".txt") 
open(114,file="movie_source/EA_unit_TD"//trim(adjustl(time))//".txt") 
open(117,file="movie_source/EH_unit_TD"//trim(adjustl(time))//".txt") 
open(119,file="movie_source/mag_unit_TD"//trim(adjustl(time))//".txt") 

do iz=nunitzbd(1),nunitzbd(2)
   do iy=nunitybd(2),nunitybd(1),-1
      do ix=nunitxbd(1),nunitxbd(2)
         i=iunit(ix,iy,iz)
         write(108,"(4E20.10e2)") dble(ix-nunitxbd(1)+0.5)*ax/10.0d0,dble(iy-nunitybd(1)+0.5)*ax/10.0d0,dble(iz-nunitzbd(1)+0.5)*cx/10.0d0,&
              engrunit(i,4)+engfeunit(i,4)
         write(111,"(4E20.10e2)") dble(ix-nunitxbd(1)+0.5)*ax/10.0d0,dble(iy-nunitybd(1)+0.5)*ax/10.0d0,dble(iz-nunitzbd(1)+0.5)*cx/10.0d0,&
              engrunit(i,1)+engfeunit(i,1)
         write(114,"(4E20.10e2)") dble(ix-nunitxbd(1)+0.5)*ax/10.0d0,dble(iy-nunitybd(1)+0.5)*ax/10.0d0,dble(iz-nunitzbd(1)+0.5)*cx/10.0d0,&
              engrunit(i,2)+engfeunit(i,2)
         write(117,"(4E20.10e2)") dble(ix-nunitxbd(1)+0.5)*ax/10.0d0,dble(iy-nunitybd(1)+0.5)*ax/10.0d0,dble(iz-nunitzbd(1)+0.5)*cx/10.0d0,&
              engrunit(i,3)+engfeunit(i,3)
         write(119,"(4E20.10e2)") dble(ix-nunitxbd(1)+0.5)*ax/10.0d0,dble(iy-nunitybd(1)+0.5)*ax/10.0d0,dble(iz-nunitzbd(1)+0.5)*cx/10.0d0,&
              magunit(i,3)
      end do
   end do
end do

end subroutine plot_engdens_unit
