subroutine plot_engdens_unit(itime,nunit,engrunit,engfeunit,magunit,nunitx,nunity,nunitz,iunit,ax,cx)
implicit none
integer(4) :: nunit,itime,nunitx,nunity,nunitz,ix,iy,iz,ir,iunit(nunitx,nunity,nunitz),i
real(8) :: engrunit(nunit,4),engfeunit(nunit,4),ax,cx,magunit(nunit,3)
character(5) :: time

write(time,'(I5)') itime
open(108,file="movie_source/Esum_unit_TD"//trim(adjustl(time))//".txt") 
open(111,file="movie_source/Eex_unit_TD"//trim(adjustl(time))//".txt") 
open(114,file="movie_source/EA_unit_TD"//trim(adjustl(time))//".txt") 
open(117,file="movie_source/EH_unit_TD"//trim(adjustl(time))//".txt") 
open(119,file="movie_source/mag_unit_TD"//trim(adjustl(time))//".txt") 

do iz=1,nunitz
   do iy=nunity,1,-1
      do ix=1,nunitx
         i=iunit(ix,iy,iz)
         write(108,"(4E20.10e2)") dble(ix-1+0.5)*ax/10.0d0,dble(iy-1+0.5)*ax/10.0d0,dble(iz-1+0.5)*cx/10.0d0,&
              engrunit(i,4)+engfeunit(i,4)
         write(111,"(4E20.10e2)") dble(ix-1+0.5)*ax/10.0d0,dble(iy-1+0.5)*ax/10.0d0,dble(iz-1+0.5)*cx/10.0d0,&
              engrunit(i,1)+engfeunit(i,1)
         write(114,"(4E20.10e2)") dble(ix-1+0.5)*ax/10.0d0,dble(iy-1+0.5)*ax/10.0d0,dble(iz-1+0.5)*cx/10.0d0,&
              engrunit(i,2)+engfeunit(i,2)
         write(117,"(4E20.10e2)") dble(ix-1+0.5)*ax/10.0d0,dble(iy-1+0.5)*ax/10.0d0,dble(iz-1+0.5)*cx/10.0d0,&
              engrunit(i,3)+engfeunit(i,3)
         write(119,"(4E20.10e2)") dble(ix-1+0.5)*ax/10.0d0,dble(iy-1+0.5)*ax/10.0d0,dble(iz-1+0.5)*cx/10.0d0,&
              magunit(i,3)
      end do
   end do
end do

end subroutine plot_engdens_unit
