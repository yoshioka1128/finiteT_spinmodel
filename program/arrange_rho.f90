program arrange_rho
real(8) :: ax,ay,rho(5)
integer(4) :: nx,ny,ix,iy,i,len,status
character(20) :: filename,dir

call getenv( 'filename', dir )
write(*,*) dir

open(10,file=""//trim(adjustl(dir))//".rhosum_110") ! from edit lapw5.def for unit 9 and execute "lapw5 lapw5.def"
open(11,file=""//trim(adjustl(dir))//".arrrho_editsum_110")

call plot

open(10,file=""//trim(adjustl(dir))//".rhoup_110")
open(11,file=""//trim(adjustl(dir))//".arrrho_editup_110")
call plot

open(10,file=""//trim(adjustl(dir))//".rhodn_110")
open(11,file=""//trim(adjustl(dir))//".arrrho_editdn_110")
call plot

! 001_1
open(10,file=""//trim(adjustl(dir))//".rhosum_001_1") ! from edit lapw5.def for unit 9 and execute "lapw5 lapw5.def"
open(11,file=""//trim(adjustl(dir))//".arrrho_editsum_001_1")

call plot

open(10,file=""//trim(adjustl(dir))//".rhoup_001_1")
open(11,file=""//trim(adjustl(dir))//".arrrho_editup_001_1")
call plot

open(10,file=""//trim(adjustl(dir))//".rhodn_001_1")
open(11,file=""//trim(adjustl(dir))//".arrrho_editdn_001_1")
call plot

! 001_1_center
open(10,file=""//trim(adjustl(dir))//".rhosum_001_1_center") ! from edit lapw5.def for unit 9 and execute "lapw5 lapw5.def"
open(11,file=""//trim(adjustl(dir))//".arrrho_editsum_001_1_center")
call plot_center

open(10,file=""//trim(adjustl(dir))//".rhoup_001_1_center")
open(11,file=""//trim(adjustl(dir))//".arrrho_editup_001_1_center")
call plot_center

open(10,file=""//trim(adjustl(dir))//".rhodn_001_1_center")
open(11,file=""//trim(adjustl(dir))//".arrrho_editdn_001_1_center")
call plot_center

! 001_2
open(10,file=""//trim(adjustl(dir))//".rhosum_001_2") ! from edit lapw5.def for unit 9 and execute "lapw5 lapw5.def"
open(11,file=""//trim(adjustl(dir))//".arrrho_editsum_001_2")

call plot

open(10,file=""//trim(adjustl(dir))//".rhoup_001_2")
open(11,file=""//trim(adjustl(dir))//".arrrho_editup_001_2")
call plot

open(10,file=""//trim(adjustl(dir))//".rhodn_001_2")
open(11,file=""//trim(adjustl(dir))//".arrrho_editdn_001_2")
call plot

! 100_2_center
open(10,file=""//trim(adjustl(dir))//".rhosum_100_2") ! from edit lapw5.def for unit 9 and execute "lapw5 lapw5.def"
open(11,file=""//trim(adjustl(dir))//".arrrho_editsum_100_2")
call plot

open(10,file=""//trim(adjustl(dir))//".rhoup_100_2")
open(11,file=""//trim(adjustl(dir))//".arrrho_editup_100_2")
call plot

open(10,file=""//trim(adjustl(dir))//".rhodn_100_2")
open(11,file=""//trim(adjustl(dir))//".arrrho_editdn_100_2")
call plot

! 100_1_center
open(10,file=""//trim(adjustl(dir))//".rhosum_100_1_center") ! from edit lapw5.def for unit 9 and execute "lapw5 lapw5.def"
open(11,file=""//trim(adjustl(dir))//".arrrho_editsum_100_1_center")
call plot_center

open(10,file=""//trim(adjustl(dir))//".rhoup_100_1_center")
open(11,file=""//trim(adjustl(dir))//".arrrho_editup_100_1_center")
call plot_center

open(10,file=""//trim(adjustl(dir))//".rhodn_100_1_center")
open(11,file=""//trim(adjustl(dir))//".arrrho_editdn_100_1_center")
call plot_center



end program arrange_rho


subroutine plot
integer(4) :: ix,iy,i,j,nx,ny
real(8) :: rho(5),ax,ay
real(8),parameter :: Bohr=0.529177d0

read(10,*) nx,ny,ax,ay
ax=ax*Bohr
ay=ay*Bohr
do ix=1,nx
   iy=0
   do i=1,ny/5
      read(10,*) rho(1),rho(2),rho(3),rho(4),rho(5)
      do j=1,5
         iy=iy+1
!         write(6,*) ax*dble(ix-1)/dble(nx-1),ay*dble(iy-1)/dble(ny-1),rho(j)
         write(11,*) ax*dble(ix-1)/dble(nx-1),ay*dble(iy-1)/dble(ny-1),rho(j)
!         write(11,*) dble(ix-1)/dble(nx-1),dble(iy-1)/dble(ny-1),rho(j)
      end do
   end do
   write(11,*) 
end do 

end subroutine plot


subroutine plot_center
integer(4) :: ix,iy,i,j,nx,ny,ix0,iy0
real(8) :: rho(5),ax,ay
real(8),parameter :: Bohr=0.529177d0
real(8),allocatable :: box(:,:)

read(10,*) nx,ny,ax,ay
allocate(box(nx,ny))
ax=ax*Bohr
ay=ay*Bohr
do ix0=1,nx
   iy0=0
   do i=1,ny/5
      read(10,*) rho(1),rho(2),rho(3),rho(4),rho(5)
      do j=1,5
         iy0=iy0+1
         if(ix0.le.nx/2) then
            ix=ix0+nx/2
            if(iy0.le.ny/2) then
               iy=iy0+ny/2
               box(ix,iy)=rho(j)
            else
               iy=iy0-ny/2
               box(ix,iy)=rho(j)
            end if
         else
            ix=ix0-nx/2
            if(iy0.le.ny/2) then
               iy=iy0+ny/2
               box(ix,iy)=rho(j)
            else
               iy=iy0-ny/2
               box(ix,iy)=rho(j)
            end if
         end if
      end do
   end do
end do 

do ix=1,nx
   do iy=1,ny
      write(11,*) ax*dble(ix-1)/dble(nx-1),ay*dble(iy-1)/dble(ny-1),box(ix,iy)
   end do
   write(11,*)
end do

end subroutine plot_center

