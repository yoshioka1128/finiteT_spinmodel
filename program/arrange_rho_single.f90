program arrange_rho
real(8) :: ax,ay,rho(5)
integer(4) :: nx,ny,ix,iy,i,len,status
character(30) :: filename,dir

!write(6,*) "directry name ?"
!read(5,*) dir
call getenv( 'filename', dir )
write(*,*) dir

open(10,file=""//trim(adjustl(dir))//".rho") ! from edit lapw5.def for unit 9 and execute "lapw5 lapw5.def"
open(11,file=""//trim(adjustl(dir))//".rho_arrange")

call plot

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
         write(11,*) ax*dble(ix-1)/dble(nx-1),ay*dble(iy-1)/dble(ny-1),rho(j)
      end do
   end do
   write(11,*) 
end do 

end subroutine plot


