subroutine rdmatrix(n4f,rdcoff,L1)
implicit none
integer(4) :: n4f,m,mm1,k
real(8) :: abc,w3jsym,rdcoff(0:6),L1

! lz of signel particle
if(n4f.lt.7) mm1=4-n4f
if(n4f.gt.7) mm1=11-n4f

! reduced matrix element
do k=2,6,2
   abc=0.0d0
   do m=mm1,3 ! lz of single particle
      abc=abc+(-1.0d0)**m*7.0d0*w3jsym(3.0d0,dble(k),3.0d0,0.0d0,0.0d0)*w3jsym(3.0d0,dble(k),3.0d0,-dble(m),0.0d0)
   end do
   rdcoff(k)=abc/w3jsym(L1,dble(k),L1,-L1,0.0d0)
end do ! for k

return
end subroutine rdmatrix
