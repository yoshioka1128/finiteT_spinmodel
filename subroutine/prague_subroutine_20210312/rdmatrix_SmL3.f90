subroutine rdmatrix_SmL3(n4f,rdcoff,L1)
implicit none
integer(4) :: n4f,m,k
real(8) :: abc,def,w3jsym,rdcoff(0:6),L1

! lz of signel particle
if(n4f.ne.5) then
   write(6,*) "error"
   stop
end if

! reduced matrix element
do k=2,6,2
   abc=0.0d0
   def=0.0d0
   m=-3
   abc=abc+(-1.0d0)**m*7.0d0*&
        w3jsym(3.0d0,dble(k),3.0d0,0.0d0,0.0d0)*w3jsym(3.0d0,dble(k),3.0d0,-dble(m),0.0d0)
   do m=0,3 ! lz of single particle
      abc=abc+(-1.0d0)**m*7.0d0*&
           w3jsym(3.0d0,dble(k),3.0d0,0.0d0,0.0d0)*w3jsym(3.0d0,dble(k),3.0d0,-dble(m),0.0d0)
   end do
   abc=2.0d0/3.0d0*abc
   do m=-2,-1
      def=def+(-1.0d0)**m*7.0d0*&
           w3jsym(3.0d0,dble(k),3.0d0,0.0d0,0.0d0)*w3jsym(3.0d0,dble(k),3.0d0,-dble(m),0.0d0)
   end do
   do m=1,3
      def=def+(-1.0d0)**m*7.0d0*&
           w3jsym(3.0d0,dble(k),3.0d0,0.0d0,0.0d0)*w3jsym(3.0d0,dble(k),3.0d0,-dble(m),0.0d0)
   end do
   def=1.0d0/3.0d0*def
   rdcoff(k)=(abc+def)/w3jsym(L1,dble(k),L1,-L1,0.0d0)
end do ! for k

write(6,*) "reduced matrix elements"
do k=2,6,2
   write(6,*) k,rdcoff(k)
end do

return
end subroutine rdmatrix_SmL3
