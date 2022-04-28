subroutine xifactor(xi,L1,J1,S1,rdcoff)
implicit none
integer(4),parameter :: NDIM=100
integer(4) :: ll,i,IER6
real(8) :: xi(6),dl,L1,J1,S1,SIXCOF(NDIM),L1MIN,L1MAX,rdcoff(0:6),abc,def,w6jsym,rdjjp

! xi factor
xi=0.0d0
do ll=2,6,2
   dl=dble(ll)
   ! factorial
   abc=1.0d0
   do i=int(2.0d0*J1+dl+1),int(2.0d0*J1-dl+2),-1
      abc=abc*dble(i)
   end do
   rdjjp=(-1)**(int(S1+L1+J1+ll))*dsqrt((2.0d0*J1+1.0d0)*(2.0d0*J1+3.0d0))*w6jsym(L1,J1,S1,J1+1.0d0,L1,dl)
   xi(ll)=rdjjp*2.0d0**ll*dsqrt((L1+1.0d0)/S1)*dsqrt((2.0d0*J1+dl+2.0d0)/((2.0d0*J1+3.0d0)*dl*(dl+1.0d0)*abc))*rdcoff(ll)
end do

end subroutine xifactor
