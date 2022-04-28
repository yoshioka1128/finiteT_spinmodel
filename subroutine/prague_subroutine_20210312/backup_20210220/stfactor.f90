subroutine stfactor(Stevens,L1,J0,S1,rdcoff)
implicit none
integer(4),parameter :: NDIM=100
integer(4) :: ll,i,IER6
real(8) :: Stevens(6),dl,ite,L1,J0,S1,SIXCOF(NDIM),L1MIN,L1MAX,rdcoff(0:6),abc,def,w6jsym

! Stevens factor
Stevens=0.0d0
do ll=2,6,2
   dl=dble(ll)
   ! factorial
   abc=1.0d0
   do i=1,int(2.0d0*J0-dl)
      abc=abc*dble(i)
   end do
   def=1.0d0
   do i=1,int(2.0d0*J0+dl+1.0d0)
      def=def*dble(i)
   end do
   do ite=0,int(J0*2.0d0)
      if(J0.gt.L1+S1.or.J0.lt.abs(L1-S1)) cycle
      if(J0.gt.J0+dl.or.J0.lt.abs(J0-dl)) cycle
      Stevens(ll)=(-1.0d0)**(L1+S1+J0)*2.0d0**ll*(2.0d0*J0+1.0d0)*dsqrt(abc/def)*w6jsym(L1,L1,dl,J0,J0,S1)*rdcoff(ll)
   end do
end do

end subroutine stfactor
