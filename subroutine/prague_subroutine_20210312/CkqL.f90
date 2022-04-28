subroutine CkqL(L1,S1,J0,Ust,imax,nJex)
implicit none
integer(4),parameter :: NDIM=100
real(8),parameter :: pi=dacos(-1.0d0)
complex(kind(0d0)),parameter :: im=(0.0d0,1.0d0)
integer(4) :: k,i,j, ite, ll, mm, imax, i0(0:100), i1, i2, l, m, jj1, jj2, nJex
real(8) :: L1, L4, M1, M2, M2MIN, M2MAX
real(8) :: J0, J1, J2, dl, S1, abc,def,w3jsym,w6jsym, Ust(imax,imax,6,-6:6)

i0=0
do i=1,nJex+1
   i0(i)=i*int(J0*2.0d0+1.0d0)+i*(i-1) ! 0,(2*J0+1),2*(2J0+1)+2*1,3*(2J0+1)+2*3,.....,i*(2*J0+1)+i*(i-1)
end do

! matrix element of spherical tensor operator
Ust=0.0d0
do jj1=0,nJex ! for excited J multiplet
   J1=J0+dble(jj1)
   do jj2=0,nJex ! for excited Jmultiplet
      J2=J0+dble(jj2)
      do ll=1,6
         dl=dble(ll)
         abc=w6jsym(L1,L1,dl,J1,J2,S1)
         do ite=0,int(J1*2.0d0)
            M1=-J1+dble(ite)
            M2MIN=MAX(-J2,-dl+M1)
            M2MAX=MIN(J2,dl+M1)
            do i=1,int(M2MAX-M2MIN)+1 ! for M2
               M2=M2MIN+dble(i-1)
               mm=int(M1-M2)
               i1=int(M1+J1+1)+i0(jj1) ! M1=-J1,-J1+1,...,J1 (i1=1,2,...,2J+1)
               i2=int(M2+J2+1)+i0(jj2)
               Ust(i1,i2,ll,mm)=(-1.0d0)**(L1+S1-M1+J1-J2)*dsqrt((2.0d0*J1+1.0d0)*(2.0d0*J2+1.0d0))*w3jsym(J1,J2,dl,-M1,M2)*abc
            end do
         end do
      end do
   end do
end do

return

end subroutine CkqL
