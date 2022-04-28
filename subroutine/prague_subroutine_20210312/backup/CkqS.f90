subroutine CkqS(L1,S1,J0,Vst,imax,nJex)
implicit none
real(8),parameter :: pi=dacos(-1.0d0)
complex(kind(0d0)),parameter :: im=(0.0d0,1.0d0)
integer(4) :: k,i,j, ite,  mm, imax, i0(0:100), i1, i2, l, m, jj1, jj2, nJex
real(8) :: L1, L4, M1, M2, M2MIN, M2MAX
real(8) :: J0, J1, J2, dl, S1, abc,def,w3jsym,w6jsym, Vst(imax,imax,-1:1)

i0=0
do i=1,nJex+1
   i0(i)=i*int(J0*2.0d0+1.0d0)+i*(i-1)
end do
i0(0)=0
! matrix element of spherical tensor operator
Vst=0.0d0
do jj1=0,nJex
   J1=J0+dble(jj1)
   do jj2=0,nJex
      J2=J0+dble(jj2)
      dl=1.0d0
      do ite=0,int(J1*2.0d0)
         M1=-J1+dble(ite)
         M2MIN=MAX(-J2,-dl+M1)
         M2MAX=MIN(J2,dl+M1)
         do i=1,int(M2MAX-M2MIN)+1
            M2=M2MIN+dble(i-1)
            mm=int(M1-M2)
            i1=int(M1+J1+1)+i0(jj1) ! M1=-J1,-J1+1,...,J1 (i1=1,2,...,2J+1)
            i2=int(M2+J2+1)+i0(jj2)
            Vst(i1,i2,mm)=(-1.0d0)**(L1+S1-M1)*dsqrt((2.0d0*J1+1.0d0)*(2.0d0*J2+1.0d0))*&
                 w3jsym(J1,J2,dl,-M1,M2)*w6jsym(S1,S1,dl,J1,J2,L1)
         end do
      end do
   end do
end do

end subroutine CkqS
