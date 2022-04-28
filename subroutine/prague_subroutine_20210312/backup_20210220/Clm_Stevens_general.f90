subroutine Clm_Stevens_general(L1,S1,J0,ZUst,imax,stfct,nJex,n4f,Stevens)
implicit none
integer(4),parameter :: NDIM=100
real(8),parameter :: pi=dacos(-1.0d0)
complex(kind(0d0)),parameter :: im=(0.0d0,1.0d0)
integer(4) :: k,i,j, ite, ll, mm, IER3, IER6, imax, i0(0:100), i1, i2, l, m, jj1, jj2, imaxj,icount
integer(4) :: nJex,mm1,n4f
real(8) :: L1, L2, L3, L4, L5, L6, M1, M2, M2MIN, M2MAX, THRCOF(NDIM),SIXCOF(NDIM), dlkl(3)
real(8) :: J0, J1, J2, dl, dm, S1, abc,def
real(8) :: L1MIN, L1MAX, gamma, stfct(6,-6:6),Stevens(6),w6jsym,w3jsym
real(8) :: cclm(-6:6,-6:6),coff(6,-6:6),rdcoff(0:6),factorial(2),Ust(imax,imax,6,-6:6)
complex(kind(0d0)) :: ZUst(imax,imax,6,-6:6)

call rdmatrix(n4f,rdcoff,L1)

! set clm table
cclm=1.0d0
cclm(2,0)=dsqrt(5.0d0/pi)/4.0d0
cclm(2,1)=dsqrt(15.0d0/pi)/2.0d0
cclm(2,2)=dsqrt(15.0d0/pi)/4.0d0
cclm(4,0)=dsqrt(1.0d0/pi)*3.0d0/16.0d0
cclm(4,1)=dsqrt(10.0d0/pi)*3.0d0/8.0d0
cclm(4,2)=dsqrt(5.0d0/pi)*3.0d0/8.0d0
cclm(4,3)=dsqrt(70.0d0/pi)*3.0d0/8.0d0
cclm(4,4)=dsqrt(35.0d0/pi)*3.0d0/16.0d0
cclm(6,0)=dsqrt(13.0d0/pi)/32.0d0
cclm(6,1)=dsqrt(273.0d0/pi)/16.0d0
cclm(6,2)=dsqrt(2730d0/pi)/64.0d0
cclm(6,3)=dsqrt(2730d0/pi)/32.0d0
cclm(6,4)=dsqrt(13.0d0/7.0d0/pi)*21.0d0/32.0d0
cclm(6,5)=dsqrt(2002.0d0/pi)*3.0d0/32.0d0
cclm(6,6)=dsqrt(26.0d0/231.0d0/pi)*231.0d0/64.0d0

do l=2,6,2
   do m=0,6
      cclm(l,-m)=cclm(l,m)
   end do
end do

! Yamada (1988) Eq. (19)
stfct=0.0d0
do l=2,6,2
   do m=-6,6
      if(m.eq.0) then
         stfct(l,m)=1.0d0/cclm(l,m)*dsqrt((2.0d0*dble(l)+1.0d0)/(4.0d0*pi))*rdcoff(l)
      else
         stfct(l,m)=1.0d0/cclm(l,m)*dsqrt((2.0d0*dble(l)+1.0d0)/(8.0d0*pi))*rdcoff(l)
      end if
   end do
end do

i0=0
do i=1,nJex+1
   i0(i)=i*int(J0*2.0d0+1.0d0)+i*(i-1)
end do
i0(0)=0

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
      M1=-J0+dble(ite)
      call DRC6J(L1, dl, J0, J0, S1, L1MIN, L1MAX, SIXCOF, NDIM, IER6)
      Stevens(ll)=(-1.0d0)**(L1+S1+J0)*2**ll*(2.0d0*J0+1.0d0)*dsqrt(abc/def)*SIXCOF(int(L1-L1MIN+1))*rdcoff(ll)
   end do
end do

! matrix element of spherical tensor operator
Ust=0.0d0
do jj1=0,nJex
   J1=J0+dble(jj1)
   do jj2=0,nJex
      J2=J0+dble(jj2)
      do ll=1,6
         dl=dble(ll)
         abc=w6jsym(L1,L1,dl,J1,J2,S1)
         do ite=0,int(J1*2.0d0)
            M1=-J1+dble(ite)
            M2MIN=MAX(-J2,-dl+M1)
            M2MAX=MIN(J2,dl+M1)
            do i=1,int(M2MAX-M2MIN)+1
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

ZUst=(0.0d0,0.0d0)
do l=2,6,2
   do m=-l,l
      do i2=1,imax
         do i1=1,imax
            if(m.eq.0) then
               ZUst(i1,i2,l,m)=Ust(i1,i2,l,m)
            else if(m.lt.0) then
               ZUst(i1,i2,l,m)=im*(Ust(i1,i2,l,-abs(m))-(-1)**m*Ust(i1,i2,l,abs(m)))
            else if(m.gt.0) then
               ZUst(i1,i2,l,m)=Ust(i1,i2,l,-abs(m))+(-1)**m*Ust(i1,i2,l,abs(m))
            end if
         end do
      end do
   end do
end do

end subroutine Clm_Stevens_general
