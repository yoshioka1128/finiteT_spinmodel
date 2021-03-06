subroutine Clm_Sm(L1,S1,J0,ZUst,imax,stfct,nJex)
implicit none
integer(4),parameter :: NDIM=100
real(8),parameter :: pi=dacos(-1.0d0)
complex(kind(0d0)),parameter :: im=(0.0d0,1.0d0)
integer(4) :: i,j, ite, ll, mm,mm2, IER3, IER6, imax, i0(0:100), i1, i2, l, m, jj1, jj2, imaxj,icount
integer(4) :: nJex
real(8) :: L1, L2, L3, L4, L5, L6, M1, M2, M2MIN, M2MAX, THRCOF(NDIM),SIXCOF(NDIM)
real(8) :: J0, J1, J2, dl, dm, S1, abc
real(8) :: L1MIN, L1MAX, gamma, stfct(6,-6:6)
real(8) :: cclm(-6:6,-6:6),coff(6,-6:6),rdcoff(0:6),factorial(2),Ust(imax,imax,6,-6:6)
complex(kind(0d0)) :: ZUst(imax,imax,6,-6:6)

! <||U||> table
! Sm
rdcoff=0.0d0
rdcoff(2)=dsqrt(2.0d0*11.0d0*13.0d0/15.0d0)/3.0d0
rdcoff(4)=dsqrt(2.0d0*13.0d0/11.0d0)*2.0d0/3.0d0
rdcoff(6)=-10.0d0/3.0d0*dsqrt(5.0d0*17.0d0/(3.0d0*11.0d0*13.0d0))

! set alm table
cclm=1.0d0
cclm(2,0)=dsqrt(5.0d0/pi)/4.0d0
cclm(2,2)=dsqrt(15.0d0/pi)/4.0d0
cclm(4,0)=dsqrt(1.0d0/pi)*3.0d0/16.0d0
cclm(4,2)=dsqrt(5.0d0/pi)*3.0d0/8.0d0
cclm(4,4)=dsqrt(35.0d0/pi)*3.0d0/16.0d0
cclm(6,0)=dsqrt(13.0d0/pi)/32.0d0
cclm(6,2)=dsqrt(2730d0/pi)/64.0d0
cclm(6,4)=dsqrt(13.0d0/7.0d0/pi)*21.0d0/32.0d0
cclm(6,6)=dsqrt(26.0d0/231.0d0/pi)*231.0d0/64.0d0
do l=2,6,2
   do m=0,6,2
      cclm(l,-m)=cclm(l,m)
   end do
end do

stfct=0.0d0
do l=2,6,2
   do m=-6,6,2
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

Ust=0.0d0
do jj1=0,nJex
!   write(6,*) "excited state",jj1
   J1=J0+dble(jj1)
   do jj2=0,nJex
      J2=J0+dble(jj2)
!      write(6,*) "J1,J2=",J1,J2
      do ll=2,6,2
         dl=dble(ll)
         do ite=0,int(J1*2.0d0)
            M1=-J1+dble(ite)
! <J1,M1|U|J2,M2>
!***PURPOSE  Evaluate the 3j symbol g(M2) = (L1 L2   L3  )
!                                           (M1 M2 -M1-M2)
!subroutine DRC3JM (L1, L2, L3, M1, M2MIN, M2MAX, THRCOF, NDIM, IER)
!     THRCOF :OUT Set of 3j coefficients generated by evaluating the
!                 3j symbol for all allowed values of M2.  THRCOF(I)
!                 will contain g(M2MIN+I-1), I=1,2,...,M2MAX-M2MIN+1.
            call DRC3JM (J1, J2, dl, -M1, M2MIN, M2MAX, THRCOF, NDIM, IER3)
!***PURPOSE  Evaluate the 6j symbol h(L1) = {L1 L2 L3}
!                                           {L4 L5 L6}
!subroutine DRC6J (L2, L3, L4, L5, L6, L1MIN, L1MAX, SIXCOF, NDIM, IER)
!     SIXCOF :OUT Set of 6j coefficients generated by evaluating the
!                 6j symbol for all allowed values of L1.  SIXCOF(I)
!                 will contain h(L1MIN+I-1), I=1,2,...,L1MAX-L1MIN+1.
            call DRC6J(L1, dl, J1, J2, S1, L1MIN, L1MAX, SIXCOF, NDIM, IER6)
            if(IER3.eq.2) cycle
            do i=1,int(M2MAX-M2MIN)+1
               M2=M2MIN+dble(i-1)
               mm=int(M1-M2)
               if(mod(mm,2).ne.0) cycle
!               write(6,*) "l,m",ll,mm
!               write(6,*) "M1,M2",M1,M2
               i1=int(M1+J1+1)+i0(jj1) ! M1=-J1,-J1+1,...,J1 (i1=1,2,...,2J+1)
               i2=int(M2+J2+1)+i0(jj2)
!               write(6,*) "(i1,i2,mm)=",i1,i2,mm
               Ust(i1,i2,ll,mm)=(-1.0d0)**(L1+S1-M1+J1-J2)*dsqrt((2.0d0*J1+1.0d0)*(2.0d0*J2+1.0d0))*&
                    THRCOF(i)*SIXCOF(int(L1-L1MIN+1))
!               write(6,*) i1,i2,Ust(i1,i2,ll,mm)
            end do
         end do
      end do
   end do
end do

ZUst=(0.0d0,0.0d0)
do l=2,6,2
   do m=0,l,2
      do i2=1,imax
         do i1=1,imax
            if(m.eq.0) then
               ZUst(i1,i2,l,m)=Ust(i1,i2,l,m)
            else if(mod(m/2,2).ne.0) then
               ZUst(i1,i2,l,m)=im*(Ust(i1,i2,l,-m)-Ust(i1,i2,l,m))
            else
               ZUst(i1,i2,l,m)=Ust(i1,i2,l,-m)+Ust(i1,i2,l,m)
            end if
         end do
      end do
   end do
end do

end subroutine Clm_Sm
