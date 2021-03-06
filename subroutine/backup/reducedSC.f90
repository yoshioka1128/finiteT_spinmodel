subroutine reducedSC(L1,S1,J1,ZUst,imax,stfct,nJex,n4f,rdsc,HexT,nr,nrmax,lambda)
implicit none
integer(4),parameter :: NDIM=100
real(8),parameter :: pi=dacos(-1.0d0)
complex(kind(0d0)),parameter :: im=(0.0d0,1.0d0)
integer(4) :: k,i,j, ite, ll, mm, IER3, IER6, imax, i1, i2, l, m, jj1, jj2, imaxj,icount
integer(4) :: nJex,mm1,n4f,nrmax,nr,ir
real(8) :: L1, L2, L3, L4, L5, L6, M1, M2, M2MIN, M2MAX, THRCOF(NDIM),SIXCOF(NDIM), dlkl(3)
real(8) :: J1, J2, dl, dm, S1, abc,def,lambda
real(8) :: L1MIN, L1MAX, gamma, stfct(6,-6:6),rdsc(6,nrmax),HexT(nrmax)
real(8) :: cclm(-6:6,-6:6),coff(6,-6:6),rdcoff(0:6),factorial(2),Ust(imax,imax,6,-6:6)
complex(kind(0d0)) :: ZUst(imax,imax,6,-6:6)

! lz of signel particle
if(n4f.lt.7) mm1=4-n4f
if(n4f.gt.7) mm1=11-n4f

! reduced matrix element
do k=2,6,2
   abc=0.0d0
   do m=mm1,3 ! lz of single particle
! (3, k, 3)
! (0, 0, 0)
      call DRC3JM (3.0d0, dble(k), 3.0d0, 0.0d0, M2MIN, M2MAX, THRCOF, NDIM, IER3) 
      do i=1,int(M2MAX-M2MIN)+1
         M2=M2MIN+dble(i-1)
         if(M2.eq.0.0d0) dlkl(1)=THRCOF(i)
      end do
! (3, k, 3)
! (-m, 0, m)
      call DRC3JM (3.0d0, dble(k), 3.0d0, -dble(m), M2MIN, M2MAX, THRCOF, NDIM, IER3)
      do i=1,int(M2MAX-M2MIN)+1
         M2=M2MIN+dble(i-1)
         if(M2.eq.0.0d0) dlkl(2)=THRCOF(i)
      end do
      abc=abc+(-1.0d0)**m*7.0d0*dlkl(1)*dlkl(2)
   end do
! (L1,  k, L1)
! (-L1, 0, L1)
   call DRC3JM (L1, dble(k), L1, -L1, M2MIN, M2MAX, THRCOF, NDIM, IER3)
   do i=1,int(M2MAX-M2MIN)+1
      M2=M2MIN+dble(i-1)
      if(M2.eq.0.0d0) dlkl(3)=THRCOF(i)
   end do
   rdcoff(k)=abc/dlkl(3)
end do ! for k

! factor
rdsc=0.0d0
do ir=1,nr
   do ll=2,6,2
      abc=2.0d0*J1+dble(ll)+2.0d0
      do i=1,2*ll
         abc=abc*(2.0d0*J1+dble(ll)+2.0d0-dble(i))
      end do
      dl=dble(ll)
      ! {L1 J1+1 S1}
      ! {J1  L1  ll}
      call DRC6J(J1+1, S1, J1, L1, dl, L1MIN, L1MAX, SIXCOF, NDIM, IER6)

      rdsc(ll,ir)=2.0d0*HexT(ir)/(lambda*(J1+1.0d0))*dsqrt(dble(ll*(ll+1)))/(2.0d0*dl+1.0d0)*2.0d0**ll
      rdsc(ll,ir)=rdsc(ll,ir)*dsqrt(S1*(L1+1.0d0)*(2.0d0*J1+1.0d0))/(J1+1.0d0)/dsqrt(abc)
      rdsc(ll,ir)=rdsc(ll,ir)*SIXCOF(int(L1-L1MIN+1))*rdcoff(ll) 

!      rdsc(ll,ir)=2.0d0*HexT(ir)/(lambda*(J1+1.0d0))*dsqrt(dble(ll*(ll+1)))/(2.0d0*dl+1.0d0)*2.0d0**ll*dsqrt(2.0d0)
!      rdsc(ll,ir)=rdsc(ll,ir)*dsqrt(S1*(L1+1.0d0)*(2.0d0*J1+1.0d0)/(J1+1.0d0)) ! <J||S^1||J>
!      rdsc(ll,ir)=rdsc(ll,ir)*dsqrt((2.0d0*J1+3.0d0)*(2.0d0*J1+1.0d0))*SIXCOF(int(L1-L1MIN+1))*rdcoff(ll) !<J+1||C^n||J> check OK !
!      rdsc(ll,ir)=rdsc(ll,ir)/dsqrt((2.0d0*J1+1.0d0)*(2.0d0*J1+2.0d0)*(2.0d0*J1+3.0d0))/dsqrt(abc)
   end do
end do

end subroutine reducedSC
