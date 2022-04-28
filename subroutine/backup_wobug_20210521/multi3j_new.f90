subroutine multi3j(L1,S1,J1,ZUst,imax,stfct,nJex,n4f,rdscn,HexT,nr,nrmax,lambda,d3j,t3j,q3j,eigen,g,JM,temp)
implicit none
integer(4),parameter :: NDIM=100
real(8),parameter :: pi=dacos(-1.0d0)
complex(kind(0d0)),parameter :: im=(0.0d0,1.0d0)
integer(4) :: k,i,j, ii, in,ite, ll, mm, IER3, IER6, imax, i1, i2, l, m, jj1, jj2, imaxj,icount
integer(4) :: nJex,mm1,n4f,nrmax,nr,ir,nJ1
real(8) :: L1, L2, L3, L4, L5, L6, M1, M2, M2MIN, M2MAX, THRCOF(NDIM),SIXCOF(NDIM), dlkl(3)
real(8) :: J1, J2, dm, S1, Jin, abc,def,ghi,jkl,mno,lambda,dz,delta(2)
real(8) :: L1MIN, L1MAX, gamma, stfct(6,-6:6),rdscn(6,nrmax,4),HexT(nrmax),d3j(6,nrmax),t3j(6,nrmax)
real(8) :: q3j(6,nrmax),d3jsym
real(8) :: cclm(-6:6,-6:6),coff(6,-6:6),factorial(2),Ust(imax,imax,6,-6:6)
real(8) :: eigen(imax),g(100),JM(100,2),temp
real(8) :: rdcll(0:6),rdcjj(0:6,2),rdcjjp(0:6,2),rdsss,rdsjj(2),rdsjjp(2)
complex(kind(0d0)) :: ZUst(imax,imax,6,-6:6)

! lz of signel particle
if(n4f.lt.7) mm1=4-n4f
if(n4f.gt.7) mm1=11-n4f

delta(1)=lambda*(J1+1.0d0)
delta(2)=lambda*(J1+2.0d0)

! reduced matrix element
do k=2,6,2
   abc=0.0d0
   do m=mm1,3 ! lz of single particle
! (3, k, 3)
! (0, 0, 0)
      call DRC3JM (3.0d0, dble(k), 3.0d0, 0.0d0, M2MIN, M2MAX, THRCOF, NDIM, IER3) 
      if(IER3.eq.2) THRCOF=0.0d0
      do i=1,int(M2MAX-M2MIN)+1
         M2=M2MIN+dble(i-1)
         if(M2.eq.0.0d0) dlkl(1)=THRCOF(i)
      end do
! (3, k, 3)
! (-m, 0, m)
      call DRC3JM (3.0d0, dble(k), 3.0d0, -dble(m), M2MIN, M2MAX, THRCOF, NDIM, IER3)
      if(IER3.eq.2) THRCOF=0.0d0
      do i=1,int(M2MAX-M2MIN)+1
         M2=M2MIN+dble(i-1)
         if(M2.eq.0.0d0) dlkl(2)=THRCOF(i)
      end do
      abc=abc+(-1.0d0)**m*7.0d0*dlkl(1)*dlkl(2)
   end do
! (L1,  k, L1)
! (-L1, 0, L1)
   call DRC3JM (L1, dble(k), L1, -L1, M2MIN, M2MAX, THRCOF, NDIM, IER3)
   if(IER3.eq.2) THRCOF=0.0d0
   do i=1,int(M2MAX-M2MIN)+1
      M2=M2MIN+dble(i-1)
      if(M2.eq.0.0d0) dlkl(3)=THRCOF(i)
   end do
   rdcll(k)=abc/dlkl(3)
end do ! for k
rdsss=dsqrt(S1*(1.0d0+S1)*(1.0d0+2.0d0*S1))

nJ1=int(J1*2.0d0)+1
do i=1,2 
   Jin=J1+dble(i-1)
   ! <J||Sz||J>
   call DRC6J(Jin, L1, Jin, S1, 1.0d0, L1MIN, L1MAX, SIXCOF, NDIM, IER6)
   if(IER6.eq.2) SIXCOF=0.0d0
   rdsjj(i)=(-1.0d0)**(int(L1+S1+Jin)+1)*(2.0d0*Jin+1.0d0)*SIXCOF(int(S1-L1MIN+1.0d0))*rdsss
   ! <J||Sz||J+1>
   call DRC6J(Jin, L1, Jin+1.0d0, S1, 1.0d0, L1MIN, L1MAX, SIXCOF, NDIM, IER6)
   if(IER6.eq.2) SIXCOF=0.0d0
   rdsjjp(i)=(-1.0d0)**(int(L1+S1+Jin))*dsqrt((2.0d0*Jin+1.0d0)*(2.0d0*Jin+3.0d0))*SIXCOF(int(S1-L1MIN+1.0d0))*rdsss
   do in=2,6,2
      ! <J||Cn||J>
      call DRC6J(Jin, S1, Jin, L1, dble(in), L1MIN, L1MAX, SIXCOF, NDIM, IER6)
      if(IER6.eq.2) SIXCOF=0.0d0
      rdcjj(in,i)=(-1.0d0)**(int(L1+S1+Jin)+in)*(2.0d0*Jin+1.0d0)*SIXCOF(int(L1-L1MIN+1.0d0))*rdcll(in)
      ! <J||Cn||J+1>
      call DRC6J(Jin, S1, Jin+1.0d0, L1, dble(in), L1MIN, L1MAX, SIXCOF, NDIM, IER6)
      if(IER6.eq.2) SIXCOF=0.0d0
      rdcjjp(in,i)=(-1.0d0)**(int(L1+S1+Jin)+in)*dsqrt((2.0d0*Jin+1.0d0)*(2.0d0*Jin+3.0d0))*SIXCOF(int(L1-L1MIN+1.0d0))*rdcll(in)
   end do
end do


! 2nd order
d3j=0.0d0
do ir=1,nr
   do ii=1,imax
      eigen(ii)=2.0d0*(g(ii)-1.0d0)*HexT(ir)*JM(ii,2)
   end do
   do in=2,6,2
      abc=4.0d0*HexT(ir)/delta(1)*rdsjjp(1)*(-rdcjjp(in,1))
      dz=0.0d0
      do ii=1,nJ1
         dz=dz+dexp(-eigen(ii)/temp)

!         call DRC3JM(J1,      1.0d0,   J1+1.0d0,-JM(ii,2),M2MIN, M2MAX, THRCOF, NDIM, IER3)
!         if(IER3.eq.2) THRCOF=0.0d0
!         def=THRCOF(int(-M2MIN+1))
         def=d3jsym(J1,1.0d0,J1+1.0d0,-JM(ii,2),0.0d0)
!         call DRC3JM(J1+1.0d0,dble(in),J1      ,-JM(ii,2),M2MIN, M2MAX, THRCOF, NDIM, IER3)
!         if(IER3.eq.2) THRCOF=0.0d0
!         def=def*THRCOF(int(-M2MIN+1))
         def=def*d3jsym(J1+1.0d0,dble(in),J1,-JM(ii,2),0.0d0)
         ! if M2=0, i=int(M2MIN)+1
!         d3j(in,ir)=d3j(in,ir)+abc*def
         if(temp.eq.0.0d0) then
            d3j(in,ir)=abc*def
            exit
         else
            d3j(in,ir)=d3j(in,ir)+abc*def*dexp(-eigen(ii)/temp)
         end if
      end do ! do ii
      if(temp.ne.0.0d0) d3j(in,ir)=d3j(in,ir)/dz

!      write(6,*) "d3j",d3j(in,ir),abc*def,dexp(-eigen(ii)/temp)
   end do ! do in
end do ! do ir
!stop


! 3rd order
rdscn=0.0d0
t3j=0.0d0
do ir=1,nr
   do ii=1,imax
      eigen(ii)=2.0d0*(g(ii)-1.0d0)*HexT(ir)*JM(ii,2)
   end do
   do in=2,6,2
      jkl=(2.0d0*HexT(ir)/delta(1))**2*2.0d0*(-rdsjjp(1)*rdcjjp(in,1)) ! 1-2
      mno=(2.0d0*HexT(ir)/delta(1))**2*(-rdsjjp(1)**2) ! 3-4
      rdscn(in,ir,1)=jkl*rdsjj(1)
      rdscn(in,ir,2)=jkl*rdsjj(2)
      rdscn(in,ir,3)=mno*rdcjj(in,1)
      rdscn(in,ir,4)=mno*rdcjj(in,2)

      dz=0.0d0
      do ii=1,nJ1
         dz=dz+dexp(-eigen(ii)/temp)
         ! 1-2
         call DRC3JM(J1,      1.0d0,   J1+1.0d0,-JM(ii,2),M2MIN, M2MAX, THRCOF, NDIM, IER3)
         if(IER3.eq.2) THRCOF=0.0d0
         abc=THRCOF(int(-M2MIN+1))
         call DRC3JM(J1+1.0d0,dble(in),J1      ,-JM(ii,2),M2MIN, M2MAX, THRCOF, NDIM, IER3)
         if(IER3.eq.2) THRCOF=0.0d0
         abc=abc*THRCOF(int(-M2MIN+1))
         ! 3-4
         call DRC3JM(J1,      1.0d0,   J1+1.0d0,-JM(ii,2),M2MIN, M2MAX, THRCOF, NDIM, IER3)
         if(IER3.eq.2) THRCOF=0.0d0
         def=THRCOF(int(-M2MIN+1))**2

         ! 1
         call DRC3JM(J1,1.0d0,J1,-JM(ii,2),M2MIN, M2MAX, THRCOF, NDIM, IER3)
         if(IER3.eq.2) THRCOF=0.0d0
         ghi=abc*rdscn(in,ir,1)*THRCOF(int(-M2MIN+1))
         ! 2
         call DRC3JM(J1+1.0d0,1.0d0,J1+1.0d0,-JM(ii,2),M2MIN, M2MAX, THRCOF, NDIM, IER3)
         if(IER3.eq.2) THRCOF=0.0d0
         ghi=ghi+abc*rdscn(in,ir,2)*THRCOF(int(-M2MIN+1))
         ! 3
         call DRC3JM(J1,dble(in),J1,-JM(ii,2),M2MIN, M2MAX, THRCOF, NDIM, IER3)
         if(IER3.eq.2) THRCOF=0.0d0
         ghi=ghi+def*rdscn(in,ir,3)*THRCOF(int(-M2MIN+1))
         ! 4
         call DRC3JM(J1+1.0d0,dble(in),J1+1.0d0,-JM(ii,2),M2MIN, M2MAX, THRCOF, NDIM, IER3)
         if(IER3.eq.2) THRCOF=0.0d0
         ghi=ghi+def*rdscn(in,ir,4)*THRCOF(int(-M2MIN+1))

         if(temp.eq.0.0d0) then
            t3j(in,ir)=(-1.0d0)**(int(J1-JM(ii,2)))*ghi
            exit
         else
            t3j(in,ir)=t3j(in,ir)+(-1.0d0)**(int(J1-JM(ii,2)))*ghi*dexp(-eigen(ii)/temp)
         end if
         
      end do ! do ii
      if(temp.ne.0.0d0) then
         t3j(in,ir)=t3j(in,ir)/dz
      end if

   end do ! do in
end do ! do ir



! 4th order
rdscn=0.0d0
q3j=0.0d0
do ir=1,nr
   do ii=1,imax
      eigen(ii)=2.0d0*(g(ii)-1.0d0)*HexT(ir)*JM(ii,2)
   end do
   do in=2,6,2
      jkl=-3.0d0/2.0d0*(2.0d0*HexT(ir))**3/delta(1)**2/delta(2) !1-2
      mno=3.0d0*(2.0d0*HexT(ir))**3/delta(1)**3 !3-4
      rdscn(in,ir,1)=jkl*rdcjjp(in,1)*(-rdsjjp(1))*rdsjjp(2)**2
      rdscn(in,ir,2)=jkl*rdcjjp(in,2)*(-rdsjjp(2))*rdsjjP(1)**2
      rdscn(in,ir,3)=mno*rdcjjp(in,1)*(-rdsjjp(1))*rdsjjp(1)**2

      dz=0.0d0
      do ii=1,nJ1
         dz=dz+dexp(-eigen(ii)/temp)
         ! 1
         call DRC3JM(J1,1.0d0,J1+1.0d0,-JM(ii,2),M2MIN, M2MAX, THRCOF, NDIM, IER3)
         if(IER3.eq.2) THRCOF=0.0d0
         abc=THRCOF(int(-M2MIN+1))
         call DRC3JM(J1+1.0d0,dble(in),J1,-JM(ii,2),M2MIN, M2MAX, THRCOF, NDIM, IER3)
         if(IER3.eq.2) THRCOF=0.0d0
         abc=abc*THRCOF(int(-M2MIN+1))
         call DRC3JM(J1+1.0d0,1.0d0,J1+1.0d0,-JM(ii,2),M2MIN, M2MAX, THRCOF, NDIM, IER3)
         if(IER3.eq.2) THRCOF=0.0d0
         ghi=abc*rdscn(in,ir,1)*THRCOF(int(-M2MIN+1))**2
         ! 2
         call DRC3JM(J1+1.0d0,dble(in),J1+2.0d0,-JM(ii,2),M2MIN, M2MAX, THRCOF, NDIM, IER3)
         if(IER3.eq.2) THRCOF=0.0d0
         abc=THRCOF(int(-M2MIN+1))
         call DRC3JM(J1+2.0d0,1.0d0,J1+1.0d0      ,-JM(ii,2),M2MIN, M2MAX, THRCOF, NDIM, IER3)
         if(IER3.eq.2) THRCOF=0.0d0
         abc=abc*THRCOF(int(-M2MIN+1))
         call DRC3JM(J1,1.0d0,J1+1.0d0,-JM(ii,2),M2MIN, M2MAX, THRCOF, NDIM, IER3)
         if(IER3.eq.2) THRCOF=0.0d0
         ghi=ghi+abc*rdscn(in,ir,2)*THRCOF(int(-M2MIN+1))**2
         ! 3
         call DRC3JM(J1,dble(in),J1+1.0d0,-JM(ii,2),M2MIN, M2MAX, THRCOF, NDIM, IER3)
         if(IER3.eq.2) THRCOF=0.0d0
         abc=THRCOF(int(-M2MIN+1))
         call DRC3JM(J1+1.0d0,1.0d0,J1,-JM(ii,2),M2MIN, M2MAX, THRCOF, NDIM, IER3)
         if(IER3.eq.2) THRCOF=0.0d0
         abc=abc*THRCOF(int(-M2MIN+1))
         call DRC3JM(J1,1.0d0,J1+1.0d0,-JM(ii,2),M2MIN, M2MAX, THRCOF, NDIM, IER3)
         if(IER3.eq.2) THRCOF=0.0d0
         ghi=ghi+abc*rdscn(in,ir,1)*THRCOF(int(-M2MIN+1))**2

         if(temp.eq.0.0d0) then
            q3j(in,ir)=-ghi
            exit
         else
            q3j(in,ir)=q3j(in,ir)-ghi*dexp(-eigen(ii)/temp)
         end if
         
      end do ! do ii
      if(temp.ne.0.0d0) then
         q3j(in,ir)=q3j(in,ir)/dz
      end if

   end do ! do in
end do ! do ir




end subroutine multi3j


function d3jsym(J1,J2,J3,M1,M2)
integer(4),parameter :: NDIM=100
integer(4) :: IER3
real(8) :: M2MIN, M2MAX, THRCOF(NDIM),J1,J2,J3,M1,M2,d3jsym
call DRC3JM(J1,J2,J3,M1,M2MIN, M2MAX, THRCOF, NDIM, IER3)
if(IER3.eq.2) THRCOF=0.0d0
d3jsym=THRCOF(int(M2-M2MIN)+1)
return
end function d3jsym
