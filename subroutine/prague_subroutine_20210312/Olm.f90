subroutine Olm(L1,S1,ZOst,imax,theta)
implicit none
integer(4),parameter :: NDIM=100
real(8),parameter :: pi=dacos(-1.0d0)
complex(kind(0d0)),parameter :: im=(0.0d0,1.0d0)
integer(4) :: i, ite, ll, mm,mm2, IER, imax, i1, i2, l, m
real(8) :: L1, L2, L3, L4, L5, L6, M1, M2, M2MIN, M2MAX, THRCOF(NDIM),SIXCOF(NDIM)
real(8) :: J1, J2, dl, dm, S1, abc
real(8) :: L1MIN, L1MAX, gamma, theta(6)
real(8) :: clm(-6:6,-6:6),coff(6,-6:6),rdcoff(0:6),factorial(2)
complex(kind(0d0)) :: ZCst(100,100,6,-6:6),ZOst(imax,imax,6,-6:6)

theta=0.0d0
ZOst=0.0d0
do ll=2,6,2
write(6,*) "ll=",ll
!Nd
dl=dble(ll)
L1=6.0d0
S1=3.0d0/2.0d0
J1=abs(L1-S1)
imax=int(J1*2.0d0)+1
!Sm
dl=dble(ll)
L1=5.0d0
S1=5.0d0/2.0d0
J1=abs(L1-S1)
imax=int(J1*2.0d0)+1
! <||U||> table
rdcoff=0.0d0
rdcoff(2)=dsqrt(2.0d0*11.0d0*13.0d0/15.0d0)/3.0d0
rdcoff(4)=dsqrt(2.0d0*13.0d0/11.0d0)*2.0d0/3.0d0
rdcoff(6)=-10.0d0/3.0d0*dsqrt(5.0d0*17.0d0/(3.0d0*11.0d0*13.0d0))

! set alm table
clm=1.0d0
clm(2,0)=dsqrt(5.0d0/pi)/4.0d0
clm(2,2)=dsqrt(15.0d0/pi)/4.0d0
clm(4,0)=dsqrt(1.0d0/pi)*3.0d0/16.0d0
clm(4,2)=dsqrt(5.0d0/pi)*3.0d0/8.0d0
clm(4,4)=dsqrt(35.0d0/pi)*3.0d0/16.0d0
clm(6,0)=dsqrt(13.0d0/pi)/32.0d0
clm(6,2)=dsqrt(2730d0/pi)/64.0d0
clm(6,4)=dsqrt(13.0d0/7.0d0/pi)*21.0d0/32.0d0
clm(6,6)=dsqrt(26.0d0/231.0d0/pi)*231.0d0/64.0d0
do l=2,6,2
   do m=0,6,2
      clm(l,-m)=clm(l,m)
   end do
end do

coff=1.0d0
do l=2,6,2
   do m=-6,6,2
      if(m.eq.0) then
         coff(l,m)=1.0d0
      else
         coff(l,m)=1.0d0/dsqrt(2.0d0)
      end if
   end do
end do


ZCst=0.0d0
do ite=0,imax-1
   M1=-J1+dble(ite)
   call DRC3JM (J1, J1, dl, -M1, M2MIN, M2MAX, THRCOF, NDIM, IER)
!   if(IER.eq.2) write(6,*) "IER=2",THRCOF
   do i=1,int(M2MAX-M2MIN)+1
      M2=M2MIN+dble(i-1)
      mm=int(M1-M2)
      i1=int(M1+J1+1)
      i2=int(M2+J1+1)
      if(2.0d0*J1+dl+1.0d0.lt.0.0d0.or.2.0d0*J1-dl.lt.0.0d0) cycle
      factorial(1)=gamma(2.0d0*J1+dl+2.0d0)
      factorial(2)=gamma(2.0d0*J1-dl+1.0d0)
      ZCst(i1,i2,ll,mm)=coff(ll,mm)*1.0d0/clm(ll,mm)*dsqrt((2.0d0*dl+1.0d0)/(4.0d0*pi))/2.0d0**ll*&
           (-1.0d0)**(J1+M1)*dsqrt(factorial(1)/factorial(2))*THRCOF(i)
   end do
end do

do mm2=-ll,ll,2
   do i1=1,imax
      do i2=1,imax

         if(mm2.eq.0) then
            ZOst(i1,i2,ll,mm2)=ZCst(i1,i2,ll,mm2)
         else
            ZOst(i1,i2,ll,-abs(mm2))=im*(ZCst(i1,i2,ll,-abs(mm2))-(-1.0d0)**mm2*ZCst(i1,i2,ll,mm2))
            ZOst(i1,i2,ll,abs(mm2))=ZCst(i1,i2,ll,-abs(mm2))+(-1.0d0)**mm2*ZCst(i1,i2,ll,mm2)
         end if

      end do
   end do
end do


do mm2=-ll,ll,2
   do i1=1,imax
      do i2=1,i1
         
         if(abs(ZOst(i1,i2,ll,mm2)).ne.0.0d0) then
            write(6,*) ll,mm2,i1,i2
            write(6,*) ZOst(i1,i2,ll,mm2),ZOst(i2,i1,ll,mm2)
         end if

      end do
   end do
end do


call DRC6J(L1, dl, J1, J1, S1, L1MIN, L1MAX, SIXCOF, NDIM, IER)
!write(6,*) "six coff",L1MAX
do i=1,L1MAX-L1MIN+1
!   write(6,*) L1MIN+i-1,SIXCOF(i)
   if(L1MIN+i-1.eq.L1) abc=SIXCOF(i)
end do
i=L1-L1MIN+1
write(6,*) "Stevens factor"
factorial(1)=gamma(2.0d0*J1+dl+2.0d0)
factorial(2)=gamma(2.0d0*J1-dl+1.0d0)
theta(ll)=(-1.0d0)**(L1+S1+J1)*2.0d0**ll*(2.0d0*J1+1.0d0)*&
     dsqrt(factorial(2)/factorial(1))*abc*rdcoff(ll)
if(2.0d0*J1+dl+1.0d0.lt.0.0d0.or.2.0d0*J1-dl.lt.0.0d0) theta(ll)=0.0d0
end do

end subroutine Olm
