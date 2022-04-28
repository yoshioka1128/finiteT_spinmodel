subroutine multi3j2(L1,S1,J1,ZUst,imax,stfct,nJex,n4f,rdscn,HexT,nr,nrmax,lambda,d3j,d3j2,t3j,q3j,eigen,EJM,g,JM,temp)
implicit none
integer(4),parameter :: NDIM=100
real(8),parameter :: pi=dacos(-1.0d0)
complex(kind(0d0)),parameter :: im=(0.0d0,1.0d0)
integer(4) :: k,i,j, ii, in,ite, ll, mm, IER3, IER6, imax, i1, i2, l, m, jj1, jj2, imaxj,icount
integer(4) :: nJex,mm1,n4f,nrmax,nr,ir,nJ1,nJ2,nJ3,i1st,i2nd
real(8) :: L1, L2, L3, L4, L5, L6, M1, M2, M2MIN, M2MAX, THRCOF(NDIM),SIXCOF(NDIM), dlkl(3)
real(8) :: J1, J2, J3, dm, S1, Jin, abc,def,ghi,jkl,mno,lambda,dz,delta(2),delta12
real(8) :: L1MIN, L1MAX, gamma, stfct(6,-6:6),rdscn(6,nrmax,5),HexT(nrmax),d3j(6,nrmax),d3j2(6,nrmax),t3j(6,nrmax)
real(8) :: q3j(6,nrmax),w3jsym,w6jsym,EJM(imax,nrmax)
real(8) :: cclm(-6:6,-6:6),coff(6,-6:6),factorial(2),Ust(imax,imax,6,-6:6)
real(8) :: eigen(imax),g(100),JM(100,2),temp
real(8) :: rdlcl(0:6),rdjcj(0:6,2),rdjcjp(0:6,2),rdjcj2p(0:6,2),rdsss,rdjsj(2),rdjsjp(2)
!real(8) :: jcj(100,2),jcjp(100,2),jsj(100,2),jsjp(100,2)
complex(kind(0d0)) :: ZUst(imax,imax,6,-6:6)

! lz of signel particle
if(n4f.lt.7) mm1=4-n4f
if(n4f.gt.7) mm1=11-n4f

nJ1=int(J1*2.0d0)+1
nJ2=int((J1+1.0d0)*2.0d0)+1
nJ3=int((J1+2.0d0)*2.0d0)+1

i1st=nJ1+1
i2nd=nJ1+nJ2+1

delta(1)=lambda*(J1+1.0d0)
delta(2)=lambda*(J1+2.0d0)
delta12=lambda*0.5*((J1+2.0d0)*(J1+3.0d0)-(J1+1.0d0)*J1)

ir=1
delta(1)=EJM(i1st+nJ2/2,ir)-EJM(nJ1/2,ir)
delta(2)=EJM(i2nd+nJ3/2,ir)-EJM(i1st+nJ2/2,ir)
delta12=EJM(i2nd+nJ3/2,ir)-EJM(nJ1/2,ir)


! reduced matrix element
do k=2,6,2
   abc=0.0d0
   do m=mm1,3 ! lz of single particle
      abc=abc+(-1.0d0)**m*7.0d0*w3jsym(3.0d0,dble(k),3.0d0,0.0d0,0.0d0)*&
           w3jsym(3.0d0,dble(k),3.0d0,-dble(m),0.0d0)
   end do
   rdlcl(k)=abc/w3jsym(L1,dble(k),L1,-L1,0.0d0)
end do ! for k
rdsss=dsqrt(S1*(1.0d0+S1)*(1.0d0+2.0d0*S1)) ! <L||Cn||L>

do i=1,2 
   Jin=J1+dble(i-1)
   rdjsj(i)=(-1.0d0)**(int(L1+S1+Jin)+1)*(2.0d0*Jin+1.0d0)*&
        w6jsym(S1,Jin,L1,Jin,S1,1.0d0)*rdsss   ! <J||Sz||J>
   rdjsjp(i)=(-1.0d0)**(int(L1+S1+Jin))*dsqrt((2.0d0*Jin+1.0d0)*(2.0d0*Jin+3.0d0))*&
        w6jsym(S1,Jin,L1,Jin+1.0d0,S1,1.0d0)*rdsss   ! <J||Sz||J+1>
   do in=2,6,2
      rdjcj(in,i)=(-1.0d0)**(int(L1+S1+Jin)+in)*(2.0d0*Jin+1.0d0)*&
           w6jsym(L1,Jin,S1,Jin,L1,dble(in))*rdlcl(in)      ! <J||Cn||J>
      rdjcjp(in,i)=(-1.0d0)**(int(L1+S1+Jin)+in)*dsqrt((2.0d0*Jin+1.0d0)*(2.0d0*Jin+3.0d0))*&
           w6jsym(L1,Jin,S1,Jin+1.0d0,L1,dble(in))*rdlcl(in)      ! <J||Cn||J+1>
      rdjcj2p(in,i)=(-1.0d0)**(int(L1+S1+Jin)+in)*dsqrt((2.0d0*Jin+1.0d0)*(2.0d0*Jin+5.0d0))*&
           w6jsym(L1,Jin,S1,Jin+2.0d0,L1,dble(in))*rdlcl(in)      ! <J||Cn||J+2>
   end do
end do

do in=-2,2,1
   do ii=1,nJ1
!      write(6,*) J1,2.0d0,J1
!      write(6,*) -JM(ii,2),dble(in),JM(ii,2)-dble(in)
      if(abs(JM(ii,2)-dble(in)).gt.J1) write(6,*) w3jsym(J1,2.0d0,J1,-JM(ii,2),dble(in))
   end do
end do
stop
!write(6,*) w3jsym(J1,2.0d0,J1,J1,1.0d0)
write(6,*) J1,2.0d0,J1
write(6,*) J1,1.0d0,-J1-1.0d0
call DRC3JM(J1,2.0d0,J1,J1,M2MIN, M2MAX, THRCOF, NDIM, IER3)
!if(IER3.eq.2) THRCOF=0.0d0
abc=w3jsym(5.0d0/2.0d0,2.0d0,5.0d0/2.0d0,5.0d0/2.0d0,1.0d0)
!write(6,*) THRCOF(int(1.0d0-M2MIN)+1),abc

stop

! nonlinear term
d3j2=0.0d0
do ir=1,nr
   do in=2,6,2
      abc=4.0d0*HexT(ir)*rdjsjp(1)*(-rdjcjp(in,1))
      dz=0.0d0
      do ii=1,nJ1
         dz=dz+dexp(-EJM(ii,ir)/temp)
         def=w3jsym(J1,1.0d0,J1+1.0d0,-JM(ii,2),0.0d0)*w3jsym(J1+1.0d0,dble(in),J1,-JM(ii,2),0.0d0)
         if(temp.eq.0.0d0) then
            d3j(in,ir)=abc*def/(EJM(i1st+ii,ir)-EJM(ii,ir))
            write(6,*) "2nd",d3j(in,ir)
            exit
         else
            d3j(in,ir)=d3j(in,ir)+abc*def*dexp(-EJM(ii,ir)/temp)/(EJM(i1st+ii,ir)-EJM(ii,ir))
         end if
      end do ! do ii
      if(temp.ne.0.0d0) d3j(in,ir)=d3j(in,ir)/dz
   end do ! do in
end do ! do ir

! 2nd order
d3j=0.0d0
do ir=1,nr
   do in=2,6,2
      abc=4.0d0*HexT(ir)*rdjsjp(1)*(-rdjcjp(in,1))
      dz=0.0d0
      do ii=1,nJ1
         dz=dz+dexp(-EJM(ii,ir)/temp)
         def=w3jsym(J1,1.0d0,J1+1.0d0,-JM(ii,2),0.0d0)*w3jsym(J1+1.0d0,dble(in),J1,-JM(ii,2),0.0d0)
         if(temp.eq.0.0d0) then
            d3j(in,ir)=abc*def/(EJM(i1st+ii,ir)-EJM(ii,ir))
            write(6,*) "2nd",d3j(in,ir)
            exit
         else
            d3j(in,ir)=d3j(in,ir)+abc*def*dexp(-EJM(ii,ir)/temp)/(EJM(i1st+ii,ir)-EJM(ii,ir))
         end if
      end do ! do ii
      if(temp.ne.0.0d0) d3j(in,ir)=d3j(in,ir)/dz
   end do ! do in
end do ! do ir
!stop

! 3rd order
rdscn=0.0d0
t3j=0.0d0
do ir=1,nr
   do in=2,6,2
      jkl=-(2.0d0*HexT(ir))**2*2.0d0*(rdjsjp(1)*rdjcjp(in,1)) ! 1-2
      mno=-(2.0d0*HexT(ir))**2*(rdjsjp(1)**2) ! 3-4
      rdscn(in,ir,1)=0.0d0
      rdscn(in,ir,2)=0.0d0

      rdscn(in,ir,3)=mno*rdjcj(in,1)
      rdscn(in,ir,4)=mno*rdjcj(in,2)
      rdscn(in,ir,5)=-2.0d0*(2.0d0*HexT(ir))**2/3.0d0*&
           rdjsjp(1)*rdjsjp(2)*rdjcj2p(in,2) ! 5
      dz=0.0d0
      do ii=1,nJ1
         dz=dz+dexp(-EJM(ii,ir)/temp)
         ! 3-4
         def=w3jsym(J1,1.0d0,J1+1.0d0,-JM(ii,2),0.0d0)**2
         ghi=def*rdscn(in,ir,3)*w3jsym(J1,dble(in),J1,-JM(ii,2),0.0d0)/(EJM(i1st+ii,ir)-EJM(ii,ir))**2 ! 3
         ghi=ghi+def*rdscn(in,ir,4)*w3jsym(J1+1.0d0,dble(in),J1+1.0d0,-JM(ii,2),0.0d0)/(EJM(i1st+ii,ir)-EJM(ii,ir))**2 ! 4
         ! 5
         ghi=ghi+rdscn(in,ir,5)*w3jsym(J1,1.0d0,J1+1.0d0,-JM(ii,2),0.0d0)*&
              w3jsym(J1+1.0d0,1.0d0,J1+2.0d0,-JM(ii,2),0.0d0)*w3jsym(J1,dble(in),J1+2.0d0,-JM(ii,2),0.0d0)*&
              (-1.0d0/(delta12*delta(1)))
              
         if(temp.eq.0.0d0) then
            t3j(in,ir)=(-1.0d0)**(int(J1-JM(ii,2)))*ghi
            write(6,*) "3rd",d3j(in,ir)
            exit
         else
            t3j(in,ir)=t3j(in,ir)+(-1.0d0)**(int(J1-JM(ii,2)))*ghi*dexp(-EJM(ii,ir)/temp)
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
   do in=2,6,2
      jkl=1.0d0/2.0d0*(2.0d0*HexT(ir))**3/delta(2)**2/delta(1) 
      jkl=jkl-3.0d0/4.0d0*(2.0d0*HexT(ir))**3/delta(1)**2/delta(2) !1-2
      mno=-1.0d0/2.0d0*(2.0d0*HexT(ir))**3 !3-4
      rdscn(in,ir,1)=jkl*rdjcjp(in,1)*rdjsjp(1)*rdjsjp(2)**2
      rdscn(in,ir,2)=jkl*rdjcjp(in,2)*rdjsjp(2)*rdjsjp(1)**2
      rdscn(in,ir,3)=mno*rdjcjp(in,1)*rdjsjp(1)**3
      dz=0.0d0
      do ii=1,nJ1
         dz=dz+dexp(-EJM(ii,ir)/temp)
!         ghi=rdscn(in,ir,1)/(EJM(i1st+ii,ir)-EJM(ii,ir))/(EJM(i2nd+ii,ir)-EJM(i1st+ii,ir))**2*&         ! 1
         ghi=rdscn(in,ir,1)*&
              w3jsym(J1,1.0d0,J1+1.0d0,-JM(ii,2),0.0d0)*&
              w3jsym(J1,dble(in),J1+1.0d0,-JM(ii,2),0.0d0)*&
              w3jsym(J1+1.0d0,1.0d0,J1+2.0d0,-JM(ii,2),0.0d0)**2
         ghi=ghi+rdscn(in,ir,2)*&         ! 2
              w3jsym(J1+1.0d0,dble(in),J1+2.0d0,-JM(ii,2),0.0d0)*&
              w3jsym(J1,1.0d0,J1+1.0d0,-JM(ii,2),0.0d0)**2*&
              w3jsym(J1+1.0d0,1.0d0,J1+2.0d0,-JM(ii,2),0.0d0)
         ghi=ghi+rdscn(in,ir,3)/(EJM(i1st+ii,ir)-EJM(ii,ir))**3*&         ! 3
              w3jsym(J1,dble(in),J1+1.0d0,-JM(ii,2),0.0d0)*&
              w3jsym(J1,1.0d0,J1+1.0d0,-JM(ii,2),0.0d0)**3
         if(temp.eq.0.0d0) then
            q3j(in,ir)=ghi
            write(6,*) "4th",d3j(in,ir)
            exit
         else
            q3j(in,ir)=q3j(in,ir)+ghi*dexp(-EJM(ii,ir)/temp)
         end if
      end do ! do ii
      if(temp.ne.0.0d0) then
         q3j(in,ir)=q3j(in,ir)/dz
      end if

   end do ! do in
end do ! do ir

end subroutine multi3j2

