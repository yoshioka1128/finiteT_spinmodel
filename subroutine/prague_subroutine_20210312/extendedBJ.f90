subroutine extendedBJ(imax,nr,free,theta,phi,zmagS,g,SO,zmoment,zmomentS,extBJ,&
     zmatrix,ZUst,Blm,zmag2,HexT,Hmag,temp,eigen,part,zunit,RWORK,WORK,J1,JM)
implicit none
integer(4),parameter :: NDIM=100,nrmax=30
integer(4) :: imax,i,ii,jj,ir,nr,m,LWORK,INFO,l,nJ1,in,IER3
real(8) :: free(nrmax),part(nrmax),theta,phi,g(100),SO(100),Blm(6,0:6,nrmax)
real(8) :: HexT(nrmax),Hmag,temp,eigen(imax),RWORK(3*imax-2)
real(8) :: J1,JM(100,2),abc,def,extBJ(7,nrmax),M2MIN,M2MAX,M2, THRCOF(NDIM),Brillouin
complex(kind(0d0)) :: zmagS(imax,imax),WORK(2*imax-1),zunit(imax,imax,nr),zmoment(imax,imax,3),zmomentS(imax,imax,3)
complex(kind(0d0)) :: zmatrix(imax,imax),ZUst(imax,imax,6,-6:6),zmag2(imax,imax)

! number of M
nJ1=int(J1*2.0d0)+1 
if(temp.eq.0.0d0) then
   do ir=1,nr
      do in=1,7
         extBJ(in,ir)=1.0d0
         do i=0,in-1
            extBJ(in,ir)=extBJ(in,ir)*(2.0d0*J1-dble(i))
         end do
         extBJ(in,ir)=(-1.0d0)**in*extBJ(in,ir)/(2.0d0**in)
      end do
   end do
else
part=0.0d0
extBJ=0.0d0
do ir=1,nr ! inequivalent R site
   eigen=0.0d0
   do ii=1,imax ! energy
      eigen(ii)=2.0d0*(g(1)-1.0d0)*HexT(ir)*JM(ii,2)
   end do
   do ii=1,nJ1 ! partition function
      part(ir)=part(ir)+dexp(-eigen(ii)/temp)
   end do
   
   do in=1,7 ! rank
      do ii=1,nJ1 ! M
         ! factorical
         abc=2.0d0*J1+dble(in)+1.0d0
         do jj=1,2*in 
            abc=abc*(2.0d0*J1+dble(in)+1.0d0-dble(jj))
         end do
         ! (J1, in, J1)
         ! (-M, 0, M)
         call DRC3JM (J1, dble(in), J1, -JM(ii,2), M2MIN, M2MAX, THRCOF, NDIM, IER3)
         do i=1,int(M2MAX-M2MIN)+1
            M2=M2MIN+dble(i-1)
            if(M2.eq.0.0d0) def=THRCOF(i)
         end do
         extBJ(in,ir)=extBJ(in,ir)+(-1.0d0)**(J1-JM(ii,2))*dsqrt(abc)/(2.0d0**in)*def*dexp(-eigen(ii)/temp)
!         write(6,*) ir,in,(-1.0d0)**(J1-JM(ii,2))*dsqrt(abc)/(2.0d0**in)*def*dexp(-eigen(ii)/temp)
         !      write(6,*) in,ii,extBJ(in,ir)
!         if(ir.eq.1) write(6,*) in,ii,ir,extBJ(in,ir) 
      end do ! do ii, M
      extBJ(in,ir)=extBJ(in,ir)/part(ir)
   end do ! do in, rank
!   write(6,*) extBJ(1,ir)/J1,Brillouin(J1,2.0d0*(g(1)-1.0d0)*abs(HexT(ir))*J1/temp)
!   stop
end do ! for rare earth site ir
!write(6,*) extBJ(2,1),extBJ(2,2)
!stop ! delete temp=300.0d0
end if

end subroutine extendedBJ

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!function Brillouin(J,x)
!real(8) :: x,J,Brillouin
!
!Brillouin=(2.0d0*J+1.0d0)/(2.0d0*J)/dtanh((2.0d0*J+1.0d0)*x/(2.0d0*J))-1.0d0/(2.0d0*J)/dtanh(x/(2.0d0*J))
!
!end function Brillouin
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
