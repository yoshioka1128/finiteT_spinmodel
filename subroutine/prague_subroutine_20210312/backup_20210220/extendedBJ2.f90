subroutine extendedBJ2(imax,extBJ,HexT,temp,J1,JM,S1)
implicit none
integer(4) :: imax,i,ii,jj,m,l,nJ1,in
real(8) :: part,HexT
real(8) :: temp,eigen(100),w3jsym
real(8) :: S1,J1,JM(100,2),abc,def,extBJ(-1:8)

eigen=0.0d0
do ii=1,imax ! energy
   eigen(ii)=-2.0d0*S1/(J1+1.0d0)*HexT*JM(ii,2) ! HexT <0, g(1)-1 <0, <S> > 0
end do

! number of M
extBJ=0.0d0
nJ1=int(J1*2.0d0)+1 
if(temp.eq.0.0d0) then
      do in=1,8 ! rank
         extBJ(in)=1.0d0
         do i=0,in-1
            extBJ(in)=extBJ(in)*(2.0d0*J1-dble(i))
         end do
         extBJ(in)=(-1.0d0)**in*extBJ(in)/(2.0d0**in) ! (-1)^in is from "odd commutation" or "(-1)^(J-M)"
      end do
else
   part=0.0d0
      do ii=1,nJ1 ! partition function
         part=part+dexp(-eigen(ii)/temp)
      end do
      
      do in=1,8 ! rank
         do ii=1,nJ1 ! M
            ! factorical
            abc=2.0d0*J1+dble(in)+1.0d0
            do jj=1,2*in 
               abc=abc*(2.0d0*J1+dble(in)+1.0d0-dble(jj))
            end do
            extBJ(in)=extBJ(in)+(-1.0d0)**(J1-JM(ii,2))*dsqrt(abc)/(2.0d0**in)*dexp(-eigen(ii)/temp)*&
                 w3jsym(J1,dble(in),J1,-JM(ii,2),0.0d0)
         end do ! do ii, M
         extBJ(in)=extBJ(in)/part
      end do ! do in, rank
end if

! positive define
   do in=1,7
      extBJ(in)=dabs(extBJ(in))
!      extBJ(in)=extBJ(in)*(-1)**in
   end do
   extBJ(0)=1.0d0
   extBJ(-1)=0.0d0

end subroutine extendedBJ2
