subroutine dispersion(itemp,imax,nr,g1,disp,HexT,temp,eigen,J1,JM)
implicit none
real(8),parameter :: pi=dacos(-1.0d0)
integer(4),parameter :: NDIM=100,nrmax=30,itempmax=1000
integer(4) :: imax,i,ii,jj,ir,nr,m,l,nJ1,in,itemp,im
real(8) :: part(nrmax),g1,rdmatJ(7)
real(8) :: HexT(nrmax,0:itempmax),temp,eigen(imax,nrmax),w3jsym
real(8) :: J1,JM(100,2),abc,def,disp(7,nrmax),Jz

do ir=1,nr
   do ii=1,imax ! energy
      eigen(ii,ir)=2.0d0*(g1-1.0d0)*HexT(ir,itemp)*JM(ii,2) ! +2S*HexT  HexT <0, g(1)-1 <0, <S> > 0
   end do
end do

! reduced matrix element
rdmatJ=0.0d0
do in=1,6
   abc=2.0d0*J1+dble(in)+1.0d0
   do jj=1,2*in 
      abc=abc*(2.0d0*J1+dble(in)+1.0d0-dble(jj))
   end do
   rdmatJ(in)=dsqrt(abc)/(2.0d0**in)
end do

disp=0.0d0
! number of M
nJ1=int(J1*2.0d0)+1 
if(temp.ne.0.0d0) then
   part=0.0d0
   do ir=1,nr ! inequivalent R site
      do ii=1,nJ1 ! partition function
         part(ir)=part(ir)+dexp(-eigen(ii,ir)/temp)
      end do
      
      do in=2,2 ! rank
         do im=0,0
            do ii=1,nJ1 ! Jz=-J1+(ii-1)
               def=(-1.0d0)**(J1-JM(ii,2))*rdmatJ(in)*w3jsym(J1,dble(in),J1,-JM(ii,2),dble(im)) ! JM2=JM-im
               disp(in,ir)=disp(in,ir)+def**2*dexp(-eigen(ii,ir)/temp)
            end do ! do ii, M
            disp(in,ir)=disp(in,ir)/part(ir)
!            write(6,*) in,im,ir,CC(in,im,ir)
         end do
      end do ! do in, rank
   end do ! for rare earth site ir
!   write(6,*) "disp check finite-T"
!   write(6,*) disp(2,1),((3.0d0*J1**2-J1*(J1+1.0d0))/2.0d0)**2
end if


end subroutine dispersion
