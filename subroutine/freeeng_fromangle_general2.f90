subroutine freeeng_fromangle_general2(imax,nr,freeloc,theta,phi,zmagS,g,SO,zmoment,zmomentS,&
     zmatrix,ZUst,Clm,zmag2,HexT,Hmag,temp,eigen,part,zunit,RWORK,WORK,ir,zH0)
implicit none
integer(4),parameter :: nrmax=30
integer(4) :: imax,ii,jj,ir,nr,m,LWORK,INFO,l
real(8) :: freeloc(nrmax),part(nrmax),theta,phi,g(100),SO(100),Clm(6,-6:6,nrmax),HexT(nrmax),Hmag,temp,eigen(imax),RWORK(3*imax-2)
complex(kind(0d0)) :: zmagS(imax,imax),WORK(2*imax-1),zunit(imax,imax,nr),zmoment(imax,imax,3),zmomentS(imax,imax,3)
complex(kind(0d0)) :: zmatrix(imax,imax),ZUst(imax,imax,6,-6:6),zmag2(imax,imax),zH0(imax,imax,nr)


   LWORK=2*imax-1
   call zheev('V','L',imax,zmatrix,imax,eigen,WORK,LWORK,RWORK,INFO)
   if(temp.eq.0.0d0) then
      freeloc(ir)=eigen(1) 
   else
      do  ii=1,imax
         part(ir)=part(ir)+dexp(-(eigen(ii)-eigen(1))/temp)
      end do
      freeloc(ir)=eigen(1)-temp*dlog(part(ir))
!      freeloc(ir)=freeloc(ir)+eigen(1)-temp*dlog(part(ir))
   end if

   do jj=1,imax
      do ii=1,imax
         zunit(ii,jj,ir)=zmatrix(ii,jj)
      end do
   end do

!end do ! for rare earth site ir

end subroutine freeeng_fromangle_general2
