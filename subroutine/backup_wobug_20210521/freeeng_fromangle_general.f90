subroutine freeeng_fromangle_general(imax,nr,free,theta,phi,zmagS,g,SO,zmoment,zmomentS,&
     zmatrix,ZUst,Clm,zmag2,HexT,Hmag,temp,eigen,part,zunit,RWORK,WORK)
implicit none
integer(4),parameter :: nrmax=30
integer(4) :: imax,ii,jj,ir,nr,m,LWORK,INFO,l
real(8) :: free(nrmax),part(nrmax),theta,phi,g(100),SO(100),Clm(6,-6:6,nrmax),HexT(nrmax),Hmag,temp,eigen(imax),RWORK(3*imax-2)
complex(kind(0d0)) :: zmagS(imax,imax),WORK(2*imax-1),zunit(imax,imax,nr),zmoment(imax,imax,3),zmomentS(imax,imax,3)
complex(kind(0d0)) :: zmatrix(imax,imax),ZUst(imax,imax,6,-6:6),zmag2(imax,imax)

LWORK=2*imax-1
do jj=1,imax
   do ii=1,imax
      zmagS(ii,jj)=&
           ((g(ii)-1.0d0)*zmoment(ii,jj,3)+zmomentS(ii,jj,3))*dcos(theta)&
           +(((g(ii)-1.0d0)*zmoment(ii,jj,1)+zmomentS(ii,jj,1))*dcos(phi)&
           +((g(ii)-1.0d0)*zmoment(ii,jj,2)+zmomentS(ii,jj,2))*dsin(phi))*dsin(theta)
   end do
end do

free=0.0d0
part=0.0d0
do ir=1,nr ! sum up inequivalent R site
   zmatrix=0.0d0
   ! SOI
   do jj=1,imax
      zmatrix(jj,jj)=SO(jj)
   end do
   ! MA
   do jj=1,imax
      do ii=jj,imax
         do l=2,6,2
            do m=-l,l
               zmatrix(ii,jj)=zmatrix(ii,jj)+ZUst(ii,jj,l,m)*Clm(l,m,ir)
!               if(ir.eq.1.and.l.eq.2.and.abs(m).eq.1.and.abs(ZUst(ii,jj,l,m)).ne.0.0d0) write(6,*) m,ii,jj,ZUst(ii,jj,l,m)*Clm(l,m,ir)
            end do
         end do
      end do
   end do
!   stop
   ! Zeeman
   do jj=1,imax
      do ii=1,imax
         zmatrix(ii,jj)=zmatrix(ii,jj)+2.0d0*zmagS(ii,jj)*HexT(ir) 
         zmatrix(ii,jj)=zmatrix(ii,jj)+zmag2(ii,jj)*Hmag 
      end do
   end do

   LWORK=2*imax-1
   call zheev('V','L',imax,zmatrix,imax,eigen,WORK,LWORK,RWORK,INFO)
   if(temp.eq.0.0d0) then
      free(ir)=eigen(1) 
   else
      do  ii=1,imax
         part(ir)=part(ir)+dexp(-(eigen(ii)-eigen(1))/temp)
      end do
      free(ir)=eigen(1)-temp*dlog(part(ir))
!      free(ir)=free(ir)+eigen(1)-temp*dlog(part(ir))
   end if

   do jj=1,imax
      do ii=1,imax
         zunit(ii,jj,ir)=zmatrix(ii,jj)
      end do
   end do

end do ! for rare earth site ir

end subroutine freeeng_fromangle_general
