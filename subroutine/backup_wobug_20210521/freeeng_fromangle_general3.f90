subroutine freeeng_fromangle_general(imax,nr,free,zH,&
     temp,eigen,part,RWORK,WORK)
implicit none
integer(4),parameter :: nrmax=30
integer(4) :: imax,ii,jj,ir,nr,m,LWORK,INFO,l
real(8) :: free(nrmax),part(nrmax),temp,eigen(imax),RWORK(3*imax-2)
complex(kind(0d0)) :: WORK(2*imax-1),zH(imax,imax,nr)
complex(kind(0d0)) :: zH2(imax,imax,nrmax)

free=0.0d0
part=0.0d0
do ir=1,nr
   LWORK=2*imax-1
   call zheev('V','L',imax,zH(1,1,ir),imax,eigen,WORK,LWORK,RWORK,INFO)
!   call zheev('V','L',imax,zmatrix,imax,eigen,WORK,LWORK,RWORK,INFO)
   if(temp.eq.0.0d0) then
      free(ir)=eigen(1) 
   else
      do  ii=1,imax
         part(ir)=part(ir)+dexp(-(eigen(ii)-eigen(1))/temp)
      end do
      free(ir)=eigen(1)-temp*dlog(part(ir))
   end if

end do ! for rare earth site ir

end subroutine freeeng_fromangle_general
