subroutine freeeng_divide_fromangle(imax,nr,free,exchloc,ealoc,theta,phi,zmagS,g,SO,zmoment,zmomentS,&
     zmatrix,ZUst,Clm,zmag2,HexT,Hmag,temp,eigen,part,zunit,RWORK,WORK)
implicit none
integer(4),parameter :: nrmax=30
integer(4) :: imax,ii,jj,ir,nr,m,LWORK,INFO,l
real(8) :: free(nrmax),part(nrmax),theta,phi,g(100),SO(100),Clm(6,0:6,nrmax),HexT(nrmax),Hmag,temp,eigen(imax),RWORK(3*imax-2)
real(8) :: exchloc(nrmax),ealoc(nrmax),abc
complex(kind(0d0)) :: zmagS(imax,imax),WORK(2*imax-1),zunit(imax,imax,nr),zmoment(imax,imax,3),zmomentS(imax,imax,3)
complex(kind(0d0)) :: zmatrix(imax,imax),ZUst(imax,imax,6,-6:6),zmag2(imax,imax),zabc,zdef

LWORK=2*imax-1
do jj=1,imax
   do ii=1,imax
      zmagS(ii,jj)=((g(ii)-1.0d0)*zmoment(ii,jj,3)+zmomentS(ii,jj,3))*dcos(theta)&
           +(((g(ii)-1.0d0)*zmoment(ii,jj,1)+zmomentS(ii,jj,1))*dcos(phi)&
           +((g(ii)-1.0d0)*zmoment(ii,jj,2)+zmomentS(ii,jj,2))*dsin(phi))*dsin(theta)
   end do
end do

free=0.0d0
part=0.0d0
do ir=1,nr ! sum up inequivalent R site
   zmatrix=0.0d0
   do jj=1,imax
      zmatrix(jj,jj)=SO(jj)
   end do
   
   do jj=1,imax
      do ii=jj,imax
         do l=2,6,2
            do m=0,l,2
               zmatrix(ii,jj)=zmatrix(ii,jj)+ZUst(ii,jj,l,m)*Clm(l,m,ir)
            end do
         end do
      end do
   end do
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
   do m=1,imax
      zabc=0.0d0
      zdef=0.0d0
      do jj=1,imax
         do ii=1,imax
            zabc=zabc+conjg(zunit(ii,m,ir))*2.0d0*zmagS(ii,jj)*HexT(ir)*zunit(jj,m,ir) ! exchange
            zdef=zdef+conjg(zunit(ii,m,ir))*ZUst(ii,jj,l,m)*Clm(l,m,ir)*zunit(jj,m,ir) ! anisotropy
         end do
      end do
      if(temp.eq.0.0d0) then
         exchloc(ir)=real(zabc)
         ealoc(ir)=real(zdef)
         exit
      else
         exchloc(ir)=exchloc(ir)+real(zabc)*dexp(-(eigen(ii)-eigen(1))/temp)/part(ir)
         ealoc(ir)=ealoc(ir)+real(zdef)*dexp(-(eigen(ii)-eigen(1))/temp)/part(ir)
      end if
   end do
end do ! for rare earth site ir

end subroutine freeeng_divide_fromangle
