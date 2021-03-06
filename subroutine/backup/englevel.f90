subroutine englevel(imax,nr,free,theta,phi,zmagS,g,SO,zmoment,zmomentS,&
     zmatrix,ZUst,Blm,zmag2,HexT,Hmag,temp,eigen,part,zunit,RWORK,WORK,ratom)
implicit none
integer(4) :: imax,ii,jj,ir,nr,m,LWORK,INFO,l
real(8) :: free(nr),part(nr),theta,phi,g(100),SO(100),Blm(6,0:6,4),HexT(nr),Hmag,temp,eigen(imax),RWORK(3*imax-2),eng0
complex(kind(0d0)) :: zmagS(imax,imax),WORK(2*imax-1),zunit(imax,imax,nr),zmoment(imax,imax,3),zmomentS(imax,imax,3)
complex(kind(0d0)) :: zmatrix(imax,imax),ZUst(imax,imax,6,-6:6),zmag2(imax,imax)
character(15) :: ratom

LWORK=2*imax-1
do jj=1,imax
   do ii=1,imax
      zmagS(ii,jj)=((g(ii)-1.0d0)*zmoment(ii,jj,3)+zmomentS(ii,jj,3))*dcos(theta)&
           +(((g(ii)-1.0d0)*zmoment(ii,jj,1)+zmomentS(ii,jj,1))*dcos(phi)&
           +((g(ii)-1.0d0)*zmoment(ii,jj,2)+zmomentS(ii,jj,2))*dsin(phi))*dsin(theta)
   end do
end do
open(17,file="englevel_"//trim(adjustl(ratom))//".txt")

free=0.0d0
part=0.0d0
do ir=1,nr ! sum up inequivalent R site


! set zero eng
   zmatrix=0.0d0
   do jj=1,imax
      zmatrix(jj,jj)=SO(jj)
   end do
   
   do jj=1,imax
      do ii=jj,imax
         do l=2,6,2
            do m=0,l,2
!               zmatrix(ii,jj)=zmatrix(ii,jj)+ZUst(ii,jj,l,m)*Blm(l,m,ir)
            end do
         end do
      end do
   end do
   do jj=1,imax
      do ii=1,imax
         zmatrix(ii,jj)=zmatrix(ii,jj)
         zmatrix(ii,jj)=zmatrix(ii,jj)+zmag2(ii,jj)*Hmag 
      end do
   end do
   LWORK=2*imax-1
   call zheev('V','L',imax,zmatrix,imax,eigen,WORK,LWORK,RWORK,INFO)
   free(ir)=eigen(1) 
   eng0=eigen(1)




! LS coupling only
   zmatrix=0.0d0
   do jj=1,imax
      zmatrix(jj,jj)=SO(jj)
   end do
   
   do jj=1,imax
      do ii=jj,imax
         do l=2,6,2
            do m=0,l,2
!               zmatrix(ii,jj)=zmatrix(ii,jj)+ZUst(ii,jj,l,m)*Blm(l,m,ir)
            end do
         end do
      end do
   end do
   do jj=1,imax
      do ii=1,imax
         zmatrix(ii,jj)=zmatrix(ii,jj)
         zmatrix(ii,jj)=zmatrix(ii,jj)+zmag2(ii,jj)*Hmag 
      end do
   end do
   LWORK=2*imax-1
   call zheev('V','L',imax,zmatrix,imax,eigen,WORK,LWORK,RWORK,INFO)
   free(ir)=eigen(1) 
   do ii=1,imax
      write(17,*) 1,eigen(ii)-eng0,ii
      write(17,*) 2,eigen(ii)-eng0,ii
      write(17,*)
   end do



! LS+2SHm
   zmatrix=0.0d0
   do jj=1,imax
      zmatrix(jj,jj)=SO(jj)
   end do

   do jj=1,imax
      do ii=jj,imax
         do l=2,6,2
            do m=0,l,2
               zmatrix(ii,jj)=zmatrix(ii,jj) !+ZUst(ii,jj,l,m)*Blm(l,m,ir)
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
   free(ir)=eigen(1) 
   do ii=1,imax
      write(17,*) 3,eigen(ii)-eng0,ii
      write(17,*) 4,eigen(ii)-eng0,ii
      write(17,*)
   end do



   zmatrix=0.0d0
   do jj=1,imax
      zmatrix(jj,jj)=SO(jj)
   end do
   
   do jj=1,imax
      do ii=jj,imax
         do l=2,6,2
            do m=0,l,2
               zmatrix(ii,jj)=zmatrix(ii,jj)+ZUst(ii,jj,l,m)*Blm(l,m,ir)
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
   free(ir)=eigen(1) 
   do ii=1,imax
      write(17,*) 5,eigen(ii)-eng0,ii
      write(17,*) 6,eigen(ii)-eng0,ii
      write(17,*)
   end do


   stop

   do jj=1,imax
      do ii=1,imax
         zunit(ii,jj,ir)=zmatrix(ii,jj)
      end do
   end do

end do ! for rare earth site ir

end subroutine englevel
