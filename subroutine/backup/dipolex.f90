!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine dipolec(diprfe,nunit,nunitx,nunity,nunitz,ax,cx,nr,rr,nfe,rfe,nx,ny,nz)
implicit none
integer(4) :: nunit,nr,nfe,nunitx,nunity,nunitz,ncycle
integer(4) :: i,i2,ir,ir2,ife,ife2,ix,iy,iz,ix2,iy2,iz2,ixd,iyd,izd,nx,ny,nz
integer(4) :: ipt1,ipt2,idim,idim2,icycley
real(8) :: rfe(nfe,3),rr(nr,3),ax,cx,delta(3,3),dis,dr(3),dfe(3),drfe(3),dfer(3),abc(3)
real(8) :: diprfe(-nunitx:nunitx,-nunity:nunity,-nunitz:nunitz,nr,nfe,3,3)

abc(1)=ax
abc(2)=ax
abc(3)=cx

delta=0.0d0
do i=1,3
   delta(i,i)=1.0d0
end do

do izd=1-nunitz,nunitz-1
   do iyd=1-nunity,nunity-1
      do ixd=1-nunitx,nunitx-1
         do i2=1,nfe
            do i=1,nr

               drfe(1)=dble(ixd-(2*nunitx-1)*nx)*ax+rr(i,1)-rfe(i2,1)
               drfe(2)=dble(iyd-(2*nunity-1)*ny)*ax+rr(i,2)-rfe(i2,2)
               drfe(3)=dble(izd-(2*nunitz-1)*nz)*cx+rr(i,3)-rfe(i2,3)
               
               dis=0.0d0
               do idim=1,3
                  dis=dis+drfe(idim)**2
               end do
               dis=dsqrt(dis)
               
               do idim2=1,3
                  do idim=1,3
                     diprfe(ixd,iyd,izd,i,i2,idim,idim2)=diprfe(ixd,iyd,izd,i,i2,idim,idim2)+&
                          (delta(idim,idim2)-3.0d0*drfe(idim)*drfe(idim2)/dis**2)/dis**3
                  end do
               end do

            end do
         end do
      end do
   end do
end do

end subroutine dipolec
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

