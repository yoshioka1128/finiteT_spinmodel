subroutine momentR_exp(temp,imax,zunitr,zspin,zorbit,partr,eigen0,stspinr,stspinSr,stspinLr)
implicit none
real(8),parameter :: pi=dacos(-1.0d0)
integer(4) :: ll,m,imax,jj,ii
real(8) :: stspinr(3),stspinSr(3),stspinLr(3),eigen0(imax),partr,temp
real(8) :: def,ghi
complex(kind(0d0)) :: zunitr(imax,imax),zspin(imax,imax,3),zorbit(imax,imax,3)

! angular momentum
stspinSr=0.0d0
stspinLr=0.0d0
stspinr=0.0d0
do ll=1,3
   do m=1,imax
      def=0.0d0
      ghi=0.0d0
      do jj=1,imax
         do ii=1,imax
            def=def-2.0d0*dreal(conjg(zunitr(ii,m))*zspin(ii,jj,ll)*zunitr(jj,m))
            ghi=ghi-dreal(conjg(zunitr(ii,m))*zorbit(ii,jj,ll)*zunitr(jj,m))
         end do
      end do
      if(temp.eq.0.0d0) then
         stspinSr(ll)=def
         stspinLr(ll)=ghi
         exit
      else
         stspinSr(ll)=stspinSr(ll)+def*dexp(-(eigen0(m)-eigen0(1))/temp)/partr
         stspinLr(ll)=stspinLr(ll)+ghi*dexp(-(eigen0(m)-eigen0(1))/temp)/partr
      end if
   end do ! m energy level
   stspinr(ll)=stspinSr(ll)+stspinLr(ll)
end do ! for spin component

end subroutine momentR_exp
