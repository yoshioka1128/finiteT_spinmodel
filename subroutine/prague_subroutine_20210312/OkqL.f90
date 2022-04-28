subroutine OkqL(ZUst,Ust,imax)
implicit none
integer(4),parameter :: nrmax=30
real(8) :: Ust(imax,imax,6,-6:6)
integer(4) :: l,m,i1,i2,imax
complex(kind(0d0)),parameter :: im=(0.0d0,1.0d0)
complex(kind(0d0)) :: ZUst(imax,imax,6,-6:6)


ZUst=(0.0d0,0.0d0)
do l=2,6,2
   do m=-l,l
      do i2=1,imax
         do i1=1,imax
            if(m.eq.0) then
               ZUst(i1,i2,l,m)=Ust(i1,i2,l,m)
            else if(m.lt.0) then
               ZUst(i1,i2,l,m)=im*(Ust(i1,i2,l,-abs(m))-(-1)**m*Ust(i1,i2,l,abs(m)))
            else if(m.gt.0) then
               ZUst(i1,i2,l,m)=Ust(i1,i2,l,-abs(m))+(-1)**m*Ust(i1,i2,l,abs(m))
            end if
         end do
      end do
   end do
end do

end subroutine OkqL
