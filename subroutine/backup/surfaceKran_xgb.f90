subroutine surfaceKran_xgb(Klcr,Ku,ran,percent,irpt,nr,nunit,nunitx,nunity,nunitz,icount,irvx)
integer(4) :: nunit,nunitx,nunity,nunitz,nr
integer(4) :: irvx,irpt(nr,nunitx,nunity,nunitz),ix,iy,iz
integer,allocatable :: seed(:)
integer(4) :: i1,i2,i3,i4,icount,seedsize,ic
real(8) :: Klcr(9,nr*nunit),Ku,ran(nunity*nunitz*4),percent

! random
call system_clock(count=ic)
call random_seed(size=seedsize)
allocate(seed(seedsize))
call random_seed(get=seed)
seed=492983
call random_seed(put=seed)
call random_number(ran)

if(percent.le.1.0d0) then
   icount=0
   do iz=1,nunitz
      do iy=1,nunity
         i1=iy+(iz-1)*nunity
         i2=iy+(iz-1)*nunity+nunitz*nunity
         i3=iy+(iz-1)*nunity+2*nunitz*nunity
         i4=iy+(iz-1)*nunity+3*nunitz*nunity
         if(ran(i1).lt.percent) then
            icount=icount+1
            Klcr(1,irpt(4,irvx,iy,iz))=-Ku
            write(6,*) i1,ran(i1)
         end if
         if(ran(i2).lt.percent) then
            icount=icount+1
            Klcr(1,irpt(5,irvx,iy,iz))=-Ku
            write(6,*) i2,ran(i2)
         end if
         if(ran(i3).lt.percent) then
            icount=icount+1
            Klcr(1,irpt(2,irvx,iy,iz))=-Ku
            write(6,*) i3,ran(i3)
         end if
         if(ran(i4).lt.percent) then
            icount=icount+1
            Klcr(1,irpt(7,irvx,iy,iz))=-Ku
            write(6,*) i4,ran(i4)
         end if
      end do
   end do
else
   do iz=1,nunitz
      do iy=1,nunity
         Klcr(1,irpt(4,irvx,iy,iz))=-Ku
         Klcr(1,irpt(5,irvx,iy,iz))=-Ku
         Klcr(1,irpt(2,irvx,iy,iz))=-Ku
         Klcr(1,irpt(7,irvx,iy,iz))=-Ku
      end do
   end do
   icount=nunity*nunitz*4
   do iz=1,nunitz
      do iy=1,nunity
         i1=iy+(iz-1)*nunity
         i2=iy+(iz-1)*nunity+nunitz*nunity
         i3=iy+(iz-1)*nunity+2*nunitz*nunity
         i4=iy+(iz-1)*nunity+3*nunitz*nunity
         if(ran(i1).lt.percent-1) then
            icount=icount+1
            Klcr(1,irpt(1,irvx,iy,iz))=-Ku
            write(6,*) i1,ran(i1)
         end if
         if(ran(i2).lt.percent-1) then
            icount=icount+1
            Klcr(1,irpt(3,irvx,iy,iz))=-Ku
            write(6,*) i2,ran(i2)
         end if
         if(ran(i3).lt.percent-1) then
            icount=icount+1
            Klcr(1,irpt(6,irvx,iy,iz))=-Ku
            write(6,*) i3,ran(i3)
         end if
         if(ran(i4).lt.percent-1) then
            icount=icount+1
            Klcr(1,irpt(8,irvx,iy,iz))=-Ku
            write(6,*) i4,ran(i4)
         end if
      end do
   end do
end if
write(6,*) "adsorbed oxigen"
write(6,*) dble(icount)/dble(nunity*nunitz*4)*100.0d0,"%"

end subroutine surfaceKran_xgb
