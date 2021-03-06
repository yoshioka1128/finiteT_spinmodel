subroutine surface_SmFe12(dK1,dK2,dK3,dK1surf,dK2surf,dK3surf,Ku,defectrate,nr,nunit,icount,irunit,irsub,nunitx,nunity,nunitz,nrtot,ranfix,iran)
integer(4) :: nunit,nr,irunit(nr*nunit,3),irsub(nr*nunit)
integer(4) :: ix,iy,iz,nrtot,ir0,c,iran
integer,allocatable :: seed(:)
integer(4) :: i1,i2,i3,i4,icount,seedsize,ic,nunitx,nunity,nunitz
real(8) :: Ku,defectrate
real(8) :: dK1(7,nr*nunit),dK2(7,nr*nunit),dK3(7,nr*nunit)
real(8) :: dK1surf(7,nr),dK2surf(7,nr),dK3surf(7,nr)
real(8),parameter :: mev=1.60218d-21 ! (Merg/mev) 
real(8),parameter :: kb=8.6173303d-2*mev ! (Merg/K)
real(8) :: x
character(10) :: ranfix

! random
call random_seed(size=seedsize)
allocate(seed(seedsize))
if(ranfix.eq."on") then
   call random_seed(put=(/iran/))
else
   call system_clock(count=c)
   call random_seed(put=(/c/))
end if

do ir=1,nrtot
   i=irsub(ir)
   ix=irunit(ir,1)
   iy=irunit(ir,2)
   iz=irunit(ir,3)
   if(ix.eq.1) then
      if(i.le.nr) then
         do ik=1,7
            dK1(ik,ir)=dK1surf(ik,i)
            dK2(ik,ir)=dK2surf(ik,i)
            dK3(ik,ir)=dK3surf(ik,i)
         end do
      end if
   end if
end do

open(93,file="R_defect_K1_in_gb_phase.txt")
icount=0
do ir=1,nrtot
   if(dK1(1,ir).lt.0.0d0) then
      write(93,*) irsub(ir),irunit(ir,1),irunit(ir,2),irunit(ir,3),dK1(1,ir)/kb
      icount=icount+1
   end if
end do

end subroutine surface_SmFe12


subroutine surface_Nd2Fe14B(dK1,dK2,dK3,dK1surf,dK2surf,dK3surf,Ku,defectrate,nr,nunit,icount,irunit,irsub,nunitx,nunity,nunitz,nrtot,ranfix,iran)
integer(4) :: nunit,nr,irunit(nr*nunit,3),irsub(nr*nunit)
integer(4) :: ix,iy,iz,nrtot,ir0,c,iran
integer,allocatable :: seed(:)
integer(4) :: i1,i2,i3,i4,icount,seedsize,ic,nunitx,nunity,nunitz
real(8) :: Ku,defectrate
real(8) :: dK1(7,nr*nunit),dK2(7,nr*nunit),dK3(7,nr*nunit)
real(8) :: dK1surf(7,nr),dK2surf(7,nr),dK3surf(7,nr)
real(8),parameter :: mev=1.60218d-21 ! (Merg/mev) 
real(8),parameter :: kb=8.6173303d-2*mev ! (Merg/K)
real(8) :: x
character(10) :: ranfix

! random
call random_seed(size=seedsize)
allocate(seed(seedsize))
if(ranfix.eq."on") then
   call random_seed(put=(/iran/))
else
   call system_clock(count=c)
   call random_seed(put=(/c/))
end if

do ir=1,nrtot
   i=irsub(ir)
   ix=irunit(ir,1)
   iy=irunit(ir,2)
   iz=irunit(ir,3)
   if(iz.eq.nunitz) then
      if(i.le.nr) then
         do ik=1,7
            dK1(ik,ir)=dK1surf(ik,i)
            dK2(ik,ir)=dK2surf(ik,i)
            dK3(ik,ir)=dK3surf(ik,i)
         end do
      end if
   end if
end do

open(93,file="R_defect_K1_in_gb_phase.txt")
icount=0
do ir=1,nrtot
   if(dK1(1,ir).lt.0.0d0) then
      write(93,*) irsub(ir),irunit(ir,1),irunit(ir,2),irunit(ir,3),dK1(1,ir)/kb
      icount=icount+1
   end if
end do

end subroutine surface_Nd2Fe14B




subroutine defectR_edge(Klcr,Ku,nr,nunit,icount,nunitx,nunity,nunitz,irunit,irsub,nrtot)
integer(4) :: nunit,nr,nunitx,nunity,nunitz
integer(4) :: ix,iy,iz,nrtot,ir0
integer,allocatable :: seed(:)
integer(4) :: i1,i2,i3,i4,icount,seedsize,ic,irunit(nr*nunit,3),irsub(nr*nunit)
real(8) :: Klcr(9,nr*nunit),Ku
real(8),parameter :: mev=1.60218d-21 ! (Merg/mev) 
real(8),parameter :: kb=8.6173303d-2*mev ! (Merg/K)

iz=nunitz
iy=1
ix=nunitx
do ir=1,nrtot
   if(irunit(ir,3).eq.iz) then
      if(irunit(ir,2).ge.iy-1.and.irunit(ir,2).le.iy) then
         if(irunit(ir,1).ge.ix-1.and.irunit(ir,1).le.ix) then
            Klcr(1,ir)=-Ku
            do ik=2,9
               Klcr(ik,ir)=0.0d0
            end do
         end if
      end if
   end if
end do

open(93,file="R_defect_K1_in_gb_phase.txt")
icount=0
do ir=1,nrtot
   if(Klcr(1,ir).lt.0.0d0) then
      write(93,*) irsub(ir),irunit(ir,1),irunit(ir,2),irunit(ir,3),Klcr(1,ir)/kb
      icount=icount+1
   end if
end do

end subroutine defectR_edge


subroutine defectR_clusterx4(Klcr,Ku,nr,nunit,nunitx,nunity,nunitz,icount,irunit,irsub,nrtot)
integer(4) :: nunit,nunitx,nunity,nunitz,nr
integer(4) :: irunit(nr*nunit,3),irsub(nr*nunit),nrtot,ir0
integer,allocatable :: seed(:)
integer(4) :: i1,i2,i3,i4,icount,seedsize,ic
real(8) :: Klcr(9,nr*nunit),Ku
real(8),parameter :: mev=1.60218d-21 ! (Merg/mev) 
real(8),parameter :: kb=8.6173303d-2*mev ! (Merg/K)
do ir=1,nrtot
   if(irunit(ir,3).ge.nunitz/2.and.irunit(ir,3).le.nunitz/2+1) then
      if(irunit(ir,2).ge.nunity/2.and.irunit(ir,2).le.nunity/2+1) then
         if(irunit(ir,1).eq.nunitx) then
            Klcr(1,ir)=-Ku
            do ik=2,9
               Klcr(ik,ir)=0.0d0
            end do
         end if
      end if
   end if
end do

open(93,file="R_defect_K1_in_gb_phase.txt")
icount=0
do ir=1,nrtot
   if(Klcr(1,ir).lt.0.0d0) then
      write(93,*) irsub(ir),irunit(ir,1),irunit(ir,2),irunit(ir,3),Klcr(1,ir)/kb
      icount=icount+1
   end if
end do

end subroutine defectR_clusterx4

subroutine defectR_clusterz4(Klcr,Ku,nr,nunit,nunitx,nunity,nunitz,icount,irunit,irsub,nrtot)
integer(4) :: nunit,nunitx,nunity,nunitz,nr
integer(4) :: irunit(nr*nunit,3),irsub(nr*nunit),nrtot,ir0
integer,allocatable :: seed(:)
integer(4) :: i1,i2,i3,i4,icount,seedsize,ic
real(8) :: Klcr(9,nr*nunit),Ku
real(8),parameter :: mev=1.60218d-21 ! (Merg/mev) 
real(8),parameter :: kb=8.6173303d-2*mev ! (Merg/K)

do ir=1,nrtot
   if(irunit(ir,1).ge.nunitx/2.and.irunit(ir,1).le.nunitx/2+1) then
      if(irunit(ir,2).ge.nunity/2.and.irunit(ir,2).le.nunity/2+1) then
         if(irunit(ir,3).eq.nunitz) then
            Klcr(1,ir)=-Ku
            do ik=2,9
               Klcr(ik,ir)=0.0d0
            end do
         end if
      end if
   end if
end do

open(93,file="R_defect_K1_in_gb_phase.txt")
icount=0
do ir=1,nrtot
   if(Klcr(1,ir).lt.0.0d0) then
      write(93,*) irsub(ir),irunit(ir,1),irunit(ir,2),irunit(ir,3),Klcr(1,ir)/kb
      icount=icount+1
   end if
end do

end subroutine defectR_clusterz4



subroutine defectR_randomx(Klcr,Ku,defectrate,nr,nunit,icount,irunit,irsub,nunitx,nunity,nunitz,nrtot,ranfix,iran)
integer(4) :: nunit,nr,irunit(nr*nunit,3),irsub(nr*nunit)
integer(4) :: ix,iy,iz,nrtot,ir0,c,iran
integer,allocatable :: seed(:)
integer(4) :: i1,i2,i3,i4,icount,seedsize,ic,nunitx,nunity,nunitz
real(8) :: Klcr(9,nr*nunit),Ku,defectrate
real(8),parameter :: mev=1.60218d-21 ! (Merg/mev) 
real(8),parameter :: kb=8.6173303d-2*mev ! (Merg/K)
real(8) :: x
character(10) :: ranfix

! random
call random_seed(size=seedsize)
allocate(seed(seedsize))
if(ranfix.eq."on") then
   call random_seed(put=(/iran/))
else
   call system_clock(count=c)
   call random_seed(put=(/c/))
end if

if(defectrate.le.50.0d0) then
   do ir=1,nrtot
      i=irsub(ir)
      ix=irunit(ir,1)
      iy=irunit(ir,2)
      iz=irunit(ir,3)
      if(ix.eq.nunitx) then
         if(i.eq.2.or.i.eq.4.or.i.eq.5.or.i.eq.7) then
            call random_number(x)
            if(x*50.0d0.lt.defectrate) then
               Klcr(1,ir)=-Ku
               do ik=2,9
                  Klcr(ik,ir)=0.0d0
               end do
            end if
         end if

      end if
   end do
else
   do ir=1,nrtot
      i=irsub(ir)
      ix=irunit(ir,1)
      iy=irunit(ir,2)
      iz=irunit(ir,3)
      if(ix.eq.nunitx) then
         if(i.eq.2.or.i.eq.4.or.i.eq.5.or.i.eq.7) then
            Klcr(1,ir)=-Ku
         else
            call random_number(x)
            if(x*50.0d0.lt.defectrate-50.0d0) then
               Klcr(1,ir)=-Ku
               do ik=2,9
                  Klcr(ik,ir)=0.0d0
               end do
            end if
         end if
      end if
   end do
end if

open(93,file="R_defect_K1_in_gb_phase.txt")
icount=0
do ir=1,nrtot
   if(Klcr(1,ir).lt.0.0d0) then
      write(93,*) irsub(ir),irunit(ir,1),irunit(ir,2),irunit(ir,3),Klcr(1,ir)/kb
      icount=icount+1
   end if
end do

end subroutine defectR_randomx


subroutine defectR_randomz(Klcr,Ku,defectrate,nr,nunit,icount,irunit,irsub,nunitx,nunity,nunitz,nrtot,ranfix,iran)
integer(4) :: nunit,nr,irunit(nr*nunit,3),irsub(nr*nunit)
integer(4) :: ix,iy,iz,nrtot,ir0,c,iran
integer,allocatable :: seed(:)
integer(4) :: i1,i2,i3,i4,icount,seedsize,ic,nunitx,nunity,nunitz
real(8) :: Klcr(9,nr*nunit),Ku,defectrate
real(8),parameter :: mev=1.60218d-21 ! (Merg/mev) 
real(8),parameter :: kb=8.6173303d-2*mev ! (Merg/K)
real(8) :: x
character(10) :: ranfix

! random
call random_seed(size=seedsize)
allocate(seed(seedsize))
if(ranfix.eq."on") then
   call random_seed(put=(/iran/))
else
   call system_clock(count=c)
   call random_seed(put=(/c/))
end if

if(defectrate.le.100.0d0) then
   do ir=1,nrtot
      i=irsub(ir)
      ix=irunit(ir,1)
      iy=irunit(ir,2)
      iz=irunit(ir,3)
      if(iz.eq.nunitz) then
         if(i.eq.1.or.i.eq.2.or.i.eq.3.or.i.eq.4) then
            call random_number(x)
            if(x*100.0d0.lt.defectrate) then
               Klcr(1,ir)=-Ku
               do ik=2,9
                  Klcr(ik,ir)=0.0d0
               end do
            end if
         end if

      end if
   end do
else
   do ir=1,nrtot
      i=irsub(ir)
      ix=irunit(ir,1)
      iy=irunit(ir,2)
      iz=irunit(ir,3)
      if(iz.eq.nunitz) then

         if(i.eq.1.or.i.eq.2.or.i.eq.3.or.i.eq.4) then
            Klcr(1,ir)=-Ku
            do ik=2,9
               Klcr(ik,ir)=0.0d0
            end do
         else
            call random_number(x)
            if(x*100.0d0.lt.defectrate-100.0d0) then
               Klcr(1,ir)=-Ku
               do ik=2,9
                  Klcr(ik,ir)=0.0d0
               end do
            end if
         end if

      end if
   end do
end if

open(93,file="R_defect_K1_in_gb_phase.txt")
icount=0
do ir=1,nrtot
   if(Klcr(1,ir).lt.0.0d0) then
      write(93,*) irsub(ir),irunit(ir,1),irunit(ir,2),irunit(ir,3),Klcr(1,ir)/kb
      icount=icount+1
   end if
end do

end subroutine defectR_randomz


subroutine defectR_uniformz(Klcr,Ku,defectrate,nr,nunit,icount,irunit,irsub,nunitx,nunity,nunitz,nrtot,ranfix)
integer(4) :: nunit,nr,irunit(nr*nunit,3),irsub(nr*nunit)
integer(4) :: ix,iy,iz,nrtot,ir0,c
integer,allocatable :: seed(:)
integer(4) :: i1,i2,i3,i4,icount,seedsize,ic,nunitx,nunity,nunitz
real(8) :: Klcr(9,nr*nunit),Ku,defectrate
real(8),parameter :: mev=1.60218d-21 ! (Merg/mev) 
real(8),parameter :: kb=8.6173303d-2*mev ! (Merg/K)
real(8) :: x

do ir=1,nrtot
   i=irsub(ir)
   ix=irunit(ir,1)
   iy=irunit(ir,2)
   iz=irunit(ir,3)
   if(iz.eq.nunitz) then
      Klcr(1,ir)=Ku*(1.0d0-defectrate/100.0d0)
!      Klcr(1,ir)=-Ku*defectrate/200.0d0
      do ik=2,9
         Klcr(ik,ir)=0.0d0
      end do
   end if
end do

open(93,file="R_defect_K1_in_gb_phase.txt")
icount=0
do ir=1,nrtot
   if(Klcr(1,ir).lt.0.0d0) then
      write(93,*) irsub(ir),irunit(ir,1),irunit(ir,2),irunit(ir,3),Klcr(1,ir)/kb
      icount=icount+1
   end if
end do

end subroutine defectR_uniformz




subroutine surface001(Klcr,Klcr001,Ku,defectrate,nr,nunit,icount,irunit,irsub,nunitx,nunity,nunitz,nrtot,ranfix,iran)
integer(4) :: nunit,nr,irunit(nr*nunit,3),irsub(nr*nunit)
integer(4) :: ix,iy,iz,nrtot,ir0,c,iran
integer,allocatable :: seed(:)
integer(4) :: i1,i2,i3,i4,icount,seedsize,ic,nunitx,nunity,nunitz
real(8) :: Klcr(9,nr*nunit),Klcr001(9,nr),Ku,defectrate
real(8),parameter :: mev=1.60218d-21 ! (Merg/mev) 
real(8),parameter :: kb=8.6173303d-2*mev ! (Merg/K)
real(8) :: x
character(10) :: ranfix

! random
call random_seed(size=seedsize)
allocate(seed(seedsize))
if(ranfix.eq."on") then
   call random_seed(put=(/iran/))
else
   call system_clock(count=c)
   call random_seed(put=(/c/))
end if

   do ir=1,nrtot
      i=irsub(ir)
      ix=irunit(ir,1)
      iy=irunit(ir,2)
      iz=irunit(ir,3)
      if(iz.eq.nunitz) then
         if(i.le.8) then
            do ik=1,9
               Klcr(ik,ir)=Klcr001(ik,i)
            end do
         end if
      end if
   end do

open(93,file="R_defect_K1_in_gb_phase.txt")
icount=0
do ir=1,nrtot
   if(Klcr(1,ir).lt.0.0d0) then
      write(93,*) irsub(ir),irunit(ir,1),irunit(ir,2),irunit(ir,3),Klcr(1,ir)/kb
      icount=icount+1
   end if
end do

end subroutine surface001


