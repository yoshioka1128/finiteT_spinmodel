program Intopencore
implicit none
integer, parameter :: iatommax=10,inoneqmax=10
real(8) :: a,b,c,d,w(195*4+1),r4f(781),x(0:781),clm(-6:6,-6:6)
real(8) :: dr(781+1),rmt,r0=0.00001d0,dx
real(8) :: Alm(49,iatommax),Almav(-6:6,0:6,2)
real(8) :: pi=dacos(-1.0d0),abc
integer(4) :: i,fi,j,k,l,m,lmax,iabc,inum,iatom(inoneqmax),inoneq,itot
character(len=32) atom,filevcoul
character(32) rmtch

! set alm table
clm=0.0d0
clm(2,0)=dsqrt(5.0d0/pi)/4.0d0

clm(2,2)=dsqrt(15.0d0/pi)/4.0d0
clm(-2,2)=dsqrt(15.0d0/pi)/4.0d0

clm(4,0)=dsqrt(1.0d0/pi)*3.0d0/16.0d0

clm(4,2)=dsqrt(5.0d0/pi)*3.0d0/8.0d0
clm(-4,2)=dsqrt(5.0d0/pi)*3.0d0/8.0d0

clm(4,3)=dsqrt(70.0d0/pi)*3.0d0/8.0d0
clm(-4,3)=dsqrt(70.0d0/pi)*3.0d0/8.0d0

clm(4,4)=dsqrt(35.0d0/pi)*3.0d0/16.0d0
clm(-4,4)=dsqrt(35.0d0/pi)*3.0d0/16.0d0

clm(6,0)=dsqrt(13.0d0/pi)/32.0d0

clm(6,2)=dsqrt(2730d0/pi)/64.0d0
clm(-6,2)=dsqrt(2730d0/pi)/64.0d0

clm(6,3)=dsqrt(2730d0/pi)/32.0d0
clm(-6,3)=dsqrt(2730d0/pi)/32.0d0

clm(6,4)=dsqrt(13.0d0/7.0d0/pi)*21.0d0/32.0d0
clm(-6,4)=dsqrt(13.0d0/7.0d0/pi)*21.0d0/32.0d0

clm(6,6)=dsqrt(26.0d0/231.0d0/pi)*231.0d0/64.0d0
clm(-6,6)=dsqrt(26.0d0/231.0d0/pi)*231.0d0/64.0d0

write(6,*) "for filename"
write(6,*) "rare earth ion ?"
read(5,*) atom
write(6,*) "Rmt?"
read(5,*) rmt
write(rmtch,'(f5.1)') rmt
write(6,*) "filename for vcoul"
read(5,*) filevcoul
write(6,*) "number of non-equivalent R ion"
read(5,*) inoneq

do i=1,inoneq
   write(6,*) "number of each R non-equivalent ion",i
   read(5,*) iatom(i)
end do

dx=(1.0d0/780.0d0)*log(rmt/r0)
! read file
open(10,file=''//trim(adjustl(filevcoul))//'.vcoul')
open(fi,file=''//trim(adjustl(atom))//''//trim(adjustl(rmtch))//'_R4f.dat')
! write file
open(11,file='Alm_for_#of'//trim(adjustl(atom))//'rmt'//trim(adjustl(rmtch))//'.txt')
open(12,file='Alm_for_lm_'//trim(adjustl(atom))//'rmt'//trim(adjustl(rmtch))//'.txt')
open(15,file='Alm_ave_'//trim(adjustl(atom))//'rmt'//trim(adjustl(rmtch))//'.txt')
open(16,file='cfp'//trim(adjustl(atom))//'opencore.txt')
!open(13,file='integrant.txt')
!open(14,file='Coulpot_RFeN66.txt')

x(0)=0.0d0
do i=1,781
   read(fi,*) abc,r4f(i)
 ! set x,dr
   x(i)=r0*dexp(dx*dble(i-1))
   dr(i)=r0*(dexp(dx*dble(i-1))-dexp(dx*dble(i-2)))
end do
dr(1)=r0
dr(782)=0.0d0

!do k=1,49
itot=0.0d0
do i=1,inoneq
   itot=itot+iatom(i)
end do

write(6,*) "total number of R ions",itot

do inum=1,itot
   read(10,*)
   read(10,*)
   read(10,*)
   read(10,*)
   read(10,*)

   write(11,*) 
   write(11,*) 'number of R',inum
   iabc=0
   do lmax=0,6
      l=lmax
      m=0
      iabc=iabc+1
      call integral(Alm,r4f,dr,x,pi,iabc,inum,l,m,clm)
      do m=1,lmax,1
         do k=1,-1,-2
            l=lmax*k
            iabc=iabc+1
            call integral(Alm,r4f,dr,x,pi,iabc,inum,l,m,clm)
         end do
      end do
   end do
   read(10,*)
   read(10,*)
   read(10,*)
end do

if(inoneq.eq.1) then
   iabc=0
   do lmax=0,6
      l=lmax
      m=0
      write(12,*) 'l,m=',l,m
      write(12,*) 'number of R,      Alm<r^l> (K)'
      iabc=iabc+1
      
      abc=0.0d0
      do inum=1,iatom(1)
         write(12,*) inum,Alm(iabc,inum)
         abc=abc+dabs(Alm(iabc,inum))
      end do
      abc=abc/dble(iatom(1))
      if(clm(l,m).ne.0.0d0) write(15,*) 'l,m=',l,m
      if(clm(l,m).ne.0.0d0) write(15,*) 'R',abc*Alm(iabc,1)/dabs(Alm(iabc,1))
      
      do m=1,lmax,1
         do k=1,-1,-2
            l=lmax*k
            write(12,*) 'l,m=',l,m
            write(12,*) 'number of R,      Alm<r^l> (K)'
            iabc=iabc+1
            abc=0.0d0
            do inum=1,iatom(1)
               write(12,*) inum,Alm(iabc,inum)
               abc=abc+dabs(Alm(iabc,inum))
            end do
            abc=abc/dble(iatom(1))
            if(clm(l,m).ne.0.0d0) write(15,*) 'l,m=',l,m
            if(clm(l,m).ne.0.0d0) write(15,*) 'R',abc*Alm(iabc,1)/dabs(Alm(iabc,1))
         end do
      end do
   end do

else if(inoneq.eq.2) then
   iabc=0
   do lmax=0,6
      l=lmax
      m=0
      write(12,*) 'l,m=',l,m
      write(12,*) 'number of R,      Alm<r^l> (K)'
      iabc=iabc+1
      
      abc=0.0d0
!      do inum=5,8
      do inum=itot/inoneq+1,itot
         write(12,*) inum,Alm(iabc,inum)
         abc=abc+dabs(Alm(iabc,inum))
      end do
      abc=abc/dble(iatom(2))
      if(clm(l,m).ne.0.0d0) write(15,*) 'l,m=',l,m
      if(clm(l,m).ne.0.0d0) write(15,*) '4f1',abc*Alm(iabc,itot)/dabs(Alm(iabc,itot))
      Almav(l,m,2)=abc*Alm(iabc,itot)/dabs(Alm(iabc,itot))
      if(clm(l,m).ne.0.0d0) write(15,*) '4f2',abc*Alm(iabc,itot/inoneq+2)/dabs(Alm(iabc,itot/inoneq+2))
      abc=0.0d0
!      do inum=1,4
      do inum=1,itot/inoneq
         write(12,*) inum,Alm(iabc,inum)
         abc=abc+dabs(Alm(iabc,inum))
      end do
      abc=abc/dble(iatom(1))
      if(clm(l,m).ne.0.0d0) write(15,*) '4g1',abc*Alm(iabc,itot/inoneq)/dabs(Alm(iabc,itot/inoneq))
      Almav(l,m,1)=abc*Alm(iabc,itot/inoneq)/dabs(Alm(iabc,itot/inoneq))
      if(clm(l,m).ne.0.0d0) write(15,*) '4g2',abc*Alm(iabc,2)/dabs(Alm(iabc,2))
      
      
      do m=1,lmax,1
         do k=1,-1,-2
            l=lmax*k
            write(12,*) 'l,m=',l,m
            write(12,*) 'number of R,      Alm<r^l> (K)'
            iabc=iabc+1
            abc=0.0d0
!            do inum=5,8
            do inum=itot/inoneq+1,itot
               write(12,*) inum,Alm(iabc,inum)
               abc=abc+dabs(Alm(iabc,inum))
            end do
            abc=abc/dble(iatom(2))
            if(clm(l,m).ne.0.0d0) write(15,*) 'l,m=',l,m
            if(clm(l,m).ne.0.0d0) write(15,*) '4f1',abc*Alm(iabc,itot)/dabs(Alm(iabc,itot))
            Almav(l,m,2)=abc*Alm(iabc,itot)/dabs(Alm(iabc,itot))
            if(clm(l,m).ne.0.0d0) write(15,*) '4f2',abc*Alm(iabc,itot/inoneq+2)/dabs(Alm(iabc,itot/inoneq+2))
            
            abc=0.0d0
!            do inum=1,4
            do inum=1,itot/inoneq
               write(12,*) inum,Alm(iabc,inum)
               abc=abc+dabs(Alm(iabc,inum))
            end do
            abc=abc/dble(iatom(1))
            if(clm(l,m).ne.0.0d0) write(15,*) '4g1',abc*Alm(iabc,itot/inoneq)/dabs(Alm(iabc,itot/inoneq))
            Almav(l,m,1)=abc*Alm(iabc,itot/inoneq)/dabs(Alm(iabc,itot/inoneq))
            if(clm(l,m).ne.0.0d0) write(15,*) '4g2',abc*Alm(iabc,2)/dabs(Alm(iabc,2))
            
         end do
      end do
   end do
end if

do i=1,2
   write(16,'(a5, i2)') atom,i
   m=0
   do l=2,6,2
   write(16,*) Almav(l,0,i),l,m
      do m=2,l,2
         if(mod(m/2,2).eq.0) then
            write(16,*) Almav(l,m,i),l,m
         else
            write(16,*) Almav(-l,m,i),l,-m
         end if
      end do
   end do
end do

stop
end program Intopencore


subroutine integral(Alm,r4f,dr,x,pi,iabc,inum,l,m,clm)
implicit none
real(8) :: Alm(49,8),r4f(781),w(195*4+1),dr(781+1),x(0:781),clm(-6:6,-6:6)
real(8) :: a,b,c,d,pi
integer(4) :: i,j,iabc,inum,l,m

      write(11,*) 'l,m=',l,m
!      write(13,*)
!      write(13,'(a,i1,a,i3,i3)') 'R',inum,',   l,m=',l,m
!      write(14,*)
!      write(14,'(a,i1,a,i3,i3)') 'R',inum,',   l,m=',l,m
      read(10,*)
      read(10,*)
      read(10,*)
      read(10,*)
      do i=1,195
         read(10,'(E22.12e2,E19.12e2,E19.12e2,E19.12e2)') a,b,c,d
         j=(i-1)*4+1
         w(j)=a
         w(j+1)=b
         w(j+2)=c
         w(j+3)=d
      end do
      read(10,'(E22.12e2)') w(195*4+1)
      
      a=0.0d0
      do i=1,781
         a=a+w(i)*r4f(i)*(dr(i)+dr(i+1))/2.0d0/x(i)**2
!         write(13,*) w(i)*r4f(i)/x(i)**2
!         write(14,*) x(i),w(i)/x(i)**2
      end do
    ! Ry=13.6057ev, ev=116049K
      Alm(iabc,inum)=a*13.6057d0*11604.9d0*clm(l,m)
      write(11,*) Alm(iabc,inum)
! ----------------------------------------

end subroutine integral
