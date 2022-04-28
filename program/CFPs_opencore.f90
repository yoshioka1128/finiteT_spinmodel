program Intopencore
implicit none
real(8) :: a,b,c,d,w(195*4+1),r4f(781),x(0:781),clm(-6:6,-6:6)
real(8) :: dr(781+1),rmt,r0=0.00001d0,dx 
real(8) :: Almav(-6:6,0:6,2)
!real(8) :: Alm(iLMmax,iatommax)
real(8),allocatable :: Alm(:,:)
real(8) :: pi=dacos(-1.0d0),abc
integer(4) :: i,fi,j,k,l,m,lmax,iabc,inum,iatommax,iatnum,iLMmax
character(len=32) :: atom,filename
character(32) :: rmtch,ach,bch,cch,dch,val
character(len=255) :: dir,dirname

call getcwd(dirname)       ! get PATH of current diretory
i=SCAN(dirname,'/',.TRUE.) ! final position of '/'
j=LEN_TRIM(dirname)        ! number of character
dir=dirname(i+1:j)

! set alm table
clm=0.0d0

clm(1,0)= dsqrt(3.0d0/pi)/2.0d0
clm(1,1)= dsqrt(3.0d0/pi)/2.0d0
clm(-1,1)=dsqrt(3.0d0/pi)/2.0d0

clm(2,0)= dsqrt(5.0d0/pi)/4.0d0 ! from Ylm
clm(2,1)= dsqrt(15.0d0/pi)/2.0d0 ! from Ylm*sqrt(2)
clm(-2,1)=dsqrt(15.0d0/pi)/2.0d0 ! from Ylm*sqrt(2)
clm(2,2)= dsqrt(15.0d0/pi)/4.0d0 ! from Ylm*sqrt(2)
clm(-2,2)=dsqrt(15.0d0/pi)/4.0d0 ! from Ylm*sqrt(2)

clm(3,0)= dsqrt(7.0d0/pi)/4.0d0
clm(3,1)= dsqrt(42.0d0/pi)/8.0d0
clm(-3,1)=dsqrt(42.0d0/pi)/8.0d0
clm(3,2)= dsqrt(105.0d0/pi)/4.0d0
clm(-3,2)=dsqrt(105.0d0/pi)/4.0d0
clm(3,3)= dsqrt(70.0d0/pi)/8.0d0
clm(-3,3)=dsqrt(70.0d0/pi)/8.0d0

clm(4,0)= dsqrt(1.0d0/pi)*3.0d0/16.0d0
clm(4,1)= dsqrt(10.0d0/pi)*3.0d0/8.0d0
clm(-4,1)=dsqrt(10.0d0/pi)*3.0d0/8.0d0
clm(4,2)= dsqrt(5.0d0/pi)*3.0d0/8.0d0
clm(-4,2)=dsqrt(5.0d0/pi)*3.0d0/8.0d0
clm(4,3)= dsqrt(70.0d0/pi)*3.0d0/8.0d0
clm(-4,3)=dsqrt(70.0d0/pi)*3.0d0/8.0d0
clm(4,4)= dsqrt(35.0d0/pi)*3.0d0/16.0d0
clm(-4,4)=dsqrt(35.0d0/pi)*3.0d0/16.0d0

clm(5,0)= dsqrt(11.0d0/pi)/16.0d0
clm(5,1)= dsqrt(165.0d0/pi)/16.0d0
clm(-5,1)=dsqrt(165.0d0/pi)/16.0d0
clm(5,2)= dsqrt(1155.0d0/pi)/8.0d0
clm(-5,2)=dsqrt(1155.0d0/pi)/8.0d0
clm(5,3)= dsqrt(770.0d0/pi)/32.0d0
clm(-5,3)=dsqrt(770.0d0/pi)/32.0d0
clm(5,4)= dsqrt(385.0d0/pi)*3.0d0/16.0d0
clm(-5,4)=dsqrt(385.0d0/pi)*3.0d0/16.0d0
clm(5,5)= dsqrt(154.0d0/pi)*3.0d0/32.0d0
clm(-5,5)=dsqrt(154.0d0/pi)*3.0d0/32.0d0

clm(6,0)= dsqrt(13.0d0/pi)/32.0d0
clm(6,1)= dsqrt(273.0d0/pi)/16.0d0
clm(-6,1)=dsqrt(273.0d0/pi)/16.0d0
clm(6,2)= dsqrt(2730d0/pi)/64.0d0
clm(-6,2)=dsqrt(2730d0/pi)/64.0d0
clm(6,3)= dsqrt(2730d0/pi)/32.0d0
clm(-6,3)=dsqrt(2730d0/pi)/32.0d0
clm(6,4)= dsqrt(13.0d0/7.0d0/pi)*21.0d0/32.0d0
clm(-6,4)=dsqrt(13.0d0/7.0d0/pi)*21.0d0/32.0d0
clm(6,5)= dsqrt(2002.0d0/pi)*3.0d0/32.0d0
clm(-6,5)=dsqrt(2002.0d0/pi)*3.0d0/32.0d0
clm(6,6)= dsqrt(26.0d0/231.0d0/pi)*231.0d0/64.0d0
clm(-6,6)=dsqrt(26.0d0/231.0d0/pi)*231.0d0/64.0d0

write(6,*) "rare earth ion ?"
read(5,*) atom
write(6,*) "number of R atoms"
read(5,*) iatommax
write(6,*) "Rmt?"
read(5,*) rmt
write(rmtch,'(f5.1)') rmt
write(6,*) "valency? (di, tri, tetera)"
read(5,*) val

! write file
open(11,file=''//trim(adjustl(dir))//'_'//trim(adjustl(val))//'_Alm_rmt'//trim(adjustl(rmtch))//'_all.txt')
open(12,file=''//trim(adjustl(dir))//'_'//trim(adjustl(val))//'_Alm_rmt'//trim(adjustl(rmtch))//'_even.txt')

dx=(1.0d0/780.0d0)*log(rmt/r0)
! read file
open(10,file=''//trim(adjustl(dir))//'.vcoul')
if(val.eq."tri") then
   open(90,file='/home/yoshioka/research/RE_orbital/'&
        //trim(adjustl(atom))//''//trim(adjustl(rmtch))//'_R4f.dat')
else if(val.eq."di") then
   open(90,file='/home/yoshioka/research/RE_orbital/'&
        //trim(adjustl(atom))//''//trim(adjustl(rmtch))//'_R4f_divalent.dat')
else if(val.eq."tetra") then
   open(90,file='home/yoshioka/research/RE_orbital/'&
        //trim(adjustl(atom))//''//trim(adjustl(rmtch))//'_R4f_tetravalent.dat')
end if

! for R4f
x(0)=0.0d0
do i=1,781
   read(90,*) abc,r4f(i)
 ! set x,dr
   x(i)=r0*dexp(dx*dble(i-1))
   dr(i)=r0*(dexp(dx*dble(i-1))-dexp(dx*dble(i-2)))
end do
dr(1)=r0
dr(782)=0.0d0

do k=1,iatommax
   write(6,*)
   read(10,*)
   read(10,*)
   read(10,*)
   read(10,*) ach,bch,iatnum
   write(6,*) trim(ach)," ",trim(bch),iatnum
   write(11,*) trim(ach)," ",trim(bch),iatnum
   write(12,*) trim(ach)," ",trim(bch),iatnum
   read(10,*) ach,bch,cch,iLMmax
   write(6,*) trim(ach)," ",trim(bch)," ",trim(cch),iLMmax
   write(11,*) trim(ach)," ",trim(bch)," ",trim(cch),iLMmax
   write(12,*) trim(ach)," ",trim(bch)," ",trim(cch),iLMmax
   allocate(Alm(iLMmax,iatommax))
   write(11,*) "Alm<r^l> (K)","l   ","m   "
   write(12,*) "Alm<r^l> (K)","l   ","m   "
   write(6,*) "Alm<r^l> (K)","l   ","m   "
   do j=1,iLMmax
      read(10,*)
      read(10,*)
      read(10,*) ach,bch,cch,l,dch,m
      read(10,*)
      call integral(Alm,r4f,dr,x,pi,j,iatnum,l,m,clm,iLMmax,iatommax)
   end do
   deallocate(Alm)
   read(10,*)
   read(10,*)
   read(10,*)
end do

stop
end program Intopencore


subroutine integral(Alm,r4f,dr,x,pi,iabc,inum,l,m,clm,iLMmax,iatommax)
implicit none
integer(4) :: i,j,iabc,inum,l,m,iLMmax,iatommax,l2,m2
real(8) :: Alm(iLMmax,iatommax),r4f(781),w(195*4+1),dr(781+1),x(0:781),clm(-6:6,-6:6)
real(8) :: a,b,c,d,pi

!write(6,*) "l,m=",l,m
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
end do
! Ry=13.6057ev, ev=116049K
Alm(iabc,inum)=a*13.6057d0*11604.9d0*clm(l,m)

if(l.le.0) then
   m2=-m
   l2=-l
else
   m2=m
   l2=l
end if

write(11,*) Alm(iabc,inum),l2,m2
if(mod(l,2).eq.0) then
   write(12,*) Alm(iabc,inum),l2,m2
   write(6,*) Alm(iabc,inum),l2,m2
end if

end subroutine integral
