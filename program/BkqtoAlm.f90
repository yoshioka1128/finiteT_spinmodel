program BkqtoAlm
! from Bkq (cm^-1) to Alm<r^l> (K)
integer(4) :: l,m,iunit
real(8) :: a,b,clm(6,-6:6),fct,pi,trK,alm(6,-6:6),Blm(6,-6:6)
real(8) :: jmul,kb=1.3806488d-23,f2,f4,f6,K19(9)
real(8) :: alpha,beta,gamma,volume,aleng,bleng,cleng,angle,dsign
character(20) :: rare

open(10,file="case_Almrl.txt")

pi=dacos(-1.0d0)
trK=1.43879d0 ! from cm-1 to K

write(6,*) "rare earth"
read(5,*) rare
write(6,*) rare

! Stevens factor
if (rare.eq."Nd") then
   jmul=9.0d0/2.0d0
   alpha=-7.0d0/(3.0d0**2*11.0**2)
   beta=-2.0d0**3*17.0d0/(3.0d0**3*11.0d0**3*13.0d0)
   gamma=-5.0d0*17.0d0*19.0d0/(3.0d0**3*7.0d0*11.0d0**3*13.0d0**2)
   aleng=8.81d0
   bleng=aleng
   cleng=12.21d0
   angle=pi/2.0d0
   iunit=4
else if (rare.eq."Dy") then
   jmul=15.0d0/2.0d0
   alpha=-2.0d0/(3.0d0**2*5.0d0*7.0d0)
   beta=-2.0d0**3/(3.0d0**3*5.0d0*7.0d0*11.0d0*13.0d0)
   gamma=2.0d0**2/(3.0d0**3*7.0d0*11.0d0**2*13.0**2)
   aleng=8.76d0
   bleng=aleng
   cleng=11.99d0
   angle=pi/2.0d0
   iunit=4
else if (rare.eq."Ho") then
   jmul=8.0d0
   alpha=-1.0d0/(2.0d0*3.0d0**2*5.0d0**2)
   beta=-1.0d0/(2.0d0*3.0d0*5.0d0*7.0d0*11.0d0*13.0d0)
   gamma=-5.0d0/(3.0d0**3*7.0d0*11.0d0**2*13.0d0**2)
   aleng=8.75d0
   bleng=aleng
   cleng=11.99d0
   angle=pi/2.0d0
   iunit=4
else if (rare.eq."Pr") then
   jmul=4.0d0
   alpha=-2.0d0*13.0d0/(3.0d0**2*5.0d0**2*11.0d0)
   beta=-2.0d0**2/(3.0d0**2*5.0d0*11.0d0**2)
   gamma=2.0d0**4*17.0d0/(3.0d0**4*5.0d0*7.0d0*11.0d0**2*13.0d0)
   aleng=8.81d0
   bleng=aleng
   cleng=12.27d0
   angle=pi/2.0d0
   iunit=4
else if (rare.eq."Tb") then
   jmul=6.0d0
   alpha=-1.0/(3.0d0**2*11.0d0)
   beta=2.0d0/(3.0d0**3*5.0d0*11.0d0**2)
   gamma=-1.0d0/(3.0d0**4*7.0d0*11.0d0**2*13.0d0)
   aleng=8.77d0
   bleng=aleng
   cleng=12.05d0
   angle=pi/2.0d0
   iunit=4
else if (rare.eq."Sm") then
   jmul=5.0d0/2.0d0
   alpha=13.0d0/(3.0d0**2*5.0d0*7.0d0)
   beta=2.0d0*13.0d0/(3.0d0**3*5.0d0*7.0d0*11.0d0)
   gamma=0.0d0
   aleng=8.72d0
   bleng=aleng
   cleng=12.63d0
   angle=2.0d0*pi/3.0d0
   iunit=6
end if
   write(6,*) "J=",jmul

f2=jmul*(jmul-0.5d0)
f4=jmul*(jmul-0.5d0)*(jmul-1.0d0)*(jmul-1.5d0)
f6=jmul*(jmul-0.5d0)*(jmul-1.0d0)*(jmul-1.5d0)*(jmul-2.0d0)*(jmul-2.5d0)

! per volume (MJ/m^3)
volume=kb*dble(iunit)/(aleng*bleng*cleng*dsin(angle)*10.0d-30)*1.0d-6

! set alm table
! set alm table
clm=0.0d0
clm(2,0)=dsqrt(5.0d0/pi)/4.0d0

clm(2,2)=dsqrt(15.0d0/pi)/4.0d0
clm(2,-2)=dsqrt(15.0d0/pi)/4.0d0

clm(4,0)=dsqrt(1.0d0/pi)*3.0d0/16.0d0

clm(4,2)=dsqrt(5.0d0/pi)*3.0d0/8.0d0
clm(4,-2)=dsqrt(5.0d0/pi)*3.0d0/8.0d0

clm(4,3)=dsqrt(70.0d0/pi)*3.0d0/8.0d0
clm(4,-3)=dsqrt(70.0d0/pi)*3.0d0/8.0d0

clm(4,4)=dsqrt(35.0d0/pi)*3.0d0/16.0d0
clm(4,-4)=dsqrt(35.0d0/pi)*3.0d0/16.0d0

clm(6,0)=dsqrt(13.0d0/pi)/32.0d0

clm(6,2)=dsqrt(2730d0/pi)/64.0d0
clm(6,-2)=dsqrt(2730d0/pi)/64.0d0

clm(6,3)=dsqrt(2730d0/pi)/32.0d0
clm(6,-3)=dsqrt(2730d0/pi)/32.0d0

clm(6,4)=dsqrt(13.0d0/7.0d0/pi)*21.0d0/32.0d0
clm(6,-4)=dsqrt(13.0d0/7.0d0/pi)*21.0d0/32.0d0

clm(6,6)=dsqrt(26.0d0/231.0d0/pi)*231.0d0/64.0d0
clm(6,-6)=dsqrt(26.0d0/231.0d0/pi)*231.0d0/64.0d0

do 
   read(8,*,end=2) l,m,a,b
   fct=dsqrt(4.0d0*pi/(2.0d0*dble(l)+1.0d0))*trK

   if(mod(l,2).eq.0) then  
      if(m.eq.0) then
         alm(l,m)=clm(l,m)*fct*a
! old
!      else if(m.gt.0) then
!         alm(l,m)=clm(l,m)*fct*a*dsqrt(2.0d0)
!      else
!         alm(l,m)=clm(l,m)*fct*b*dsqrt(2.0d0)
!      end if
! end old
      else if(m.lt.0) then
         alm(l,-m)=clm(l,m)*fct*a*dsqrt(2.0d0)
         alm(l,m)=clm(l,m)*fct*b*dsqrt(2.0d0)
      end if

      write(6,*) l,m,alm(l,m)
   end if
end do

2  do l=2,6,2
      write(66,*) alm(l,0),l,0
      do m=1,l
         write(66,*) alm(l,-m),l,-m
         write(66,*) alm(l,m),l,m
      end do
   end do

   Blm(2,0)=alm(2,0)*alpha
   Blm(2,-2)=alm(2,-2)*alpha
   Blm(4,0)=alm(4,0)*beta
   Blm(4,-2)=alm(4,-2)*beta
   Blm(4,4)=alm(4,4)*beta
   Blm(6,0)=alm(6,0)*gamma
   Blm(6,-2)=alm(6,-2)*gamma
   Blm(6,4)=alm(6,4)*gamma
   Blm(6,-6)=alm(6,-6)*gamma

   dsign=1.0d0
   K19(1)=-3.0d0*f2*Blm(2,0)-40.0d0*f4*Blm(4,0)-168.0d0*f6*Blm(6,0)
   K19(2)=35.0d0*f4*Blm(4,0)+378.0d0*f6*Blm(6,0)
   K19(3)=f4*Blm(4,4)+10.0d0*f6*Blm(6,4)
   K19(4)=-231.0d0*f6*Blm(6,0)
   K19(5)=-11.0d0*f6*Blm(6,4)
   K19(6)=dsign*(f2*Blm(2,-2)+6.0d0*f4*Blm(4,-2)+15.0d0*f6*Blm(6,-2))
   K19(7)=dsign*(-7.0d0*f4*Blm(4,-2)-48.0d0*f6*Blm(6,-2))
   K19(8)=dsign*33.0d0*f6*Blm(6,-2)
   K19(9)=dsign*f6*Blm(6,-6)
   do i=1,9
      write(67,*) K19(i)*volume,"(MJ/m^3)  ","  K",i
   end do

   write(6,*) "K1=",-3.0d0*f2*alpha*alm(2,0)*volume,"(MJ/m^3)"
   write(6,*) "K1=",(-3.0d0*f2*alpha*alm(2,0)-40.0d0*f4*beta*alm(4,0))*volume,"(MJ/m^3)"
   write(6,*) "K1=",(-3.0d0*f2*alpha*alm(2,0)-40.0d0*f4*beta*alm(4,0)-168.0d0*f6*gamma*alm(6,0))*volume,"(MJ/m^3)"
   write(66,*) "K1=",-3.0d0*f2*alpha*alm(2,0)*volume,"(MJ/m^3)"
   write(66,*) "K1=",(-3.0d0*f2*alpha*alm(2,0)-40.0d0*f4*beta*alm(4,0))*volume,"(MJ/m^3)"
   write(66,*) "K1=",(-3.0d0*f2*alpha*alm(2,0)-40.0d0*f4*beta*alm(4,0)-168.0d0*f6*gamma*alm(6,0))*volume,"(MJ/m^3)"
stop

end program BkqtoAlm
