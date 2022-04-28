program K1
implicit none
real(8) J,S,Lang,alpha,beta,gamma,aleng,bleng,cleng,volume,f(6),A20
real(8),parameter :: kb=1.3806488d-23 ! (J/K)
character(4) rare
write(6,*) "R for RFe12"
read(5,*) rare
if(rare.eq."Pr") then
!Pr
   J=4.0d0
   S=1.0d0
   alpha=-2.0d0**2*13.0d0/(3.0d0**2*5.0d0**2*11.0d0)
   beta=-2.0d0**2/(3.0d0**2*5.0d0*11.0d0**2)
   gamma=2.0d0**4*17.0d0/(3.0d0**4*5.0d0*7.0d0*11.0d0**2*13.0d0)
   aleng=8.556d0 ! angstrom
   bleng=aleng
   cleng=4.773d0 ! angstrom
 
! Nd 
!   J=9.0d0/2.0d0
!   S=3.0d0/2.0d0
!   alpha=-7.0d0/(3.0d0**2*11.0d0**2)
!   aleng=8.537d0 ! angstrom
!   bleng=8.618d0
!   cleng=4.880d0 ! angstrom
else if (rare.eq."Sm") then
!Sm 
   J=5.0d0/2.0d0
   S=5.0d0/2.0d0
   alpha=13.0d0/(3.0d0**2*5.0d0*7.0d0)
   beta=2.0d0*13.0d0/(3.0d0**2*5.0d0*7.0d0*11.0d0)
   gamma=0.0d0
   aleng=8.425d0 ! angstrom
   bleng=aleng
   cleng=4.775d0 ! angstrom
else
   write(6,*) "cannot supprot",rare
   stop
end if


   Lang=J+S

   f(2)=J*(J-0.5d0)
   volume=2.0d0*kb/(aleng*bleng*cleng*1.0d-30)*1.0d-6 ! temperature to J/m^3
   write(6,*) "A20"
   read(5,*) A20
   write(6,*) -3.0d0*f(2)*alpha*A20*volume
end program K1

