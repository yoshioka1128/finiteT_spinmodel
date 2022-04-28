subroutine fitting_Kuzmin
implicit none
real(8) :: temp,Ms,s,Tc,MsT,MsT_Kuzmin
integer(4) :: i

! fitting @ 300K
open(44,file="PrFe105V15_MsT.txt")
open(45,file="PrFe105V15N_MsT.txt")
temp=300.0d0 ! for PrFe105V15
do i=0,1000,10
   temp=dble(i)
   Ms=131.4 ! 121.4 @T=300
   s=0.30d0
   Tc=625.0d0
   MsT=MsT_Kuzmin(temp,Ms,s,Tc)
   write(6,*) temp,MsT
   write(44,*) temp,MsT
end do

write(6,*) 
temp=300.0d0 ! for PrFe105V15N
do i=0,1000,10
   temp=dble(i)
   Ms=157.5 ! 142.8 @T=300K
   s=1.24d0
   Tc=820.0d0
   MsT=MsT_Kuzmin(temp,Ms,s,Tc)
   write(6,*) temp,MsT
   write(45,*) temp,MsT
end do
stop

end subroutine fitting_Kuzmin
