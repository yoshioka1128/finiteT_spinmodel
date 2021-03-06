subroutine extendedLJ2(extLJ,HexT,temp,J1,S1,Tc)
real(8) :: extLJ(-1:8),HexT,temp,J1,S1,Tc,dnu,abc,ri1,ri2,rk1,rk2,rip,rkp
integer(4) :: i

! extended Langevin function
extLJ=0.0d0
if(temp.eq.0.0d0) then
      extLJ(1)=J1
      extLJ(2)=(3.0d0*J1**2-J1**2)/2.0d0
      extLJ(3)=(5.0d0-3.0d0)*J1**3/2.0d0
      extLJ(4)=(35.0d0*J1**4-(30.0d0*J1**2)*J1**2+(3.0d0*J1**4))/8.0d0
      extLJ(5)=(63.0d0-70.0d0+15.0d0)*J1**5/8.0d0
      extLJ(6)=(231.0d0-315.0d0+105.0d0-5.0d0)*J1**6/16.0d0
      extLJ(7)=(429.0d0-693.0d0+315.0d0-35.0d0)*J1**7/16
      extLJ(8)=(6435.0d0-12012.0d0+6930.0d0-1260.0d0+35.0d0)*J1**8/128
else if(temp.lt.Tc) then
   do i=1,8
      dnu=dble(i)+0.5d0
         abc=dabs(-2.0d0*S1*J1/(J1+1.0d0)*HexT/temp) ! Toga et al. (2016) use m=g_J*J
!         write(6,*) "check1",i
         call bessik(abc,dnu,ri1,rk1,rip,rkp)
!         write(6,*) "check2",i
         call bessik(abc,0.5d0,ri2,rk2,rip,rkp)
         extLJ(i)=J1**i*ri1/ri2
   end do
else 
   extLJ=0.0d0
end if
extLJ(0)=1.0d0
extLJ(-1)=0.0d0

end subroutine extendedLJ2
