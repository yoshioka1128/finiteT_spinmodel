!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine MF_uniform_excitedJ(nfe,nr,nfe0,nr0,nff,nrf,nrf0,nrr,J1,gr,nunit,temp,magFe,magR,JRFe,JFeR,JFeFe,JRR,&
     Hex,Hex0,magR0,magFe0,magRT,magFeT,HexT,dconv,Zff,Zrf,Zrr)
implicit none
real(8),parameter :: mev=1.60218d-21 ! (Merg/mev) ! 1 (erg)=1.0d-7 (J)=6.24151d11 (ev)
real(8),parameter :: kb=8.6173303d-2*mev ! (Merg/K)
real(8),parameter :: joule=10.0d0 ! (Merg/J)
real(8),parameter :: muB=9.274d-24 ! (Merg/kOe) or (kerg/Oe)
real(8),parameter :: pi=dacos(-1.0d0)
integer(4) :: nfe,nr,nfe0,nr0,nunit,i,ite,itemp,ife,ir
real(8) :: temp,temp0,gf=2.0d0,gs=2.0d0,magR(nr*nunit),magFe(nfe*nunit),Hex
real(8) :: nff,nrr,nrf,gr,dconv,nff0,nrr0,nrf0
real(8) :: Zff,Zrf,Zrr,HexT(700)
real(8) :: alphaf,betar,Jr,Jf,theta
real(8) :: JRFe,JFeR,JFeFe,JRR,Brillouin,Langevin,muFe,mu,Hex0
real(8) :: HRT,HFeT,xR,xFe,abc,def,ghi,magFe0,magR0,magFeT,magRT,J1,magFe1,magR1,magFe2,magR2,x,magRS1

write(6,*) "nrf,nff,nrr"
write(6,*) nrf,nff,nrr
nrf0=nrf
nff0=nff
nrr0=nrr
!magR1=magR0
magFe1=magFe0

open(129,file="SmFe12_opencore_RE1_J_theta-phi_nJex2.txt")
open(130,file="SmFe12_opencore_RE1_S_theta-phi_nJex2.txt")
read(129,*)
read(130,*)

! temperature dependence
open(11,file="magnetization_Hex_Tdept.txt")
do itemp=0,700,10
   read(129,*) abc,theta,def,magR1
   magR1=magR1*dcos(theta*pi/180.0d0)
   read(130,*) abc,theta,def,magRS1
   magRS1=magRS1*dcos(theta*pi/180.0d0)
   if(abc.ne.dble(itemp)) then
      write(6,*) "read error"
   end if
   temp0=dble(itemp)
!   magR1=magR0
!   magFe1=magFe0
!   write(6,*) abc,magR1,magRS1
   nrf=nrf0*magRS1/magR1
   magFe2=magFe0
   magFe1=magFe0
   do ite=1,20000
!      write(6,"(I5,3F15.5)") ite,2.0d0*magR1,14.0d0*magFe1,magFe1*14.0d0+magR1*2.0d0
!      HRT=dconv*(dble(nr0)*nrr*magR1+dble(nfe0)*nrf*magFe1)
      HFeT=dconv*(dble(nfe0)*nff*magFe1+dble(nr0)*nrf*magR1)
!write(6,*) muB*magR0*HRT/(kb*temp0),muB*magFe0*HFeT/(kb*temp)
!      xR=muB*magR0*HRT/(kb*temp0)
      xFe=muB*magFe0*HFeT/(kb*temp0)

!      magR1=magR0*Brillouin(J1,xR)
!      magFe1=magFe0*Langevin(xFe) ! J=\infty
      magFe1=magFe0*Brillouin(1.0d0,xFe) ! J=1
      if(abs(magFe1-magFe2)/magFe0.lt.1.0d-10) exit
!      if(abs(magR1-magR2)/magR0.lt.1.0d-10.and.abs(magFe1-magFe2)/magFe0.lt.1.0d-10) exit
!      magR2=magR1
      magFe2=magFe1
   end do

   write(11,"(F10.3,4F25.17)") temp0,dble(nr0)*magR1,dble(nfe0)*magFe1,&
        magFe1*dble(nfe0)+magR1*dble(nr0),&
        dconv*(dble(nfe0)*nrf*magFe1)*muB/kb*gr/(gr-1.0d0)/2.0d0
   write(6,"(F10.3,4F25.17)") temp0,dble(nr0)*magR1,dble(nfe0)*magFe1,&
        magFe1*dble(nfe0)+magR1*dble(nr0),&
        dconv*(dble(nfe0)*nrf*magFe1)*muB/kb*gr/(gr-1.0d0)/2.0d0
   HexT(itemp)=Hex0*magFe1/magFe0
end do
stop

! for given temperature
magR1=magR0
magFe1=magFe0
do ite=1,20000
   HRT=dconv*(dble(nr0)*nrr*magR1+dble(nfe0)*nrf*magFe1)
   HFeT=dconv*(dble(nfe0)*nff*magFe1+dble(nr0)*nrf*magR1)
   xR=muB*magR0*HRT/(kb*temp)
   xFe=muB*magFe0*HFeT/(kb*temp)
   magR1=magR0*Brillouin(J1,xR)
   magFe1=magFe0*Brillouin(1.0d0,xFe) ! J=1
   if(abs(magR1-magR2)/magR0.lt.1.0d-10.and.abs(magFe1-magFe2)/magFe0.lt.1.0d-10) exit
   magR2=magR1
   magFe2=magFe1
end do
write(6,*) "T, muR (mu_B), muFe (mu_B), mutot (mu_B/f.u.)"
write(6,"(4F15.5)") temp,magR1,magFe1,magFe1*dble(nfe0)+magR1*dble(nr0)
magR=magR1*muB ! gr*J*muB
magRT=magR1*muB
magFe=magFe1*muB ! gs*S*muB
magFeT=magFe1*muB

! Interaction
JRFe=dble(nfe0)*dconv*muB/Zrf*gf*(gr/(gr-1.0d0))*nrf ! for S*S
write(6,*) dble(nfe0),dconv*muB,Zrf*gf*(gr/(gr-1.0d0)),nrf,JRFe ! for S*S
!stop
JFeFe=dble(nfe0)/2.0d0*dconv*muB/Zff*gf**2*nff ! for S*S
JRR=0.0d0 ! for S*S

!Hex=-dconv*(2.0d0*nrr*magR1+14.0d0*nrf*magFe1)*muB/kb*gr/(gr-1.0d0)/2.0d0 ! d (kOe), magFe (muB), dconv*magFe (Merg), dconv*magFe/mev/kb (K)
Hex=-((JRFe/gs)*magFe1/gs*Zrf+2.0d0*(JRR/gs)*(gr-1.0d0)/gr*magR1*Zrr)/kb

end subroutine MF_uniform_excitedJ
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
