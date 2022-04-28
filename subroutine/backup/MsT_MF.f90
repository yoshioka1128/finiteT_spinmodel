!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine MsT_MF(nfe,nr,nunit,temp,magFe,magR,JRFe,JFeR,JFeFe,JRR,Hex,Hex0,magR0,magFe0,magRT,magFeT,HexT)
implicit none
real(8),parameter :: ax=8.81d0,cx=12.21d0 ! (angstrom)
real(8),parameter :: mev=1.60218d-21 ! (Merg/mev) ! 1 (erg)=1.0d-7 (J)=6.24151d11 (ev)
real(8),parameter :: kb=8.6173303d-2*mev ! (Merg/K)
real(8),parameter :: joule=10.0d0 ! (Merg/J)
real(8),parameter :: muB=9.274d-24 ! (Merg/kOe) or (kerg/Oe)
real(8),parameter :: d=4.0d0*1.0d24/(ax**2*cx)*muB ! from f.u. to density (kOe) or (kerg/Oe/cm^3)
real(8),parameter :: pi=dacos(-1.0d0)
real(8) :: temp,temp0,gf=2.0d0,gs=2.0d0,magR(nr*nunit),magFe(nfe*nunit),Hex
real(8) :: nff,nrr,nrf,gr
real(8) :: Zff,Zrf,Zrr,HexT(700)
real(8) :: alphaf,betar,Jr,Jf
real(8) :: JRFe,JFeR,JFeFe,JRR,Brillouin,Langevin,muFe,mu,Hex0
real(8) :: HRT,HFeT,xR,xFe,abc,def,ghi,magFe0,magR0,magFeT,magRT,Stot,Ltot,Jtot,magFe1,magR1,magFe2,magR2,x
integer(4) :: nfe,nr,nunit,i,ite,itemp,ife,ir
character(10) :: atom

atom="Nd"
! Molecular field approx.
! NdFeB angular momentum
if(atom.eq."Nd") then
   Stot=3.0d0/2.0d0
   Ltot=6.0d0
   Jtot=dabs(Ltot-Stot)
else if(atom.eq."Sm") then
   Stot=5.0d0/2.0d0
   Ltot=5.0d0
   Jtot=dabs(Ltot-Stot)
end if
gr=1.5d0+(Stot*(Stot+1.0d0)-Ltot*(Ltot+1.0d0))/(2.0d0*Jtot*(Jtot+1.0d0)) ! rande=8.0d0/11.0d0 for Nd
write(6,*) atom,gr,Jtot

! magnetization at T=0
!magtot0=35.8d0 ! for Nd
!magR0=9.0d0/2.0d0*rande ! |gr*J| @T=0
!magFe0=(magtot0-magR0*2.0d0)/14.0d0 ! |gs*S| @T=0

! NdFeB molecular field coefficients by fitting from magnetization
nff=5650d0
nrf=(-Hex0*2.0d0*(gr-1.0d0)/gr*kb/muB)/14.0d0/magFe0/d ! input Hex(T=0)=260 K, Hex(T=290)=249 K
nrr=0.0d0
!nrr=100.0d0

write(6,*) "nrf,nff,nrr"
write(6,*) nrf,nff,nrr
magR1=magR0
magFe1=magFe0

! temperature dependence
open(11,file="magnetization_Hex_Tdept.txt")
do itemp=1,700
   temp0=dble(itemp)
   magR1=magR0
   magFe1=magFe0

   do ite=1,20000
!      write(6,"(I5,3F15.5)") ite,2.0d0*magR1,14.0d0*magFe1,magFe1*14.0d0+magR1*2.0d0

      HRT=d*(2.0d0*nrr*magR1+14.0d0*nrf*magFe1)
      HFeT=d*(14.0d0*nff*magFe1+2.0d0*nrf*magR1)

!write(6,*) muB*magR0*HRT/(kb*temp0),muB*magFe0*HFeT/(kb*temp)
      xR=muB*magR0*HRT/(kb*temp0)
      xFe=muB*magFe0*HFeT/(kb*temp0)

      magR1=magR0*Brillouin(Jtot,xR)
!      magFe1=magFe0*Langevin(xFe) ! J=\infty
      magFe1=magFe0*Brillouin(1.0d0,xFe) ! J=1
      if(abs(magR1-magR2)/magR0.lt.1.0d-10.and.abs(magFe1-magFe2)/magFe0.lt.1.0d-10) exit
      magR2=magR1
      magFe2=magFe1
   end do
!   write(6,"(I5,3F15.5)") ite,2.0d0*magR1,14.0d0*magFe1,magFe1*14.0d0+magR1*2.0d0

   write(11,"(F10.3,4F25.17)") temp0,&
        2.0d0*magR1,&
        14.0d0*magFe1,&
        magFe1*14.0d0+magR1*2.0d0,&
        d*(14.0d0*nrf*magFe1)*muB/kb*gr/(gr-1.0d0)/2.0d0
   HexT(itemp)=Hex0*magFe1/magFe0
!   write(6,*) itemp,HexT(itemp)
 !  write(6,*) "T, 2*muR, 14*muFe, mutot (mu_B/f.u.)"
!   write(6,"(4F15.5)") temp0,2.0d0*magR1,14.0d0*magFe1,magFe1*14.0d0+magR1*2.0d0
end do

! for given temperature
magR1=magR0
magFe1=magFe0
do ite=1,20000
!      write(6,"(I5,3F15.5)") ite,2.0d0*magR1,14.0d0*magFe1,magFe1*14.0d0+magR1*2.0d0
   HRT=d*(2.0d0*nrr*magR1+14.0d0*nrf*magFe1)
   HFeT=d*(14.0d0*nff*magFe1+2.0d0*nrf*magR1)
!write(6,*) muB*magR0*HRT/(kb*temp),muB*magFe0*HFeT/(kb*temp)
   xR=muB*magR0*HRT/(kb*temp)
   xFe=muB*magFe0*HFeT/(kb*temp)

   magR1=magR0*Brillouin(Jtot,xR)
!      magFe1=magFe0*Langevin(xFe) ! J=\infty
   magFe1=magFe0*Brillouin(1.0d0,xFe) ! J=1
   if(abs(magR1-magR2)/magR0.lt.1.0d-10.and.abs(magFe1-magFe2)/magFe0.lt.1.0d-10) exit
   magR2=magR1
   magFe2=magFe1
end do
write(6,*) "T, muR (mu_B), muFe (mu_B), mutot (mu_B/f.u.)"
write(6,"(4F15.5)") temp,magR1,magFe1,magFe1*14.0d0+magR1*2.0d0
magR=magR1*muB ! gr*J*muB
magRT=magR1*muB
magFe=magFe1*muB ! gs*S*muB
magFeT=magFe1*muB


!open(19,file="Brillouin_Langevin.txt")
!do i=1,500
!   x=dble(i)/100.0d0
!   abc=Brillouin(1.0d0,x)
!   def=Langevin(x)
!   ghi=Brillouin(Jtot,x)
!   write(19,"(4F15.5)") x,abc,ghi,def
!end do

! bond number
Zff=10.5d0
Zrf=16.0d0 ! for R
Zrr=2.5d0

! Interaction
JRFe=14.0d0*d*muB/Zrf*gf*(gr/(gr-1.0d0))*nrf ! d*muB (kOe*muB=Merg) 
JFeFe=7.0d0*d*muB/Zff*gf**2*nff
JRR=d*muB/Zrr*(gr/(gr-1.0d0))**2*nrr

!Hex=-d*(2.0d0*nrr*magR1+14.0d0*nrf*magFe1)*muB/kb*gr/(gr-1.0d0)/2.0d0 ! d (kOe), magFe (muB), d*magFe (Merg), d*magFe/mev/kb (K)
Hex=-((JRFe/gs)*magFe1/gs*Zrf+2.0d0*(JRR/gs)*(gr-1.0d0)/gr*magR1*Zrr)/kb

JRFe=JRFe*(gr-1)/(gf*gr)/muB**2 ! (Merg/muB^2)
JFeR=JRFe
JFeFe=JFeFe/(gf**2)/muB**2 ! 2*JFeFe is pair interaction
JRR=JRR*((gr-1.0d0)/gr)**2/muB**2 ! 2*JRR is pair interaction

end subroutine MsT_MF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function Brillouin(J,x)
real(8) :: x,J,Brillouin

Brillouin=(2.0d0*J+1.0d0)/(2.0d0*J)/dtanh((2.0d0*J+1.0d0)*x/(2.0d0*J))-1.0d0/(2.0d0*J)/dtanh(x/(2.0d0*J))

end function Brillouin
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function Langevin(x)
real(8) :: x,J,Langevin

Langevin=1.0d0/dtanh(x)-1.0d0/x

end function Langevin
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
