!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine MF_uniform(nfe,nr,nunit,temp,magFe,magR,JRFe,JFeR,JFeFe,JRR,Ku,KuFe)
implicit none
real(8),parameter :: ax=8.81d0,cx=12.21d0 ! (angstrom)
real(8),parameter :: mev=1.60218d-21 ! (Merg/mev) ! 1 (erg)=1.0d-7 (J)=6.24151d11 (ev)
real(8),parameter :: kb=8.6173303d-2*mev ! (Merg/K)
real(8),parameter :: joule=10.0d0 ! (Merg)
real(8),parameter :: muB=9.274d-24 ! (Merg/kOe) or (kerg/Oe)
real(8),parameter :: d=4.0d0*1.0d24/(ax**2*cx)*muB ! (kOe) or (kerg/Oe/cm^3)
real(8) :: temp,gf=2.0d0,gr=8.0d0/11.0d0,magR(nr*nunit),magFe(nfe*nunit),Hm
real(8) :: nff,nrr,nrf
real(8) :: Zff,Zrf,Zrr
real(8) :: alphaf,betar,Jr,Jf
real(8) :: JRFe,JFeR,JFeFe,JRR,Ku,KuFe
integer(4) :: nfe,nr,nunit

! Herbst et al.
if(temp.eq.200.0d0) then
   magFe=2.05d0*muB
   magR=2.61d0*muB
   KuFe=1.05d0*joule*ax**2*cx*1.0d-30/dble(nfe) !(Merg) Hirosawa
   Ku=5.30948717919229*1.0d6*(ax**2*cx*1.0d-30)*joule/8.0d0 ! (Merg) Hm=350 K
else if(temp.eq.300.0d0) then
   magFe=1.96d0*muB
   magR=1.96d0*muB
   KuFe=1.13d0*joule*ax**2*cx*1.0d-30/dble(nfe)
   Ku=5.14040457839072d0*1.0d6*(ax**2*cx*1.0d-30)*joule/8.0d0 ! (Merg)
else if(temp.eq.400.0d0) then
   magFe=1.77d0*muB
   magR=1.52d0*muB
   KuFe=0.97d0*joule*ax**2*cx*1.0d-30/dble(nfe)
   Ku=4.23364279755664*1.0d6*(ax**2*cx*1.0d-30)*joule/8.0d0 ! (Merg) Hm=350 K
else 
   stop
end if

! bond number
Zff=10.5d0
Zrf=16.0d0
Zrr=2.5d0
! molecular field coefficients
nff=5868d0
nrf=2198.0d0
nrr=335.0d0

JRFe=14.0d0*d*muB/Zrf*2.0d0*(gr/(gr-1.0d0))*nrf ! d*muB (kOe) 
JFeFe=7.0d0*d*muB/Zff*4.0d0*nff
JRR=d*muB/Zrr*(gr/(gr-1.0d0))**2*nrr



write(6,*) "JFeFe, JRFe, JRR (mev)"
write(6,*) JFeFe/mev,JRFe/mev,JRR/mev
write(6,*) "molecular field for muR (K)"
write(6,*) Hm,(2.0d0*((gr-1.0d0)/gr)**2*JRR*Zrr*magR(1)+0.5d0*(gr-1.0d0)/gr*JRFe*Zrf*magFe(1))/muB/kb
write(6,*) "molecular field for 2*(g-1)*J (K)"
write(6,*) Hm*gr/(gr-1.0d0)/2.0d0,(JRFe/gf*16.0d0*magFe(1)+2.0d0*JRR*(gr-1.0d0)/gr*2.5d0*magR(1))/2.0d0/muB/kb

JRFe=JRFe*(gr-1)/(gf*gr)/muB**2 ! (kOe/muB^2)
JFeR=JRFe
JFeFe=2.0d0*JFeFe/(gf**2)/muB**2
JRR=2.0d0*JRR*((gr-1.0d0)/gr)**2/muB**2

Hm=d*(2.0d0*nrr*magFe(1)+14.0d0*nrf*magR(1))/kb ! d (kOe), magFe (muB), d*magFe (Merg), d*magFe/mev/kb (K)
JRFe=-2.4d0*(gr-1)/(gf*gr)*mev/muB**2
JFeR=-2.4d0*(gr-1)/(gf*gr)*mev/muB**2
JFeFe=2.0d0*3.8d0/(gf**2)*mev/muB**2
JRR=2.0d0*0.2d0*((gr-1.0d0)/gr)**2*mev/muB**2


!write(6,*) JRFe,JFeFE,JRR
!write(6,*) magFe(1),magR(1)

! for Tc
Jr=9.0d0/2.0d0 ! if Nd
Jf=1.0d0
alphaf=14.0d0*(gf*Jf)**2*d*muB/kb*(Jf+1.0d0)/(3.0d0*Jf)
betar=2.0d0*(gr*Jr)**2*d*muB/kb*(Jf+1.0d0)/(3.0d0*Jf)
write(6,*) "Tc (K) R2Fe14B"
write(6,*) 0.5d0*(alphaf*nff+betar*nrr+dsqrt((alphaf*nff-betar*nrr)**2+4.0d0*alphaf*betar*nrf**2))
write(6,*) "Tc (K) Fe2"
write(6,*) alphaf*nff
write(6,*) "Tc (K) R2"
write(6,*) betar*nrr

end subroutine MF_uniform
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

