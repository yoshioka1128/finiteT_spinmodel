subroutine anisotropy_tetra(temp,nr,freeki,dtheta0,dK1,dK2,dK3,Kloc)
implicit none
integer(4),parameter :: nrmax=30
integer(4) :: ir,nr,ip
real(8) :: freeki(0:4,0:4,nrmax),dtheta0,temp,def
real(8) :: K1sub(0:4,nrmax),K2sub(0:4,nrmax),K3sub(0:4,nrmax),Kloc(11,nrmax),dKi(9),dK1(2),dK2(3),dK3(4)

Kloc=0.0d0
do ir=1,nr
   do ip=0,3 ! ip=0: phi=0, ip=1: phi=pi/8, ip=2: phi=pi/4, ip=3: phi=pi/12
      K1sub(ip,ir)=freeki(1,ip,ir)/dtheta0**2
      K2sub(ip,ir)=K1sub(ip,ir)/3.0d0+(freeki(2,ip,ir)-4.0d0*freeki(1,ip,ir))/12.0d0/dtheta0**4
      K3sub(ip,ir)=((freeki(3,ip,ir)-6.0d0*freeki(2,ip,ir)+15.0d0*freeki(1,ip,ir))/8.0d0/dtheta0**6&
           -2.0d0*K1sub(ip,ir)+30.0d0*K2sub(ip,ir))/45.0d0
   end do
end do
do ir=1,nr 
   Kloc(1,ir)=K1sub(2,ir) ! K1
   Kloc(6,ir)=K1sub(0,ir)-Kloc(1,ir) ! K11
! general version
   Kloc(2,ir)=(K2sub(0,ir)+K2sub(2,ir))/2.0d0 ! Kloc(2,ir)
   Kloc(7,ir)=dsqrt(2.0d0)*(K2sub(1,ir)-Kloc(2,ir)) ! K22
   Kloc(3,ir)=Kloc(2,ir)-K2sub(2,ir) ! K222
   Kloc(8,ir)=(K3sub(0,ir)-4.0d0*K3sub(3,ir)/dsqrt(3.0d0)&
        +dsqrt(2.0d0)*K3sub(1,ir)-(2.0d0/dsqrt(3.0d0)-1.0d0)*K3sub(2,ir))&
        /(2.0d0+dsqrt(2.0d0)-2.0d0*dsqrt(3.0d0)) ! K3
   Kloc(9,ir)=2.0d0/dsqrt(3.0d0)*(K3sub(3,ir)+0.5d0*K3sub(2,ir)-1.5d0*Kloc(8,ir)) ! K33 or K'
   Kloc(10,ir)=Kloc(8,ir)-K3sub(2,ir) ! K333 or K''
   Kloc(11,ir)=dsqrt(2.0d0)*(Kloc(8,ir)+Kloc(9,ir)/dsqrt(2.0d0)-K3sub(1,ir)) ! K3333 or K'''
!   write(6,*) "K2",K2sub(0,ir),K2sub(1,ir),K2sub(2,ir)
end do
!stop

dK1=0.0d0
dK2=0.0d0
dK3=0.0d0
do ir=1,nr
   dK1(1)=dK1(1)+K1sub(1,ir) ! K1*sin^2theta: K1sub(ip=1,ir)=freeki(1,ip=1,ir)/dtheta0**2
   dK2(1)=dK2(1)+K2sub(1,ir) ! K2*sin^4theta: 
   dK3(1)=dK3(1)+K3sub(1,ir)
end do


return
end subroutine anisotropy_tetra
