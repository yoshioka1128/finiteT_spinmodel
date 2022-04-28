subroutine Tmix_analytical(ir0,sys,method,lambda,HexT0,ratom,Tc)
implicit none
integer(4),parameter :: nrmax=30,itempmax=1000,ianmax=2*10+1
real(8),parameter :: pi=dacos(-1.0d0),Tmax=800.0d0,dT=2.0d0
integer(4) :: n4f,itemp,m,n,i,j,k,l,ii,jj,ll,imax,ir0
real(8) :: abc,temp,S1,L1,J1,Hmax,Tc,Dex0,Dso,JM(100,2),lambda
real(8) :: extBJ(-1:8),extBJ0(-1:8),extLJ(-1:8),HexT0(0:itempmax)
real(8) :: Tmix1(6),Tmix2(6),T0mix1(6),T0mix2(6),Tclmix1(6),Tclmix2(6)
character(50) :: sys,ratom,method

if(sys.eq."SmFe12_L3") then
   call rareearth_SmL3(ratom,n4f,0,imax,JM,S1,L1,J1)
else 
   call rareearth(ratom,n4f,0,imax,JM,S1,L1,J1)
end if

! Brillouin function
open(215,file=""//trim(sys)//"_"//trim(method)//"_TmixC_per_T0mix_xdept.txt")
open(220,file=""//trim(sys)//"_"//trim(method)//"_TmixB_per_T0mix_xdept.txt")
open(218,file=""//trim(sys)//"_"//trim(method)//"_TmixC_per_B0_xdept.txt")
open(222,file=""//trim(sys)//"_"//trim(method)//"_TmixB_per_B0_xdept.txt")
open(219,file=""//trim(sys)//"_"//trim(method)//"_TmixC_xdept.txt")
open(221,file=""//trim(sys)//"_"//trim(method)//"_TmixB_xdept.txt")
open(216,file=""//trim(sys)//"_"//trim(method)//"_B_per_B0_even.txt")
open(217,file=""//trim(sys)//"_"//trim(method)//"_B_per_B0_odd.txt")
open(315,file=""//trim(sys)//"_"//trim(method)//"_TmixB_per_T0mixB_xdept.txt")
open(316,file=""//trim(sys)//"_"//trim(method)//"_TmixdC_per_T0mixB_xdept.txt")
open(317,file=""//trim(sys)//"_"//trim(method)//"_TmixC_per_T0mixB_xdept.txt")
! saturated magnetization, molecular field, anisotropy constants for Fe sublattice

Dso=lambda*(J1+1.0d0)
do itemp=0,int(Tmax/dT)  ! loop for Temperature 
   temp=dble(itemp)*dT
   Dex0=abs(-2.0d0*S1/(J1+1.0d0)*HexT0(itemp))
   call extendedBJ2(extBJ,HexT0(itemp),temp,J1,JM,S1)
   call TVJl(Dex0,Dso,J1,L1,S1,Tmix1,Tmix2,HexT0(itemp),extBJ)
   call extendedLJ2(extLJ,HexT0(itemp),temp,J1,S1,Tc)
   call TVJl(Dex0,Dso,J1,L1,S1,Tclmix1,Tclmix2,HexT0(itemp),extLJ)
   if(temp.eq.0.0d0) then
      do l=1,6
         T0mix1(l)=Tmix1(l)
         T0mix2(l)=Tmix2(l)
      end do
      extBJ0=extBJ
   end if

   write(215,"(6f12.6)") temp,temp/(J1*Dex0),&
        (Tmix1(1)+Tmix2(1))/(T0mix1(1)+T0mix2(1)),&
        (Tmix1(2)+Tmix2(2))/(T0mix1(2)+T0mix2(2)),&
        (Tmix1(4)+Tmix2(4))/(T0mix1(4)+T0mix2(4)),&
        (Tmix1(6)+Tmix2(6))/(T0mix1(6)+T0mix2(6))
   write(220,"(6f12.6)") temp,temp/(J1*Dex0),&
        (Tmix1(1))/(T0mix1(1)),&
        (Tmix1(2))/(T0mix1(2)),&
        (Tmix1(4))/(T0mix1(4)),&
        (Tmix1(6))/(T0mix1(6))
   write(218,"(6f12.6)") temp,temp/(J1*Dex0),&
        (Tmix1(1)+Tmix2(1))/extBJ0(1),&
        (Tmix1(2)+Tmix2(2))/extBJ0(2),&
        (Tmix1(4)+Tmix2(4))/extBJ0(4),&
        (Tmix1(6)+Tmix2(6))/extBJ0(6)
   write(222,"(6f12.6)") temp,temp/(J1*Dex0),&
        (Tmix1(1))/extBJ0(1),&
        (Tmix1(2))/extBJ0(2),&
        (Tmix1(4))/extBJ0(4),&
        (Tmix1(6))/extBJ0(6)
   write(219,"(6f12.6)") temp,temp/(J1*Dex0),&
        (Tmix1(1)+Tmix2(1)),&
        (Tmix1(2)+Tmix2(2)),&
        (Tmix1(4)+Tmix2(4)),&
        (Tmix1(6)+Tmix2(6))

   write(221,"(6f12.6)") temp,temp/(J1*Dex0),&
        (Tmix1(1)),&
        (Tmix1(2)),&
        (Tmix1(4)),&
        (Tmix1(6))

   write(216,"(6f12.6)") temp,temp/(J1*Dex0),&
        extBJ(2)/extBJ0(2),extBJ(4)/extBJ0(4)!,extBJ(6,ir0)/extBJ0(6,ir0),extBJ(8,ir0)/extBJ0(8,ir0)
   write(217,"(6f12.6)") temp,temp/(J1*Dex0),&
        extBJ(1)/extBJ0(1),extBJ(3)/extBJ0(3),extBJ(5)/extBJ0(5)!,extBJ(7,ir0)/extBJ0(7,ir0)


   write(315,"(6f12.6)") temp,temp/(J1*Dex0),&
        Tmix1(1)/T0mix1(1),Tmix1(2)/T0mix1(2),&
        Tmix1(4)/T0mix1(4),Tmix1(6)/T0mix1(6)
   write(316,"(6f12.6)") temp,temp/(J1*Dex0),&
        Tmix2(1)/T0mix1(1),Tmix2(2)/T0mix1(2),&
        Tmix2(4)/T0mix1(4),Tmix2(6)/T0mix1(6)
   write(317,"(6f12.6)") temp,temp/(J1*Dex0),&
        (Tmix1(1)+Tmix2(1))/T0mix1(1),(Tmix1(2)+Tmix2(2))/T0mix1(2),&
        (Tmix1(4)+Tmix2(4))/T0mix1(4),(Tmix1(6)+Tmix2(6))/T0mix1(6)
end do ! do itemp

end subroutine Tmix_analytical
