subroutine Ki_analytical(sys,method,wwoK1TM,volume,lambda,Alm,rdch,K1TMT,HexT,MTMT,Tc,ratom,nr,dmult)
implicit none
integer(4),parameter :: nrmax=30,itempmax=1000,ianmax=2*10+1,mm=100
real(8),parameter :: pi=dacos(-1.0d0),Tmax=800.0d0,dT=2.0d0
real(8),parameter :: kb=1.3806488d-23 ! (J/K)
real(8),parameter :: muB=9.2740100783d-24 ! (Merg/kOe) or (kerg/Oe) or (J/T)
real(8),parameter :: mu0=4.0d0*pi*1.0d-7 ! (NA^-2) 
integer(4) :: nr,n4f,itemp,m,n,i,j,k,l,ii,jj,ll,ir,imax,ir0
real(8) :: abc,def,temp,S1,L1,J1,Hmax,MTM,Tc,Ms0,g1,dcoff,KuTM,volume,dmult,Dex(nrmax),Dso,xfp,Mtot
real(8) :: JM(mm,2),Stevens(6),xi(6),cclm(6,-6:6),lambda
real(8) :: extBJ(-1:8,nrmax),CC(7,-7:7,nrmax),disp(7,nrmax),rdcoff(0:6),Alm(6,-6:6,nrmax)
real(8) :: extBJ0(-1:8,nrmax),extLJ(-1:8,nrmax),MTMT(itempmax),HexT(nrmax,0:itempmax)
real(8) :: dKi1st(9),dKimix1(9),dKimix2(9),dKi2nd(9),dKidisp(9),K1TMT(itempmax)
real(8) :: dmag1st(2,nrmax),dmagmixex1(2,nrmax),dmagmixex2(2,nrmax)
real(8) :: dmagmixcf1(2,nrmax),dmagmixcf2(2,nrmax)
real(8) :: dmagcl1st(2,nrmax),dmagclmixcf1(2,nrmax),dmagclmixcf2(2,nrmax)
real(8) :: dmagclmixex1(2,nrmax),dmagclmixex2(2,nrmax),dclKi(9),dclKimix1(9),dclKimix2(9),eigen(ianmax,nrmax)
real(8) :: coffKi(9,6,-6:6),coffmix1(6,-1:8),coffmix2(6,-1:8),coffKi2nd(9,6,0:6),coffdisp(9,6),coffmag(6,-6:6)
real(8) :: Tmix1(6,nrmax),Tmix2(6,nrmax),T0mix1(6,nrmax),T0mix2(6,nrmax),Tclmix1(6,nrmax),Tclmix2(6,nrmax)
character(50) :: sys,ratom,method,wwoK1TM,rdch

dcoff=kb/volume*1.0d-6 ! convert from [K] to [MJ/m^3]
if(sys.eq."SmFe12_L3") then
   call rareearth_SmL3(ratom,n4f,0,imax,JM,S1,L1,J1)
else 
   call rareearth(ratom,n4f,0,imax,JM,S1,L1,J1)
end if
g1=1.0d0+(J1*(J1+1.0d0)+S1*(S1+1.0d0)-L1*(L1+1.0d0))/(2.0d0*J1*(J1+1.0d0))
Dso=lambda*(J1+1.0d0)

! matrix element of spherical tensor operators
if(sys.eq."SmFe12_L3") then
   call rdmatrix_SmL3(n4f,rdcoff,L1)
else
   call rdmatrix(n4f,rdcoff,L1) ! reduced matrix element
end if
call stfactor(Stevens,L1,J1,S1,rdcoff) ! Stevens factor
call xifactor(xi,L1,J1,S1,rdcoff) ! Stevens factor
call numfactor(cclm)
write(6,*) "Stevens factor and Xi factor"
do ll=2,6,2
   write(6,"(A,I5,2f25.15)") "l=",ll,Stevens(ll),xi(ll)
end do

! Ki, moment quantum
open(600,file=""//trim(sys)//"_"//trim(method)//"_Mtot_mix2_so_ex_rd"//trim(rdch)//".txt") 
open(160,file=""//trim(sys)//"_"//trim(method)//"_moment_mix2_so_ex_rd"//trim(rdch)//".txt") 
write(600,*) "# T [K],Mtot [K/atom]"
write(160,*) "# T [K],-2<S> [K/atom], -<L> [K/atom]"
! Ki, moment classical limit
open(55,file=""//trim(sys)//"_"//trim(method)//"_"//trim(wwoK1TM)//"_K1K2K3tot_1st_mix2_rd"//trim(rdch)//".txt")
write(54,*) "# T [K],Ms [T]"
write(55,*) "# T [K],K1 [MJ/m^3],K2 [MJ/m^3],K3 [MJ/m^3]"

do itemp=0,int(Tmax/dT)  ! loop for Temperature 
   temp=dble(itemp)*dT
   do ir=1,nr
      Dex(ir)=abs(-2.0d0*S1/(J1+1.0d0)*HexT(ir,itemp))
      call extendedBJ2(extBJ(:,ir),HexT(ir,itemp),temp,J1,JM,S1)
      call TVJl(Dex(ir),Dso,J1,L1,S1,Tmix1(:,ir),Tmix2(:,ir),HexT(ir,itemp),extBJ(:,ir))
      call extendedLJ2(extLJ(:,ir),HexT(ir,itemp),temp,J1,S1,Tc)
      call TVJl(Dex(ir),Dso,J1,L1,S1,Tclmix1(:,ir),Tclmix2(:,ir),HexT(ir,itemp),extLJ(:,ir))
   end do

   ! Fe sublattice
   MTM=MTMT(itemp) ! [muB/unitcell]
   if(wwoK1TM.eq."wK1Fe") then
      KuTM=K1TMT(itemp) ! [K/unitcell]
   else if(wwoK1TM.eq."woK1Fe") then
      KuTM=0.0d0
   else
      write(6,*) "read wwoK1TM error",wwoK1TM
      stop
   end if

   CC=0.0d0
   disp=0.0d0
   if(temp.le.Tc) then
      call interJ2nd(itemp,ianmax,nr,g1,CC,HexT,temp,eigen,J1,JM)
      call dispersion(itemp,ianmax,nr,g1,disp,HexT,temp,eigen,J1,JM)
   end if
   
   coffKi=0.0d0
   coffKi(1,2,0)=-3.0d0  ! contribute to K1
   coffKi(1,4,0)=-40.0d0
   coffKi(1,6,0)=-168.0d0
   coffKi(2,4,0)=35.0d0 ! contribute to K2
   coffKi(2,6,0)=378.0d0
   coffKi(3,6,0)=-231.0d0 ! contribute to K3
   coffKi(4,4,4)=1.0d0
   coffKi(4,6,4)=10.0d0 ! contribute to K2^1
   call tlm(coffmag,0.0d0,0.0d0)
   
   coffKi2nd=0.0d0
   coffKi2nd(1,2,1)=     -36.0d0*(dsqrt(5.0d0/(8.0d0*pi))/cclm(2,1))**2
   coffKi2nd(1,2,2)=      0.0d0
   coffKi2nd(2,2,1)=      36.0d0*(dsqrt(5.0d0/(8.0d0*pi))/cclm(2,1))**2
   coffKi2nd(2,2,2)=-9.0d0/4.0d0*(dsqrt(5.0d0/(8.0d0*pi))/cclm(2,2))**2
   if(temp.eq.0) then
      coffdisp=0.0d0
   else
      coffdisp(1,2)=-0.5d0/temp*(-12.0d0)
      coffdisp(2,2)=-0.5d0/temp*(9.0d0)
   end if
   
   dKi1st=0.0d0
   dKimix1=0.0d0
   dKimix2=0.0d0
   
   dclKi=0.0d0
   dclKimix1=0.0d0
   dclKimix2=0.0d0
   dKi2nd=0.0d0
   dKidisp=0.0d0
   
   dmagmixcf1=0.0d0
   dmagmixcf2=0.0d0
   dmagmixex1=0.0d0
   dmagmixex2=0.0d0
   do ir=1,nr ! inequivalent atom 
      abc=(2.0d0*(L1+1.0d0))/(3.0d0*(J1+1.0d0))
      ! moment w/o mixing 
      dmag1st(1,ir) =-2.0d0*S1/(J1+1.0d0)*extBJ(1,1)
      dmag1st(2,ir) =(L1+1.0d0)/(J1+1.0d0)*extBJ(1,1)
      ! moment w/ mixing 
      dmagmixex1(1,ir)=-2.0d0*abc*Tmix1(1,ir)
      dmagmixex2(1,ir)=-2.0d0*abc*Tmix2(1,ir)
      dmagmixex1(2,ir)=abc*Tmix1(1,ir)
      dmagmixex2(2,ir)=abc*Tmix2(1,ir)
      ! classical limit moment w/o mixing
      dmagcl1st(1,ir)   =-2.0d0*S1/(J1+1.0d0)*extLJ(1,1)
      dmagcl1st(2,ir)   =(L1+1.0d0)/(J1+1.0d0)*extLJ(1,1)
      ! classical limit moment w/ mixing
      dmagclmixex1(1,ir)=-2.0d0*abc*Tclmix1(1,ir)
      dmagclmixex2(1,ir)=-2.0d0*abc*Tclmix2(1,ir)
      dmagclmixex1(2,ir)=abc*Tclmix1(1,ir)
      dmagclmixex2(2,ir)=abc*Tclmix2(1,ir)
      ! 1st order perturbation
      do l=2,6,2 ! for Alm(l)
         do m=-l,l
            do j=1,4 ! K1, K2, K3, K2^1
               ! Ki w/o mixing 
               dKi1st(j)=dKi1st(j) +coffKi(j,l,m)*Stevens(l)*Alm(l,m,ir)*extBJ(l,ir) ! w/o mixing
               dclKi(j)=dclKi(j)   +coffKi(j,l,m)*Stevens(l)*Alm(l,m,ir)*extLJ(l,ir)
               ! Ki w/ J mixing from Hex*Hcf
               abc=xi(l)*dble(l*(l+1))/(dble(2*l+1))*Alm(l,m,ir)
               dKimix1(j)=dKimix1(j)     +abc*coffKi(j,l,m)*Tmix1(l,ir) ! Kuz'min (2002) and Magnani (2003) w/ mixing
               dKimix2(j)=dKimix2(j)     +abc*coffKi(j,l,m)*Tmix2(l,ir)
               dclKimix1(j)=dclKimix1(j) +abc*coffKi(j,l,m)*Tclmix1(l,ir)
               dclKimix2(j)=dclKimix2(j) +abc*coffKi(j,l,m)*Tclmix2(l,ir)
               ! magnetic moment w/ J mixing from Hex*Hcf
               def=2.0d0*S1/(J1+1.0d0)*(1.0d0/Dex(ir))*abc
               dmagmixcf1(1,ir)=dmagmixcf1(1,ir)       +def*coffmag(l,m)*Tmix1(l,ir)
               dmagmixcf2(1,ir)=dmagmixcf2(1,ir)       +def*coffmag(l,m)*Tmix2(l,ir)
               dmagmixcf1(2,ir)=dmagmixcf1(2,ir) -0.5d0*def*coffmag(l,m)*Tmix1(l,ir)
               dmagmixcf2(2,ir)=dmagmixcf2(2,ir) -0.5d0*def*coffmag(l,m)*Tmix2(l,ir)
            end do ! do j
         end do ! do m
      end do ! do l 
      do j=1,4
         do i=-2,2 ! second order perturbation
            if(i.eq.0) cycle
            dKi2nd(j)=dKi2nd(j)+coffKi2nd(j,2,abs(i))*CC(2,-i,ir)*(Stevens(2)*Alm(2,0,ir))**2
         end do
         dKidisp(j)=dKidisp(j)+coffdisp(j,2)*(disp(2,ir)-extBJ(2,ir)**2)*(Stevens(2)*Alm(2,0,ir))**2
      end do ! for j
   end do ! for ir

   abc=0.0d0
   do ir=1,nr
      do i=1,2
         abc=abc+(dmag1st(i,ir)+dmagmixex1(i,ir)+dmagmixex2(i,ir))
      end do
   end do
   Mtot=abc*dmult+MTM

   write(600,*) temp,mu0*Mtot*muB/volume


    ir=1
   ! total anisotropy constants
   write(55,"(4f25.15)") temp,&
       ((dKi1st(1)+dKimix1(1)+dKimix2(1))*dmult+KuTM)*dcoff,&
       (dKi1st(2)+dKimix1(2)+dKimix2(2))*dmult*dcoff,&
       (dKi1st(3)+dKimix1(3)+dKimix2(3))*dmult*dcoff ! K1,K2 [MJ/m^3]

   ! magnetic moment
   write(160,*) temp,& ! w/ mixng , no CF
        (dmag1st(1,ir)+dmagmixex1(1,ir)+dmagmixex2(1,ir)),&
        (dmag1st(2,ir)+dmagmixex1(2,ir)+dmagmixex2(2,ir))
   
   ! Nucleation field, FP field, MA fiel
   abc=(dKi1st(2)+dKimix1(2)+dKimix2(2))*dmult/((dKi1st(1)+dKimix1(1)+dKimix2(1))*dmult+KuTM)
   xfp=(-1.0d0+dsqrt(-3.0d0/abc-2.0d0))/3.0d0
end do ! itemp

end subroutine Ki_analytical
