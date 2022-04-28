subroutine energy_analytical(temp,sys,method,volume,lambda,Alm,KuTM,Hex,ratom,nr,dmult,rdch)
implicit none
integer(4),parameter :: nrmax=30,nJmax=21,mm=100
real(8),parameter :: pi=dacos(-1.0d0)
real(8),parameter :: kb=1.3806488d-23 ! (J/K)
real(8),parameter :: dtheta=0.001d0,dphi0=0.015d0
integer(4) :: nr,n4f,itemp,m,i,j,k,l,imax,ii,it,ip,ll,ir,itmax,ipmax
real(8) :: temp,S1,L1,J1,theta,phi,dcoff,Dso,Dex,abc,KuTM,g1,volume,dmult,dphi,lambda
real(8) :: JM(mm,2),Stevens(6),extBJ(-1:8,nrmax),rdcoff(0:6),Alm(6,-6:6,nrmax)
real(8) :: dmagRT0(3),Tmix1(6,nrmax),Tmix2(6,nrmax),free0(3,nrmax),freeR(3,nrmax),gibbs(3),Hex(nrmax)
real(8) :: engmixcf1(nrmax),engmixcf2(nrmax),engcf(nrmax),coffmag(6,-6:6),xi(6),sc(3)
real(8) :: englvcf(nJmax),eigen0(nJmax),CJst(nJmax,-1:8,-8:8),coffmix1(6,-1:8),coffmix2(6,-1:8)
real(8) :: englvmixex1(nJmax),englvmixex2(nJmax),englvmixcf1(nJmax),englvmixcf2(nJmax)
character(50) :: sys,ratom,method,Tch,rdch

dcoff=kb/volume*1.0d-6 ! convert from [K] to [MJ/m^3]
write(Tch,'(f5.1)') temp
Tch=adjustl(Tch)
if(sys.eq."SmFe12_L3") then
   call rareearth_SmL3(ratom,n4f,0,imax,JM,S1,L1,J1)
else 
   call rareearth(ratom,n4f,0,imax,JM,S1,L1,J1)
end if
g1=1.0d0+(J1*(J1+1.0d0)+S1*(S1+1.0d0)-L1*(L1+1.0d0))/(2.0d0*J1*(J1+1.0d0))

! matrix element of spherical tensor operators
if(sys.eq."SmFe12_L3") then
   call rdmatrix_SmL3(n4f,rdcoff,L1)
else
   call rdmatrix(n4f,rdcoff,L1) ! reduced matrix element
end if
call stfactor(Stevens,L1,J1,S1,rdcoff) ! Stevens factor
call xifactor(xi,L1,J1,S1,rdcoff) ! Stevens factor

! Spherical Tensor operator @ T=0
Dex=abs(-2.0d0*S1/(J1+1.0d0)*Hex(1))
Dso=lambda*(J1+1.0d0)
do ir=1,nr
   call extendedBJ2(extBJ(:,ir),Hex(ir),temp,J1,JM,S1)
   call TVJl(Dex,Dso,J1,L1,S1,Tmix1(:,ir),Tmix2(:,ir),Hex(ir),extBJ(:,ir))
end do

abc=(2.0d0*(L1+1.0d0))/(3.0d0*(J1+1.0d0))
dmagRT0=0.0d0
do ir=1,nr
   dmagRT0(1)=dmagRT0(1)+g1*extBJ(1,ir)
   dmagRT0(2)=dmagRT0(2)+g1*extBJ(1,ir)-abc*Tmix1(1,ir)
   dmagRT0(3)=dmagRT0(3)+g1*extBJ(1,ir)-abc*Tmix1(1,ir)-abc*Tmix2(1,ir)
end do

open(799,file=""//trim(sys)//"_"//trim(method)//"_T"//trim(Tch)//"_eng_an_phi_theta90_rd"//trim(rdch)//".txt")
itmax=int(0.5d0*pi/dtheta)
do it=0,itmax,itmax
   theta=0.5d0*pi*dble(it)/dble(itmax)
   if(it.ne.0) then
      dphi=dphi0/dsin(theta)
      ipmax=int(0.5d0*pi/dphi)+1
   else
      ipmax=1
   end if
   do ip=0,ipmax*2
      phi=0.5d0*pi*dble(ip)/dble(ipmax)
      call tlm(coffmag,theta,phi)
      engcf=0.0d0
      engmixcf1=0.0d0
      engmixcf2=0.0d0
      do ir=1,nr
         do l=2,6,2
            do m=-l,l
               engcf(ir)=engcf(ir)+coffmag(l,m)*Alm(l,m,ir)*Stevens(l)*extBJ(l,ir)
               abc=xi(l)*dble(l*(l+1))/(dble(2*l+1))*coffmag(l,m)*Alm(l,m,ir)
               engmixcf1(ir)=engmixcf1(ir)+abc*Tmix1(l,ir)
               engmixcf2(ir)=engmixcf2(ir)+abc*Tmix2(l,ir)
            end do
         end do
      end do
      freeR=0.0d0
      do ir=1,nr
         freeR(1,ir)=freeR(1,ir)+engcf(ir)
         freeR(2,ir)=freeR(2,ir)+engcf(ir)+engmixcf1(ir)
         freeR(3,ir)=freeR(3,ir)+engcf(ir)+engmixcf1(ir)+engmixcf2(ir)
      end do
      if(it.eq.0) free0=freeR
      gibbs=0.0d0
      do i=1,3
         abc=0.0d0
         do ir=1,nr
            freeR(i,ir)=freeR(i,ir)-free0(i,ir)
            abc=abc+freeR(i,ir)
         end do
         gibbs(i)=gibbs(i)+abc*dmult+KuTM*dsin(theta)**2
      end do
      if(it.eq.itmax) write(799,"(4f25.15)") phi,freeR(3,1),freeR(3,2),gibbs(3)*dcoff
   end do  ! ip theta
end do   ! it theta





if(sys.eq."SmFe12_L3") then
   call rdmatrix_SmL3(n4f,rdcoff,L1)
else
   call rdmatrix(n4f,rdcoff,L1) ! reduced matrix element
end if
call stfactor(Stevens,L1,J1,S1,rdcoff) ! Stevens factor
call xifactor(xi,L1,J1,S1,rdcoff) ! Stevens factor

! Spherical Tensor operator @ T=0
Dex=abs(-2.0d0*S1/(J1+1.0d0)*Hex(1))
Dso=lambda*(J1+1.0d0)
do ir=1,nr
   call extendedBJ2(extBJ(:,ir),Hex(ir),temp,J1,JM,S1)
   call TVJl(Dex,Dso,J1,L1,S1,Tmix1(:,ir),Tmix2(:,ir),Hex(ir),extBJ(:,ir))
end do

it=0
theta=0.0d0
phi=0.0d0
abc=(2.0d0*(L1+1.0d0))/(3.0d0*(J1+1.0d0))
dmagRT0=0.0d0
do ir=1,nr
   dmagRT0(1)=dmagRT0(1)+g1*extBJ(1,ir)
   dmagRT0(2)=dmagRT0(2)+g1*extBJ(1,ir)-abc*Tmix1(1,ir)
   dmagRT0(3)=dmagRT0(3)+g1*extBJ(1,ir)-abc*Tmix1(1,ir)-abc*Tmix2(1,ir)
end do

itmax=int(0.5d0*pi/dtheta)
open(800,file=""//trim(sys)//"_"//trim(method)//"_T"//trim(Tch)//"_eng_an_theta_phi22.5_rd"//trim(rdch)//".txt")
do it=0,itmax
   theta=0.5d0*pi*dble(it)/dble(itmax)
   call tlm(coffmag,theta,pi/8.0d0)   
   engcf=0.0d0
   engmixcf1=0.0d0
   engmixcf2=0.0d0
   do ir=1,nr
      do l=2,6,2
         do m=-l,l
            engcf(ir)=engcf(ir)+coffmag(l,m)*Alm(l,m,ir)*Stevens(l)*extBJ(l,ir)
            abc=xi(l)*dble(l*(l+1))/(dble(2*l+1))*coffmag(l,m)*Alm(l,m,ir)
            engmixcf1(ir)=engmixcf1(ir)+abc*Tmix1(l,ir)
            engmixcf2(ir)=engmixcf2(ir)+abc*Tmix2(l,ir)
         end do
      end do
   end do
   freeR=0.0d0
   do ir=1,nr
      freeR(1,ir)=freeR(1,ir)+engcf(ir)
      freeR(2,ir)=freeR(2,ir)+engcf(ir)+engmixcf1(ir)
      freeR(3,ir)=freeR(3,ir)+engcf(ir)+engmixcf1(ir)+engmixcf2(ir)
   end do
   if(it.eq.0) free0=freeR
   gibbs=0.0d0
   do i=1,3
      abc=0.0d0
      do ir=1,nr
         freeR(i,ir)=freeR(i,ir)-free0(i,ir)
         abc=abc+freeR(i,ir)
      end do
      gibbs(i)=gibbs(i)+abc*dmult+KuTM*dsin(theta)**2
   end do
   write(800,"(4f25.15)") theta,freeR(3,1),freeR(3,2),gibbs(3)*dcoff
end do   ! it theta

open(801,file=""//trim(sys)//"_"//trim(method)//"_T"//trim(Tch)//"_eng_an_theta_phi0_rd"//trim(rdch)//".txt")
do it=0,itmax
   theta=0.5d0*pi*dble(it)/dble(itmax)
   call tlm(coffmag,theta,0.0d0)   
   engcf=0.0d0
   engmixcf1=0.0d0
   engmixcf2=0.0d0
   do ir=1,nr
      do l=2,6,2
         do m=-l,l
            engcf(ir)=engcf(ir)+coffmag(l,m)*Alm(l,m,ir)*Stevens(l)*extBJ(l,ir)
            abc=xi(l)*dble(l*(l+1))/(dble(2*l+1))*coffmag(l,m)*Alm(l,m,ir)
            engmixcf1(ir)=engmixcf1(ir)+abc*Tmix1(l,ir)
            engmixcf2(ir)=engmixcf2(ir)+abc*Tmix2(l,ir)
         end do
      end do
   end do
   freeR=0.0d0
   do ir=1,nr
      freeR(1,ir)=freeR(1,ir)+engcf(ir)
      freeR(2,ir)=freeR(2,ir)+engcf(ir)+engmixcf1(ir)
      freeR(3,ir)=freeR(3,ir)+engcf(ir)+engmixcf1(ir)+engmixcf2(ir)
   end do
   if(it.eq.0) free0=freeR
   gibbs=0.0d0
   do i=1,3
      abc=0.0d0
      do ir=1,nr
         freeR(i,ir)=freeR(i,ir)-free0(i,ir)
         abc=abc+freeR(i,ir)
      end do
      gibbs(i)=gibbs(i)+abc*dmult+KuTM*dsin(theta)**2
   end do
   write(801,"(4f25.15)") theta,freeR(3,1),freeR(3,2),gibbs(3)*dcoff
end do   ! it theta




open(281,file=""//trim(sys)//"_"//trim(method)//"_englevel_T"//trim(Tch)//"_lwstJso_rd"//trim(rdch)//".txt")
open(275,file=""//trim(sys)//"_"//trim(method)//"_englevel_T"//trim(Tch)//"_lwstJsoex_rd"//trim(rdch)//".txt")
open(280,file=""//trim(sys)//"_"//trim(method)//"_englevel_T"//trim(Tch)//"_lwstJsoexcf_rd"//trim(rdch)//".txt")
open(276,file=""//trim(sys)//"_"//trim(method)//"_englevel_T"//trim(Tch)//"_mixex1_rd"//trim(rdch)//".txt")
open(277,file=""//trim(sys)//"_"//trim(method)//"_englevel_T"//trim(Tch)//"_mixex2_rd"//trim(rdch)//".txt")
open(278,file=""//trim(sys)//"_"//trim(method)//"_englevel_T"//trim(Tch)//"_mixexcf1_rd"//trim(rdch)//".txt")
open(279,file=""//trim(sys)//"_"//trim(method)//"_englevel_T"//trim(Tch)//"_mixexcf2_rd"//trim(rdch)//".txt")
call CkqJ(J1,CJst,JM)
! set coffmag
call tlm(coffmag,0.0d0,0.0d0)   
coffmix1=0.0d0
coffmix2=0.0d0
do l=1,6
   coffmix1(l,l-1)=(2.0d0*J1+dble(l)+1.0d0)/2.0d0
   coffmix1(l,l+1)=-2.0d0/(2.0d0*J1+dble(l)+2.0d0)
   do j=-1,1,2
      coffmix2(l,l-1+j)=coffmix2(l,l-1+j)&
           -coffmix1(l,l+j)*(dble(l+j)*(2.0d0*J1-dble(l+j)+1.0d0)*(2.0d0*J1+dble(l+j)+1.0d0)/(4.0d0*(2.0d0*dble(l+j)+1.0d0)))
      coffmix2(l,l+1+j)=coffmix2(l,l+1+j)&
           -coffmix1(l,l+j)*(dble(l+j)+1.0d0)/(2.0d0*dble(l+j)+1.0d0)
   end do
end do
coffmix2=coffmix2*(Dex/Dso)*(L1+S1+1.0d0)/((J1+2.0d0)*S1)

do ii=1,int(2.0d0*J1)+1    
   eigen0(ii)=-2.0d0*S1/(J1+1.0d0)*Hex(1)*JM(ii,2) ! Hex <0, g(1)-1 <0, <S> >0 
   do l=2,6,2
      do m=-l,l
         englvcf(ii)=englvcf(ii)+coffmag(l,m)*Alm(l,m,1)*Stevens(l)*CJst(ii,l,0)!&
      end do
   end do
   do ll=-1,8
      abc=-Dex**2/Dso*(L1+1.0d0)/(3.0d0*S1)
      englvmixex1(ii)=englvmixex1(ii)+abc*coffmix1(1,ll)*CJst(ii,ll,0)
      englvmixex2(ii)=englvmixex2(ii)+abc*coffmix2(1,ll)*(-CJst(ii,ll,0))
      do l=2,6,2 ! l for Alm
         do m=-l,l
            abc=xi(l)*dble(l*(l+1))/(dble(2*l+1))*Dex/Dso*coffmag(l,m)*Alm(l,m,1)
            englvmixcf1(ii)=englvmixcf1(ii)+abc*coffmix1(l,ll)*(-CJst(ii,ll,0))
            englvmixcf2(ii)=englvmixcf2(ii)+abc*coffmix2(l,ll)*CJst(ii,ll,0)
         end do
      end do
   end do

   write(281,*) ii,0.0d0,0
   write(281,*) ii,0.0d0,1
   write(275,*) ii,0.0d0+eigen0(ii),0
   write(275,*) ii,0.0d0+eigen0(ii),1
   write(280,*) ii,0.0d0+eigen0(ii)+englvcf(ii),0
   write(280,*) ii,0.0d0+eigen0(ii)+englvcf(ii),1
   
   write(276,*) ii,0.0d0+eigen0(ii)+englvmixex1(ii),0
   write(276,*) ii,0.0d0+eigen0(ii)+englvmixex1(ii),1
   write(278,*) ii,0.0d0+eigen0(ii)+englvcf(ii)+englvmixex1(ii)+englvmixcf1(ii),0
   write(278,*) ii,0.0d0+eigen0(ii)+englvcf(ii)+englvmixex1(ii)+englvmixcf1(ii),1
   
   write(277,*) ii,0.0d0+eigen0(ii)+englvmixex1(ii)+englvmixex2(ii),0
   write(277,*) ii,0.0d0+eigen0(ii)+englvmixex1(ii)+englvmixex2(ii),1
   write(279,*) ii,0.0d0+eigen0(ii)+englvcf(ii)+englvmixex1(ii)+englvmixex2(ii)+&
        englvmixcf1(ii)+englvmixcf2(ii),0
   write(279,*) ii,0.0d0+eigen0(ii)+englvcf(ii)+englvmixex1(ii)+englvmixex2(ii)+&
        englvmixcf1(ii)+englvmixcf2(ii),1
   write(275,*)
   write(280,*)
   write(281,*)
   write(276,*)
   write(277,*)
   write(278,*)
   write(279,*)
end do ! do ii

end subroutine energy_analytical
