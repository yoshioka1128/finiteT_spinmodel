program test
  implicit none
  integer,parameter :: jigen=26
  double precision,parameter :: pi=4*atan(1d0) 
  character :: jobz,uplo
  integer :: n,lda,lwork,info
  double complex :: matrix(jigen,jigen),CEF(jigen,jigen)
  double complex :: work(3*jigen-1)
  double precision :: Eigenvalue(jigen)
  double precision :: rwork(3*jigen)
  !zheevの初期設定
  real(8) :: w3jsym,b
  real(8) :: a(0:3) !対角成分の結晶場項
  integer(4) ::i,j,k,m
  double complex :: A20,A40,A60,A2_2,A4_2,A44,A6_2,A64,A6_6
  double complex :: Le(0:180),Se(0:180),sxe(0:180),lex(0:180),sye(0:180),ley(0:180),Sc(0:180) 
  double precision :: reduce2,reduce4,reduce6 !還元行列要素
  double precision :: T,Tc,lambda,Hm,Hm0,grandenergy,k1,yk1,mt
  double precision :: theta,phi,Z,Nw(jigen)
  integer :: lz(jigen),angle_easy
  double precision :: sz(jigen)
  complex*16,parameter :: imag=cmplx(0d0,1d0)
  double precision :: Energy(0:180),FreeEnergy(0:180),N4f(0:180)
  double precision :: V,E_d,F_easy
  jobz='V'; uplo='U'


  n=jigen
  lda=jigen
  lwork=3*jigen-1

  lz(1)=3d0; lz(2)=3d0; lz(3)=2d0; lz(4)=2d0; lz(5)=1d0; lz(6)=1d0; lz(7)=0d0
  lz(8)=0d0; lz(9)=-1d0; lz(10)=-1d0; lz(11)=-2d0; lz(12)=-2d0; lz(13)=-3d0; lz(14)=-3d0

  do i=1,jigen,2
     sz(i)=0.5
  end do
  do i=2,jigen,2
     sz(i)=-0.5
  end do






  A20=-466.43   
  A40=-22.67
  A60=5.79 
  Hm0=429.9




  lambda=916.0
  Tc=653.0
  T=0


!!!!!!!!!!!!!!!!!! parameter !!!!!!!!!!!!!!
  
  V=-2300
  E_d=-2200

  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  reduce2=-2d0*sqrt(7d0/15d0)
  reduce4=sqrt(14d0/11d0)
  reduce6=(-10d0*sqrt(7d0/(3d0*11d0*13d0)))

  !結晶場項の対角要素の計算


  !結晶場項の非対角要素の計算
  do i=1,jigen
     do j=1,jigen
        CEF(i,j)=0
     end do
  end do

  !A20の寄与
  do i=1,14
     CEF(i,i)=CEF(i,i)+2d0*reduce2*w3jsym(3d0,2d0,3d0,-1d0*lz(i),0d0)*A20*((-1.0)**(lz(i)-3d0))
  end do
  !A40の寄与
  do i=1,14
     CEF(i,i)=CEF(i,i)+8d0*reduce4*w3jsym(3d0,4d0,3d0,-1d0*lz(i),0d0)*A40*((-1.0)**(lz(i)-3d0))
  end do
  !A60の寄与
  do i=1,14
     CEF(i,i)=CEF(i,i)+16d0*reduce6*w3jsym(3d0,6d0,3d0,-1d0*lz(i),0d0)*A60*((-1.0)**(lz(i)-3d0))
  end do




  ! print*,a(0),a(1),a(2),a(3)
  print *, '#Temperature(K),K_1(MJ/m^3)'
  phi=pi/12.0*0d0
  do m=0,32
     T=m*(20.0)
     Hm=Hm0*(1-0.7*(T/Tc)**(3.0/2.0)-(1-0.7)*(T/Tc)**(5.0/2.0))**(1.0/3.0)

     F_easy=10000
     angle_easy=1000
     do k=0,180
        theta=-0.5*pi+k*pi/180.0
        do i=1,jigen
           do j=1,jigen
              matrix(i,j)=CEF(i,j) !行列要素に結晶場項のみ代入
           end do
        end do

        do i=1,jigen,2
           matrix(i,i)=matrix(i,i)+Hm*cos(theta) !交換項
        end do

        do i=2,jigen,2
           matrix(i,i)=matrix(i,i)-Hm*cos(theta) !交換項
        end do

        do i=1,jigen
           matrix(i,i)=matrix(i,i)+lambda*lz(i)*sz(i) !スピン軌道相互作用の対角項
        end do

        matrix(2,3)=lambda*sqrt(6.0)/2.0; matrix(12,13)=matrix(2,3)
        matrix(4,5)=lambda*sqrt(10.0)/2.0; matrix(10,11)=matrix(4,5)
        matrix(6,7)=lambda*sqrt(12.0)/2.0; matrix(8,9)=matrix(6,7) !スピン軌道相互作用の非対角項

        do i=1,13,2
           matrix(i,i+1)=Hm*sin(theta)*cos(phi)-Hm*imag*sin(theta)*sin(phi)
        end do


        call ligand(matrix,15,V,E_d,pi/2d0,0d0)
        call ligand(matrix,17,V,E_d,pi/2d0,pi/3d0)
        call ligand(matrix,19,V,E_d,pi/2d0,pi/3d0*2d0)
        call ligand(matrix,21,V,E_d,pi/2d0,pi)
        call ligand(matrix,23,V,E_d,pi/2d0,pi/3d0*4d0)
        call ligand(matrix,25,V,E_d,pi/2d0,pi/3d0*5d0)

        do i=15,jigen

           matrix(i,i)=E_d
        end do



        !           if((k==90).and.m==1) then
        !           print *,"A="
        !           do i=1,jigen
        !              write(6,650) (abs(matrix(i,j)),j=1,jigen)
        !650           format(26e12.4)
        !              end do
        !                     end if


        call zheev(jobz,uplo,n,matrix,lda,Eigenvalue,work,lwork,rwork,info)

     !              if((k==90).and.m==100) then
     !              print *,"V="
     !              do i=1,14
     !                 write(6,630) ((matrix(i,j)),j=1,1)
     !   630           format(26e12.4)
     !                               end do

     !                end if


        !以降解析


        grandenergy=Eigenvalue(1)
        do i=1,jigen
           Eigenvalue(i)=Eigenvalue(i)-grandenergy
        end do

        Z=0
        do i=1,jigen
           Z=Z+exp(-Eigenvalue(i)/T)
        end do

        if(T>0.1) then
           do i=1,jigen
              Nw(i)=exp(-Eigenvalue(i)/T)/Z
           end do
        else
           Nw(1)=1
           do i=2,jigen
              Nw(i)=0
           end do
        end if
        Energy(k)=0
        do i=1,jigen
           Energy(k)=Energy(k)+Eigenvalue(i)*Nw(i)
        end do
        Energy(k)=Energy(k)+grandenergy

        if(T<0.1) then
           FreeEnergy(k)=grandenergy
        else
           FreeEnergy(k)=-T*log(Z)+grandenergy
        end if

 le(k)=0
        lex(k)=0
        ley(k)=0
        Se(k)=0
        sye(k)=0
        sxe(k)=0
        Sc(k)=0

        do j=1,jigen
           Le(k)=Le(k)+(abs(matrix(1,j))**2)*3.0*Nw(j)
           Le(k)=Le(k)+(abs(matrix(2,j))**2)*3.0*Nw(j)
           Le(k)=Le(k)+(abs(matrix(3,j))**2)*2.0*Nw(j)
           Le(k)=Le(k)+(abs(matrix(4,j))**2)*2.0*Nw(j)
           Le(k)=Le(k)+(abs(matrix(5,j))**2)*1.0*Nw(j)
           Le(k)=Le(k)+(abs(matrix(6,j))**2)*1.0*Nw(j)

           Le(k)=Le(k)-(abs(matrix(9,j))**2)*1.0*Nw(j)
           Le(k)=Le(k)-(abs(matrix(10,j))**2)*(1.0)*Nw(j)   
           Le(k)=Le(k)-(abs(matrix(11,j))**2)*(2.0)*Nw(j)
           Le(k)=Le(k)-(abs(matrix(12,j))**2)*(2.0)*Nw(j)
           Le(k)=Le(k)-(abs(matrix(13,j))**2)*(3.0)*Nw(j)
           Le(k)=Le(k)-(abs(matrix(14,j))**2)*(3.0)*Nw(j)




           do i=1,14,2
              Se(k)=Se(k)+0.5*(abs(matrix(i,j))**2)*Nw(j)
           end do
           do i=2,14,2
              Se(k)=Se(k)-0.5*(abs(matrix(i,j)**2))*Nw(j)
           end do

           do i=15,jigen,2
              Sc(k)=Sc(k)+0.5*(abs(matrix(i,j))**2)*Nw(j)
           end do

           do i=16,jigen,2
              Sc(k)=Sc(k)-0.5*(abs(matrix(i,j))**2)*Nw(j)
           end do

           Lex(k)=Lex(k)+sqrt(6.0)*((matrix(1,j)*conjg(matrix(3,j)))+(conjg(matrix(1,j))*matrix(3,j)))*Nw(j)*0.5
           Lex(k)=Lex(k)+sqrt(10.0)*((matrix(3,j)*conjg(matrix(5,j)))+(conjg(matrix(3,j))*matrix(5,j)))*Nw(j)*0.5
           Lex(k)=Lex(k)+sqrt(12.0)*((matrix(5,j)*conjg(matrix(7,j)))+(conjg(matrix(5,j))*matrix(7,j)))*Nw(j)*0.5
           Lex(k)=Lex(k)+sqrt(12.0)*((matrix(7,j)*conjg(matrix(9,j)))+(conjg(matrix(7,j))*matrix(9,j)))*Nw(j)*0.5
           Lex(k)=Lex(k)+sqrt(10.0)*((matrix(9,j)*conjg(matrix(11,j)))+(conjg(matrix(9,j))*matrix(11,j)))*Nw(j)*0.5
           Lex(k)=Lex(k)+sqrt(6.0)*((matrix(11,j)*conjg(matrix(13,j)))+(conjg(matrix(11,j))*matrix(13,j)))*Nw(j)*0.5

           Lex(k)=Lex(k)+sqrt(6.0)*((matrix(2,j)*conjg(matrix(4,j)))+(conjg(matrix(2,j))*matrix(4,j)))*Nw(j)*0.5
           Lex(k)=Lex(k)+sqrt(10.0)*((matrix(4,j)*conjg(matrix(6,j)))+(conjg(matrix(4,j))*matrix(6,j)))*Nw(j)*0.5
           Lex(k)=Lex(k)+sqrt(12.0)*((matrix(6,j)*conjg(matrix(8,j)))+(conjg(matrix(6,j))*matrix(8,j)))*Nw(j)*0.5
           Lex(k)=Lex(k)+sqrt(12.0)*((matrix(8,j)*conjg(matrix(10,j)))+(conjg(matrix(8,j))*matrix(10,j)))*Nw(j)*0.5
           Lex(k)=Lex(k)+sqrt(10.0)*((matrix(10,j)*conjg(matrix(12,j)))+(conjg(matrix(10,j))*matrix(12,j)))*Nw(j)*0.5
           Lex(k)=Lex(k)+sqrt(6.0)*((matrix(12,j)*conjg(matrix(14,j)))+(conjg(matrix(12,j))*matrix(14,j)))*Nw(j)*0.5


           Ley(k)=Ley(k)+sqrt(6.0)*((matrix(1,j)*conjg(matrix(3,j)))-(conjg(matrix(1,j))*matrix(3,j)))*Nw(j)*0.5
           Ley(k)=Ley(k)+sqrt(10.0)*((matrix(3,j)*conjg(matrix(5,j)))-(conjg(matrix(3,j))*matrix(5,j)))*Nw(j)*0.5
           Ley(k)=Ley(k)+sqrt(12.0)*((matrix(5,j)*conjg(matrix(7,j)))-(conjg(matrix(5,j))*matrix(7,j)))*Nw(j)*0.5
           Ley(k)=Ley(k)+sqrt(12.0)*((matrix(7,j)*conjg(matrix(9,j)))-(conjg(matrix(7,j))*matrix(9,j)))*Nw(j)*0.5
           Ley(k)=Ley(k)+sqrt(10.0)*((matrix(9,j)*conjg(matrix(11,j)))-(conjg(matrix(9,j))*matrix(11,j)))*Nw(j)*0.5
           Ley(k)=Ley(k)+sqrt(6.0)*((matrix(11,j)*conjg(matrix(13,j)))-(conjg(matrix(11,j))*matrix(13,j)))*Nw(j)*0.5

           Ley(k)=Ley(k)+sqrt(6.0)*((matrix(2,j)*conjg(matrix(4,j)))-(conjg(matrix(2,j))*matrix(4,j)))*Nw(j)*0.5
           Ley(k)=Ley(k)+sqrt(10.0)*((matrix(4,j)*conjg(matrix(6,j)))-(conjg(matrix(4,j))*matrix(6,j)))*Nw(j)*0.5
           Ley(k)=Ley(k)+sqrt(12.0)*((matrix(6,j)*conjg(matrix(8,j)))-(conjg(matrix(6,j))*matrix(8,j)))*Nw(j)*0.5
           Ley(k)=Ley(k)+sqrt(12.0)*((matrix(8,j)*conjg(matrix(10,j)))-(conjg(matrix(8,j))*matrix(10,j)))*Nw(j)*0.5
           Ley(k)=Ley(k)+sqrt(10.0)*((matrix(10,j)*conjg(matrix(12,j)))-(conjg(matrix(10,j))*matrix(12,j)))*Nw(j)*0.5
           Ley(k)=Ley(k)+sqrt(6.0)*((matrix(12,j)*conjg(matrix(14,j)))-(conjg(matrix(12,j))*matrix(14,j)))*Nw(j)*0.5



           Sxe(k)=Sxe(k)+((matrix(1,j)*conjg(matrix(2,j)))+(conjg(matrix(1,j))*matrix(2,j)))*Nw(j)*0.5
           Sxe(k)=Sxe(k)+((matrix(3,j)*conjg(matrix(4,j)))+(conjg(matrix(3,j))*matrix(4,j)))*Nw(j)*0.5
           Sxe(k)=Sxe(k)+((matrix(5,j)*conjg(matrix(6,j)))+(conjg(matrix(5,j))*matrix(6,j)))*Nw(j)*0.5
           Sxe(k)=Sxe(k)+((matrix(7,j)*conjg(matrix(8,j)))+(conjg(matrix(7,j))*matrix(8,j)))*Nw(j)*0.5
           Sxe(k)=Sxe(k)+((matrix(9,j)*conjg(matrix(10,j)))+(conjg(matrix(9,j))*matrix(10,j)))*Nw(j)*0.5
           Sxe(k)=Sxe(k)+((matrix(11,j)*conjg(matrix(12,j)))+(conjg(matrix(11,j))*matrix(12,j)))*Nw(j)*0.5
           Sxe(k)=Sxe(k)+((matrix(13,j)*conjg(matrix(14,j)))+(conjg(matrix(13,j))*matrix(14,j)))*Nw(j)*0.5


           Sye(k)=Sye(k)+((matrix(1,j)*conjg(matrix(2,j)))-(conjg(matrix(1,j))*matrix(2,j)))*Nw(j)*0.5
           Sye(k)=Sye(k)+((matrix(3,j)*conjg(matrix(4,j)))-(conjg(matrix(3,j))*matrix(4,j)))*Nw(j)*0.5
           Sye(k)=Sye(k)+((matrix(5,j)*conjg(matrix(6,j)))-(conjg(matrix(5,j))*matrix(6,j)))*Nw(j)*0.5
           Sye(k)=Sye(k)+((matrix(7,j)*conjg(matrix(8,j)))-(conjg(matrix(7,j))*matrix(8,j)))*Nw(j)*0.5
           Sye(k)=Sye(k)+((matrix(9,j)*conjg(matrix(10,j)))-(conjg(matrix(9,j))*matrix(10,j)))*Nw(j)*0.5
           Sye(k)=Sye(k)+((matrix(11,j)*conjg(matrix(12,j)))-(conjg(matrix(11,j))*matrix(12,j)))*Nw(j)*0.5
           Sye(k)=Sye(k)+((matrix(13,j)*conjg(matrix(14,j)))-(conjg(matrix(13,j))*matrix(14,j)))*Nw(j)*0.5

           

        end do

        N4f(k)=0
        do i=1,jigen
           do j=1,14
              N4f(k)=N4f(k)+(abs(matrix(j,i))**2d0)*Nw(i)
           end do
           
        end do
        mt=((1-0.7*(T/Tc)**(1.5)-0.3*(T/Tc)**(2.5))**(1.0/3.0))
     yk1=47.07*(mt**3.0)+8.0/7.0*21.77*((mt**3.0)-(mt**10.0))+8.0/7.0*(-81.84)* &
          ((mt**3.0)-18.0/11.0*(mt**10.0)+7.0/11.0*(mt**21.0))
        if(FreeEnergy(k)+yk1*(sin(theta)**2)<F_easy) then
           F_easy=FreeEnergy(k)+yk1*(sin(theta)**2)
           angle_easy=k
           end if


        
        !        print *,theta,Eigenvalue(1)+grandenergy,Eigenvalue(2)+grandenergy,FreeEnergy(k)
        end do

     mt=((1-0.7*(T/Tc)**(1.5)-0.3*(T/Tc)**(2.5))**(1.0/3.0))
     yk1=47.07*(mt**3.0)+8.0/7.0*21.77*((mt**3.0)-(mt**10.0))+8.0/7.0*(-81.84)* &
          ((mt**3.0)-18.0/11.0*(mt**10.0)+7.0/11.0*(mt**21.0))
     k1=(FreeEnergy(92)+FreeEnergy(88)-(2d0*FreeEnergy(90)))/(pi/90.0)/(pi/90.0)*0.5

     print *,T,(k1)*10d0/62d0*1.017,yk1*10d0/62d0*1.017,real(Se(angle_easy)),real(Le(angle_easy)),real(Sc(90)) &
          ,8.327*mt
     
  end do






end program test

subroutine ligand(a,retu,V,E_d,theta_ligand,phi_ligand)
  implicit none
  integer :: i,j,k,retu
  integer,parameter :: jigen=26
  double complex :: a(jigen,jigen)
  double precision :: V,theta_ligand,phi_ligand,E_d
  complex*16,parameter :: imag=cmplx(0d0,1d0)
  a(1,retu)=-V*sqrt(5.0/16.0)*(sin(theta_ligand)**3.0)*(cos(3.0*phi_ligand)-imag*sin(3.0*phi_ligand)) !共役なので-imagとなっている
  a(3,retu)=V*sqrt(15.0/8.0)*(sin(theta_ligand)**2.0)*cos(theta_ligand)*(cos(2.0*phi_ligand)-imag*sin(2.0*phi_ligand))
  a(5,retu)=-V*sqrt(3.0/16.0)*sin(theta_ligand)*(5*cos(theta_ligand)*cos(theta_ligand)-1)*(cos(phi_ligand)-imag*sin(phi_ligand))
  a(7,retu)=V*sqrt(1.0/4.0)*(5.0*cos(theta_ligand)*cos(theta_ligand)*cos(theta_ligand)-3.0*cos(theta_ligand))
  a(9,retu)=V*sqrt(3.0/16.0)*sin(theta_ligand)*(5*cos(theta_ligand)*cos(theta_ligand)-1)*(cos(phi_ligand)+imag*sin(phi_ligand))
  a(11,retu)=V*sqrt(15.0/8.0)*(sin(theta_ligand)**2.0)*cos(theta_ligand)*(cos(2.0*phi_ligand)+imag*sin(2.0*phi_ligand))
  a(13,retu)=V*sqrt(5.0/16.0)*(sin(theta_ligand)**3.0)*(cos(3.0*phi_ligand)+imag*sin(3.0*phi_ligand))

  do i=2,14,2
     a(i,retu+1)=a(i-1,retu)
  end do
  a(retu,retu)=E_d
  a(retu+1,retu+1)=E_d
end subroutine ligand
