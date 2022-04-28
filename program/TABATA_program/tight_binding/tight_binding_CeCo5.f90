program test
  implicit none
  integer,parameter :: jigen=20,Nmax=100,emax=100,loop_max=500
  double precision,parameter :: pi=4*atan(1d0) 
  character :: jobz,uplo
  integer :: n,lda,lwork,info
  double complex :: matrix(jigen,jigen),CEF(jigen,jigen)
  double complex :: work(3*jigen-1),f2f1(14,14),f2f1_pr(14,14),f2f1_loop1(14,14)
  double precision :: Eigenvalue(jigen)
  double precision :: rwork(3*jigen)
  !zheevの初期設定
  real(8) :: w3jsym,b
  real(8) :: a(0:3) !対角成分の結晶場項
  integer(4) ::i,j,k,m,l1,l2,q,angle(12)
  double complex :: A20,A40,A60,A2_2,A4_2,A44,A6_2,A64,A6_6,y3m
  double complex :: lx_ex,ly_ex,lz_ex,sx_ex,sy_ex,sz_ex
  double precision :: reduce2,reduce4,reduce6 !還元行列要素
  double precision :: T,Tc,lambda,Hm,Hm0,grandenergy,k1,yk1,mt
  double precision :: theta,phi,Z,Nw(jigen)
  integer :: lz(jigen),nx,ny,theta_int,loop,loop2,count
  double precision :: sz(jigen),mesh=10000d0/emax,integral_d,integral_f
  double precision :: integral_ddn,integral_dup,integral_fdn,integral_fup
  complex*16,parameter :: imag=cmplx(0d0,1d0)
  double precision :: Energy(0:180),FreeEnergy(0:180)
  double precision :: Lz_expect(0:180),Sz_expect(0:180)
  double precision :: V,E_d,eps
  double precision :: r1x,r2x,r3x,r1y,r2y,r3y,kx,ky,ek1,ek2,tab,tac,tbc,U,mixing
  double precision :: E(-50*emax:50*emax),rho(-50*emax:50*emax),Prho(jigen,-50*emax:50*emax)
  double precision :: rhofup(-50*emax:50*emax),rhofdn(-50*emax:50*emax),rhodup(-50*emax:50*emax),rhoddn(-50*emax:50*emax)
  double precision :: all_eigenvalue(Nmax*Nmax*jigen),mu(0:loop_max),N4f(0:loop_max),Nm(jigen),mu_min,mu_plu,delta,delta_matrix
  double precision :: J_partition,Z_partition,eps_matrix,integral
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


  r1x=1d0
  r2x=-1d0/2d0
  r3x=-1d0/2d0
  r1y=0
  r2y=sqrt(3d0)/2d0
  r3y=-sqrt(3d0)/2d0




  !  A20=450.638563
  !  A40=-62.81006927
  !  A60=-3.354432565
  !  Hm0=490.5

  A20=352.761407771946*0    
  A2_2=848.463507683960*0 !    
  A40=-57.5547507771566*0    
  A4_2=-107.848645323908*0  !   
  A44=-153.101109490958*0  !   
  A60=-3.77177119406693*0     
  A6_2=7.37649603613955*0  !   
  A64=-37.4610114035808*0   !  
  A6_6=-6.48478297211961*0   !  
  Hm0=474.2*0



  A20=-466.43   
  A40=-22.67
  A60=5.79 
  Hm0=429.9




  lambda=916.0
  Tc=653d0


 

 

  reduce2=-2d0*sqrt(7d0/15d0)
  reduce4=sqrt(14d0/11d0)
  reduce6=(-10d0*sqrt(7d0/(3d0*11d0*13d0)))

  !結晶場項の対角要素の計算

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



  

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! main parameter !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  tbc=-10000 !cf DOS from ab intio
  V=-5000d0*0d0
  phi=pi/12d0
  loop=0
  loop2=0
  U=70000 !cf Saso

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




  
  
!  eps=Nmax**2d0/1000
  eps=0.1d0
  eps_matrix=0.0001d0
open(9,file='MAE.dat',status='replace')
 
angle(1)=90
angle(2)=92
angle(3)=100
angle(4)=110
angle(5)=120
angle(6)=130
angle(7)=140
angle(8)=150
angle(9)=160
angle(10)=170
angle(11)=180
angle(12)=190

do q=1,12,1
     theta_int=angle(q)
     print*,theta_int-90,"loop start"
     loop=0
     loop2=0
     T=10d0
     N4f(:)=0
     E_d=-10000
     mu_min=8000d0
     mu_plu=20000d0
     delta=50000
     mixing=0.5d0
     count=0.0
     f2f1(:,:)=0
     f2f1_pr(:,:)=0
     f2f1_loop1(:,:)=0
     integral=0     

  do while(abs(integral-Nmax*Nmax*4d0)>eps.and.loop2<200)

     if(loop2>99) then
        eps=0.0001
        end if
     loop2=loop2+1
     mu(:)=(mu_min+mu_plu)/2d0
     loop=0
     N4f(:)=0
     delta=10000
     delta_matrix=10000

!!!!!!!!!!!!!!!!!!恣意的な初期条件!!!!!!!!!!!!!!!!!  
!     f2f1(:,:)=0
!     f2f1_pr(:,:)=0

     !     do i=7,8,1
!        f2f1(i,i)=1
!     end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     do while(loop<loop_max.and.(delta_matrix>eps_matrix))
        loop=loop+1

        all_eigenvalue(:)=0

        f2f1_Pr(:,:)=f2f1(:,:)*mixing+f2f1_pr(:,:)*(1d0-mixing)
        f2f1(:,:)=0
        Nm(:)=0
        lx_ex=0
        ly_ex=0
        lz_ex=0
        sx_ex=0
        sy_ex=0
        sz_ex=0

      
        integral=0



        do i=-50*emax,50*emax
           E(i)=i*mesh
           rho(i)=0
           do j=1,jigen
              Prho(j,i)=0
           end do
           rhofup(i)=0
           rhofdn(i)=0
           rhodup(i)=0
           rhoddn(i)=0
        end do



        do nx=1,Nmax
           do ny=1,Nmax


              kx=4d0*pi/3d0/Nmax*(nx-0.5*ny)
              ky=4d0*pi/3d0/Nmax*sqrt(3d0)/2d0*ny







!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!局在項!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

              
              theta=-0.5*pi+theta_int*pi/180.0


              Hm=Hm0*(1-0.7*(T/Tc)**(3.0/2.0)-(1-0.7)*(T/Tc)**(5.0/2.0))**(1.0/3.0)


              matrix(:,:)=CEF(:,:)



              do i=1,14,2
                 matrix(i,i)=matrix(i,i)+Hm*cos(theta) !交換項
              end do

              do i=2,14,2
                 matrix(i,i)=matrix(i,i)-Hm*cos(theta) !交換項
              end do

              do i=1,14
                 matrix(i,i)=matrix(i,i)+lambda*lz(i)*sz(i) !スピン軌道相互作用の対角項
              end do

              matrix(2,3)=lambda*sqrt(6.0)/2.0; matrix(12,13)=matrix(2,3)
              matrix(4,5)=lambda*sqrt(10.0)/2.0; matrix(10,11)=matrix(4,5)
              matrix(6,7)=lambda*sqrt(12.0)/2.0; matrix(8,9)=matrix(6,7) !スピン軌道相互作用の非対角項

              do i=1,13,2
                 matrix(i,i+1)=Hm*sin(theta)*cos(phi)-Hm*imag*sin(theta)*sin(phi)
              end do



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!局在項!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!相互作用項!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

              do i=1,14
                 do j=1,14

                    matrix(i,i)=matrix(i,i)+f2f1_pr(j,j)*U  !Hartree term

                 end do
              end do

              do i=1,14
                 do j=1,14

                    matrix(i,j)=matrix(i,j)-f2f1_pr(j,i)*U   !Fock term

                 end do
              end do


       


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!相互作用項!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!混成項!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


              do i=1,14,2
                 matrix(i,15)=V*y3m(lz(i),pi/2d0,0)*(cos(kx*r1x+ky*r1y)+imag*sin(kx*r1x+ky*r1y)) &
                      +V*y3m(lz(i),pi/2d0,2d0*pi/3d0)*(cos(kx*r2x+ky*r2y)+imag*sin(kx*r2x+ky*r2y)) &
                      +V*y3m(lz(i),pi/2d0,4d0*pi/3d0)*(cos(kx*r3x+ky*r3y)+imag*sin(kx*r3x+ky*r3y))  !exp(-ik(A-B))=exp(ik(B-A))
                 matrix(i,17)=V*y3m(lz(i),pi/2d0,0)*(cos(kx*r1x+ky*r1y)-imag*sin(kx*r1x+ky*r1y)) &
                      +V*y3m(lz(i),pi/2d0,2d0*pi/3d0)*(cos(kx*r2x+ky*r2y)-imag*sin(kx*r2x+ky*r2y)) &
                      +V*y3m(lz(i),pi/2d0,4d0*pi/3d0)*(cos(kx*r3x+ky*r3y)-imag*sin(kx*r3x+ky*r3y))  !exp(-ik(A-C))=exp(-ik(B-A)) 
              end do


              do i=2,14,2
                 matrix(i,16)=V*y3m(lz(i),pi/2d0,0)*(cos(kx*r1x+ky*r1y)+imag*sin(kx*r1x+ky*r1y)) &
                      +V*y3m(lz(i),pi/2d0,2d0*pi/3d0)*(cos(kx*r2x+ky*r2y)+imag*sin(kx*r2x+ky*r2y)) &
                      +V*y3m(lz(i),pi/2d0,4d0*pi/3d0)*(cos(kx*r3x+ky*r3y)+imag*sin(kx*r3x+ky*r3y))
                 matrix(i,18)=V*y3m(lz(i),pi/2d0,0)*(cos(kx*r1x+ky*r1y)-imag*sin(kx*r1x+ky*r1y)) &
                      +V*y3m(lz(i),pi/2d0,2d0*pi/3d0)*(cos(kx*r2x+ky*r2y)-imag*sin(kx*r2x+ky*r2y)) &
                      +V*y3m(lz(i),pi/2d0,4d0*pi/3d0)*(cos(kx*r3x+ky*r3y)-imag*sin(kx*r3x+ky*r3y))
              end do

              matrix(15,17)=tbc*((cos(kx*r1x+ky*r1y)+imag*sin(kx*r1x+ky*r1y))+(cos(kx*r2x+ky*r2y)+imag*sin(kx*r2x+ky*r2y)) &
                   +(cos(kx*r3x+ky*r3y)+imag*sin(kx*r3x+ky*r3y)) )   !exp(-ik(B-C))
              matrix(16,18)=tbc*((cos(kx*r1x+ky*r1y)+imag*sin(kx*r1x+ky*r1y))+(cos(kx*r2x+ky*r2y)+imag*sin(kx*r2x+ky*r2y)) &
                   +(cos(kx*r3x+ky*r3y)+imag*sin(kx*r3x+ky*r3y)) )   !exp(-ik(B-C))
              matrix(15,19)=conjg(matrix(15,17)) !exp(-ik(B-A))
              matrix(16,20)=conjg(matrix(15,17)) !exp(-ik(B-A))
              matrix(17,19)=matrix(15,17)        !exp(-ik(C-A))=exp(-ik(B-C))
              matrix(18,20)=matrix(16,18)        !exp(-ik(C-A))=exp(-ik(B-C))   


              do i=1,14
                 matrix(i,i)=matrix(i,i)-E_d
              end do


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!         




              if((theta_int==900).and.nx==Nmax.and.ny==Nmax.and.loop==1.and.loop2==1) then
                 print* ,"kx=",kx,"ky=",ky
                 print *,"A="
                 do i=1,jigen
                    write(6,620) (abs(matrix(i,j)),j=1,jigen)
620                 format(26e12.4)
                 end do
              end if

              call zheev(jobz,uplo,n,matrix,lda,Eigenvalue,work,lwork,rwork,info)

              
              do i=1,jigen
                 all_eigenvalue((Nmax)*jigen*(nx-1)+jigen*(ny-1)+i)=eigenvalue(i)
              end do




              if((theta_int==900).and.nx<2.and.ny<2.and.loop==1.and.loop2==1) then
                 print *,"V="
                 do i=1,jigen
                    write(6,630) (abs(matrix(i,j)),j=1,jigen)
630                 format(26e12.4)
                 end do
                 print*,"Eigenvalue="
                 write(6,640) (Eigenvalue(i),i=1,jigen)
640              format(26e12.4)
              end if



             do i=1,jigen
                 do j=1,14
                    N4f(loop)=N4f(loop)+(abs(matrix(j,i))**2d0/Nmax/Nmax)/(1d0+exp((eigenvalue(i)-mu(loop))/T))
                    end do
do j=1,jigen
    Nm(j)=Nm(j)+(abs(matrix(j,i))**2d0/Nmax/Nmax)/(1d0+exp((eigenvalue(i)-mu(loop))/T))
            end do    
              end do

        do j=1,jigen
           lz_ex=lz_ex+(abs(matrix(1,j))**2)*3.0*1d0/(1d0+exp((eigenvalue(j)-mu(loop))/T))
           lz_ex=lz_ex+(abs(matrix(2,j))**2)*3.0*1d0/(1d0+exp((eigenvalue(j)-mu(loop))/T))
           lz_ex=lz_ex+(abs(matrix(3,j))**2)*2.0*1d0/(1d0+exp((eigenvalue(j)-mu(loop))/T))
           lz_ex=lz_ex+(abs(matrix(4,j))**2)*2.0*1d0/(1d0+exp((eigenvalue(j)-mu(loop))/T))
           lz_ex=lz_ex+(abs(matrix(5,j))**2)*1.0*1d0/(1d0+exp((eigenvalue(j)-mu(loop))/T))
           lz_ex=lz_ex+(abs(matrix(6,j))**2)*1.0*1d0/(1d0+exp((eigenvalue(j)-mu(loop))/T))

           lz_ex=lz_ex-(abs(matrix(9,j))**2)*1.0*1d0/(1d0+exp((eigenvalue(j)-mu(loop))/T))
           lz_ex=lz_ex-(abs(matrix(10,j))**2)*(1.0)*1d0/(1d0+exp((eigenvalue(j)-mu(loop))/T))   
           lz_ex=lz_ex-(abs(matrix(11,j))**2)*(2.0)*1d0/(1d0+exp((eigenvalue(j)-mu(loop))/T))
           lz_ex=lz_ex-(abs(matrix(12,j))**2)*(2.0)*1d0/(1d0+exp((eigenvalue(j)-mu(loop))/T))
           lz_ex=lz_ex-(abs(matrix(13,j))**2)*(3.0)*1d0/(1d0+exp((eigenvalue(j)-mu(loop))/T))
           lz_ex=lz_ex-(abs(matrix(14,j))**2)*(3.0)*1d0/(1d0+exp((eigenvalue(j)-mu(loop))/T))




           do i=1,14,2
              sz_ex=sz_ex+0.5*(abs(matrix(i,j))**2)*1d0/(1d0+exp((eigenvalue(j)-mu(loop))/T))
           end do
           do i=2,14,2
              sz_ex=sz_ex-0.5*(abs(matrix(i,j)**2))*1d0/(1d0+exp((eigenvalue(j)-mu(loop))/T))
           end do
           
           lx_ex=lx_ex+sqrt(6.0)*((matrix(1,j)*conjg(matrix(3,j)))+(conjg(matrix(1,j)) &
                *matrix(3,j)))*1d0/(1d0+exp((eigenvalue(j)-mu(loop))/T))*0.5
           lx_ex=lx_ex+sqrt(10.0)*((matrix(3,j)*conjg(matrix(5,j)))+(conjg(matrix(3,j))*matrix(5,j))) &
                *1d0/(1d0+exp((eigenvalue(j)-mu(loop))/T))*0.5
           lx_ex=lx_ex+sqrt(12.0)*((matrix(5,j)*conjg(matrix(7,j)))+(conjg(matrix(5,j))*matrix(7,j))) &
                *1d0/(1d0+exp((eigenvalue(j)-mu(loop))/T))*0.5
           lx_ex=lx_ex+sqrt(12.0)*((matrix(7,j)*conjg(matrix(9,j)))+(conjg(matrix(7,j))*matrix(9,j))) &
                *1d0/(1d0+exp((eigenvalue(j)-mu(loop))/T))*0.5
           lx_ex=lx_ex+sqrt(10.0)*((matrix(9,j)*conjg(matrix(11,j)))+(conjg(matrix(9,j))*matrix(11,j))) &
                *1d0/(1d0+exp((eigenvalue(j)-mu(loop))/T))*0.5
           lx_ex=lx_ex+sqrt(6.0)*((matrix(11,j)*conjg(matrix(13,j)))+(conjg(matrix(11,j))*matrix(13,j))) &
                *1d0/(1d0+exp((eigenvalue(j)-mu(loop))/T))*0.5

           lx_ex=lx_ex+sqrt(6.0)*((matrix(2,j)*conjg(matrix(4,j)))+(conjg(matrix(2,j))*matrix(4,j))) &
                *1d0/(1d0+exp((eigenvalue(j)-mu(loop))/T))*0.5
           lx_ex=lx_ex+sqrt(10.0)*((matrix(4,j)*conjg(matrix(6,j)))+(conjg(matrix(4,j))*matrix(6,j))) &
                *1d0/(1d0+exp((eigenvalue(j)-mu(loop))/T))*0.5
           lx_ex=lx_ex+sqrt(12.0)*((matrix(6,j)*conjg(matrix(8,j)))+(conjg(matrix(6,j))*matrix(8,j))) &
                *1d0/(1d0+exp((eigenvalue(j)-mu(loop))/T))*0.5
           lx_ex=lx_ex+sqrt(12.0)*((matrix(8,j)*conjg(matrix(10,j)))+(conjg(matrix(8,j))*matrix(10,j))) &
                *1d0/(1d0+exp((eigenvalue(j)-mu(loop))/T))*0.5
           lx_ex=lx_ex+sqrt(10.0)*((matrix(10,j)*conjg(matrix(12,j)))+(conjg(matrix(10,j))*matrix(12,j))) &
                *1d0/(1d0+exp((eigenvalue(j)-mu(loop))/T))*0.5
           lx_ex=lx_ex+sqrt(6.0)*((matrix(12,j)*conjg(matrix(14,j)))+(conjg(matrix(12,j))*matrix(14,j))) &
                *1d0/(1d0+exp((eigenvalue(j)-mu(loop))/T))*0.5


           ly_ex=ly_ex+sqrt(6.0)*((matrix(1,j)*conjg(matrix(3,j)))-(conjg(matrix(1,j))*matrix(3,j))) &
                *1d0/(1d0+exp((eigenvalue(j)-mu(loop))/T))*0.5
           ly_ex=ly_ex+sqrt(10.0)*((matrix(3,j)*conjg(matrix(5,j)))-(conjg(matrix(3,j))*matrix(5,j))) &
                *1d0/(1d0+exp((eigenvalue(j)-mu(loop))/T))*0.5
           ly_ex=ly_ex+sqrt(12.0)*((matrix(5,j)*conjg(matrix(7,j)))-(conjg(matrix(5,j))*matrix(7,j))) &
                *1d0/(1d0+exp((eigenvalue(j)-mu(loop))/T))*0.5
           ly_ex=ly_ex+sqrt(12.0)*((matrix(7,j)*conjg(matrix(9,j)))-(conjg(matrix(7,j))*matrix(9,j))) &
                *1d0/(1d0+exp((eigenvalue(j)-mu(loop))/T))*0.5
           ly_ex=ly_ex+sqrt(10.0)*((matrix(9,j)*conjg(matrix(11,j)))-(conjg(matrix(9,j))*matrix(11,j))) &
                *1d0/(1d0+exp((eigenvalue(j)-mu(loop))/T))*0.5
           ly_ex=ly_ex+sqrt(6.0)*((matrix(11,j)*conjg(matrix(13,j)))-(conjg(matrix(11,j))*matrix(13,j))) &
                *1d0/(1d0+exp((eigenvalue(j)-mu(loop))/T))*0.5

           ly_ex=ly_ex+sqrt(6.0)*((matrix(2,j)*conjg(matrix(4,j)))-(conjg(matrix(2,j))*matrix(4,j))) &
                *1d0/(1d0+exp((eigenvalue(j)-mu(loop))/T))*0.5
           ly_ex=ly_ex+sqrt(10.0)*((matrix(4,j)*conjg(matrix(6,j)))-(conjg(matrix(4,j))*matrix(6,j))) &
                *1d0/(1d0+exp((eigenvalue(j)-mu(loop))/T))*0.5
           ly_ex=ly_ex+sqrt(12.0)*((matrix(6,j)*conjg(matrix(8,j)))-(conjg(matrix(6,j))*matrix(8,j))) &
                *1d0/(1d0+exp((eigenvalue(j)-mu(loop))/T))*0.5
           ly_ex=ly_ex+sqrt(12.0)*((matrix(8,j)*conjg(matrix(10,j)))-(conjg(matrix(8,j))*matrix(10,j))) &
                *1d0/(1d0+exp((eigenvalue(j)-mu(loop))/T))*0.5
           ly_ex=ly_ex+sqrt(10.0)*((matrix(10,j)*conjg(matrix(12,j)))-(conjg(matrix(10,j))*matrix(12,j))) &
                *1d0/(1d0+exp((eigenvalue(j)-mu(loop))/T))*0.5
           ly_ex=ly_ex+sqrt(6.0)*((matrix(12,j)*conjg(matrix(14,j)))-(conjg(matrix(12,j))*matrix(14,j))) &
                *1d0/(1d0+exp((eigenvalue(j)-mu(loop))/T))*0.5

           

           sx_ex=sx_ex+((matrix(1,j)*conjg(matrix(2,j)))+(conjg(matrix(1,j))*matrix(2,j))) &
                *1d0/(1d0+exp((eigenvalue(j)-mu(loop))/T))*0.5
           sx_ex=sx_ex+((matrix(3,j)*conjg(matrix(4,j)))+(conjg(matrix(3,j))*matrix(4,j))) &
                *1d0/(1d0+exp((eigenvalue(j)-mu(loop))/T))*0.5
           sx_ex=sx_ex+((matrix(5,j)*conjg(matrix(6,j)))+(conjg(matrix(5,j))*matrix(6,j))) &
                *1d0/(1d0+exp((eigenvalue(j)-mu(loop))/T))*0.5
           sx_ex=sx_ex+((matrix(7,j)*conjg(matrix(8,j)))+(conjg(matrix(7,j))*matrix(8,j))) &
                *1d0/(1d0+exp((eigenvalue(j)-mu(loop))/T))*0.5
           sx_ex=sx_ex+((matrix(9,j)*conjg(matrix(10,j)))+(conjg(matrix(9,j))*matrix(10,j))) &
                *1d0/(1d0+exp((eigenvalue(j)-mu(loop))/T))*0.5
           sx_ex=sx_ex+((matrix(11,j)*conjg(matrix(12,j)))+(conjg(matrix(11,j))*matrix(12,j))) &
                *1d0/(1d0+exp((eigenvalue(j)-mu(loop))/T))*0.5
           sx_ex=sx_ex+((matrix(13,j)*conjg(matrix(14,j)))+(conjg(matrix(13,j))*matrix(14,j))) &
                *1d0/(1d0+exp((eigenvalue(j)-mu(loop))/T))*0.5

           
           sy_ex=sy_ex+((matrix(1,j)*conjg(matrix(2,j)))-(conjg(matrix(1,j))*matrix(2,j))) &
                *1d0/(1d0+exp((eigenvalue(j)-mu(loop))/T))*0.5
           sy_ex=sy_ex+((matrix(3,j)*conjg(matrix(4,j)))-(conjg(matrix(3,j))*matrix(4,j))) &
                *1d0/(1d0+exp((eigenvalue(j)-mu(loop))/T))*0.5
           sy_ex=sy_ex+((matrix(5,j)*conjg(matrix(6,j)))-(conjg(matrix(5,j))*matrix(6,j))) &
                *1d0/(1d0+exp((eigenvalue(j)-mu(loop))/T))*0.5
           sy_ex=sy_ex+((matrix(7,j)*conjg(matrix(8,j)))-(conjg(matrix(7,j))*matrix(8,j))) &
                *1d0/(1d0+exp((eigenvalue(j)-mu(loop))/T))*0.5
           sy_ex=sy_ex+((matrix(9,j)*conjg(matrix(10,j)))-(conjg(matrix(9,j))*matrix(10,j))) &
                *1d0/(1d0+exp((eigenvalue(j)-mu(loop))/T))*0.5
           sy_ex=sy_ex+((matrix(11,j)*conjg(matrix(12,j)))-(conjg(matrix(11,j))*matrix(12,j))) &
                *1d0/(1d0+exp((eigenvalue(j)-mu(loop))/T))*0.5
           sy_ex=sy_ex+((matrix(13,j)*conjg(matrix(14,j)))-(conjg(matrix(13,j))*matrix(14,j))) &
                *1d0/(1d0+exp((eigenvalue(j)-mu(loop))/T))*0.5
           
           ly_ex=imag*ly_ex
           sy_ex=imag*sy_ex
           
        end do
        
              do i=1,jigen
                rho(nint(eigenvalue(i)/mesh))=rho(nint(eigenvalue(i)/mesh))+1
                 do j=1,jigen
                    Prho(j,nint(eigenvalue(i)/mesh))=Prho(j,nint(eigenvalue(i)/mesh))+abs(matrix(j,i))**2
                 end do
              end do



              do i=1,jigen
                 integral=integral+1d0/(1d0+exp((Eigenvalue(i)-mu(loop))/T))
              end do


              do m=1,jigen
                 do i=1,14 !f2
                    do j=1,14 !f1
                       if(i==j) then
                          f2f1(i,i)=f2f1(i,i)+abs(matrix(i,m))**2d0/Nmax/Nmax/(1d0+exp((eigenvalue(m)-mu(loop))/T))
                       else
                          f2f1(i,j)=f2f1(i,j)+matrix(j,m)*conjg(matrix(i,m))/Nmax/Nmax/(1d0+exp((eigenvalue(m)-mu(loop))/T))
                       end if
                    end do
                 end do
              end do



           end do
        end do



     if((theta_int==92).and.nx<3.and.ny<3.and.loop2==1) then
        print *,"f2f1_Pr="
        do i=1,14
           write(6,660) (f2f1_pr(i,j),j=1,14)
660        format(26e12.4)
        end do
     end if

     delta_matrix=0
     do i=1,14
        do j=1,14
           delta_matrix=delta_matrix+abs(f2f1_pr(i,j)-f2f1(i,j))
        end do
     end do

!!!!!!!!!!!!!!!!!!!! 恣意的な初期条件の代入 !!!!!!!!!!!!!!!!!!!!!!


     !   f2f1(:,:)=0

     !  do i=7,8,1
     !     f2f1(i,i)=1
     !     end do


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


     if((theta_int==92).and.nx<3.and.ny<3) then
        print *,"f2f1="
        do i=1,14
           write(6,650) (f2f1(i,j),j=1,14)
650        format(28e12.4)
        end do
     end if






     call heapsort(Nmax*Nmax*jigen,all_eigenvalue)



     delta=0

!     print*,loop,"loop end,N4f=",N4f(loop),"del_matrix=",delta_matrix,"integral=",integral,"/",Nmax*Nmax*jigen,"mu=",mu(loop) &
 !         ,"mixing=",mixing
  !   write(*,'(i2,a15,e20.3,a15,e10.8,a15,e10.5)') loop,"loop end,N4f=",N4f(loop),"del_matrix=",delta_matrix,"integral=",integral 
  end do


!!!!!!!!!!!!!!!!!!!!なぜかこれを入れるとうまくいく!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  
  if(loop2==1) then
     f2f1(:,:)=f2f1_loop1(:,:)
     end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     
  if(loop<loop_max) then
     count=0
     mixing=0.5d0
     eps_matrix=0.0001d0
     if(integral>Nmax*Nmax*4d0) then
        mu_plu=mu(0)
     else
        mu_min=mu(0)
     end if
  else
     eps_matrix=0.001d0
     count=count+1
     mixing=count*0.1d0
  end if
  if(count==10) exit

!  print*,loop2,"loop2 end,","integral=",integral,"mu=",mu(loop),"mu_min=",mu_min,"mu_plu=",mu_plu,"loop=",loop &
!       ,"count=",count,"N4f=",N4f(loop)
 write(*,'(i5,a20,e15.7,a10,e15.7,a10,i5,a10,e14.5)') loop2,"loop2 end,integral=",integral,"mu=",mu(loop), &
       "loop=",loop,"N4f=",N4f(loop)

 
end do



J_partition=0

do i=1,(Nmax**2)*jigen
   if(all_eigenvalue(i)-mu(loop)<-T*100d0) then
      J_partition=J_partition-T*(-(all_eigenvalue(i)-mu(loop))/T)
   else
         J_partition=J_partition-T*log(1d0+exp(-(all_eigenvalue(i)-mu(loop))/T))
      end if
end do


print*,"N4f=",N4f(loop),"Nm=",Nm
print*,"Lz=",real(Lz_ex)/Nmax/Nmax,"Ly=",real(Ly_ex)/Nmax/Nmax,"Lx=",real(Lx_ex)/Nmax/Nmax



write(9,'(4e18.7)') (theta_int-90)*1d0,(J_partition/Nmax/Nmax+mu(loop)*integral/Nmax/Nmax),loop2*1d0,N4f(loop)
end do




open(10,file='DOS1.dat',status='replace')
open(11,file='DOS2.dat',status='replace')
do i=-10*emax,10*emax


  do j=1,14,2
     rhofup(i)=rhofup(i)+prho(j,i)
  end do

  do j=2,14,2
     rhofdn(i)=rhofdn(i)+prho(j,i)
  end do

  do j=15,20,2
     rhodup(i)=rhodup(i)+prho(j,i)
  end do


  do j=16,20,2
     rhoddn(i)=rhoddn(i)+prho(j,i)
  end do




  write(10,'(23e18.7)')E(i),rho(i)/Nmax/Nmax/mesh,(prho(j,i)/Nmax/Nmax/mesh,j=1,jigen)
  write(11,'(5e18.7)')E(i),rhofup(i)/Nmax/Nmax/mesh,rhofdn(i)/Nmax/Nmax/mesh,rhodup(i)/Nmax/Nmax/mesh,rhoddn(i)/Nmax/Nmax/mesh




  if(mod(i,50)==0) then
     !     print*,"E=",E(i),"integral_d=",integral_d,"integral_f=",integral_f
  end if
end do

close(9)
close(10)
close(11)











end program

subroutine heapsort(n,array)
implicit none
integer,intent(in) :: n
double precision,intent(inout) :: array(1:n)

integer ::i,k,j,l
double precision :: t

if(n.le.0)then
  write(6,*)"Error, at heapsort"; stop
endif
if(n.eq.1)return

l=n/2+1
k=n
do while(k.ne.1)
  if(l.gt.1)then
     l=l-1
     t=array(L)
  else
     t=array(k)
     array(k)=array(1)
     k=k-1
     if(k.eq.1) then
        array(1)=t
        exit
     endif
  endif
  i=l
  j=l+l
  do while(j.le.k)
     if(j.lt.k)then
        if(array(j).lt.array(j+1))j=j+1
     endif
     if (t.lt.array(j))then
        array(i)=array(j)
        i=j
        j=j+j
     else
        j=k+1
     endif
  enddo
  array(i)=t
enddo

return
end subroutine heapsort
