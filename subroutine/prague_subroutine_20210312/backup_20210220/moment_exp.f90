subroutine moment_exp(temp,nr,nfe,imax,zunit0,zspin,zorbit,thetamag,phimag,theta0,phi0,Femag,dmult,sys,method,nJexch,volume,part0,eigen0,wwoK1Fe,rdch)
implicit none
integer(4),parameter :: nrmax=30
real(8),parameter :: mu0=4.0d0*dacos(-1.0d0)*1.0d-7 ! (NA^-2)
real(8),parameter :: muB=9.2740100783d-24 ! (Merg/kOe) or (kerg/Oe) or (J/T)
real(8),parameter :: pi=dacos(-1.0d0)
integer(4) :: ir,nr,nfe,ll,m,imax,jj,ii
real(8) :: stspin0(3,nrmax),mspin(3,nrmax),morbit(3,nrmax),eigen0(imax),part0(nrmax),temp,volume,dmult
real(8) :: def,ghi,Rmag0,Femag,Femag0,magtot(3),thetamag,phimag,thetatot,phitot,abstot,theta0,phi0
real(8) :: thetaj0(3),phij0(3),absj0(3),thetaS0(3),phiS0(3),absS0(3),thetaL0(3),phiL0(3),absL0(3)
complex(kind(0d0)) :: zunit0(imax,imax,nr),zspin(imax,imax,3),zorbit(imax,imax,3)
character(50) :: sys,method,nJexch,wwoK1Fe,rdch

! angular momentum
stspin0=0.0d0
mspin=0.0d0
morbit=0.0d0
do ir=1,nr
   do ll=1,3
      do m=1,imax
         def=0.0d0
         ghi=0.0d0
         do jj=1,imax
            do ii=1,imax
               def=def-2.0d0*dreal(conjg(zunit0(ii,m,ir))*zspin(ii,jj,ll)*zunit0(jj,m,ir))
               ghi=ghi-dreal(conjg(zunit0(ii,m,ir))*zorbit(ii,jj,ll)*zunit0(jj,m,ir))
            end do
         end do
         if(temp.eq.0.0d0) then
            mspin(ll,ir)=def
            morbit(ll,ir)=ghi
            exit
         else
            mspin(ll,ir)=mspin(ll,ir)+def*dexp(-(eigen0(m)-eigen0(1))/temp)/part0(ir)                        
            morbit(ll,ir)=morbit(ll,ir)+ghi*dexp(-(eigen0(m)-eigen0(1))/temp)/part0(ir)                        
         end if
      end do ! m energy level
      stspin0(ll,ir)=mspin(ll,ir)+morbit(ll,ir)
   end do ! for spin component
end do ! for ir=1,2

! magnetization of RE along the applied field direction
Rmag0=0.0d0
do ir=1,nr
   Rmag0=Rmag0+(dsin(thetamag)*(stspin0(1,ir)*dcos(phimag)+stspin0(2,ir)*dsin(phimag))+stspin0(3,ir)*dcos(thetamag))
end do
Rmag0=Rmag0*dmult

! magnetization of Fe along the direction of applied field
Femag0=Femag*(dsin(thetamag)*dsin(theta0)*dcos(phimag-phi0)+dcos(thetamag)*dcos(theta0))

! total magnetization
magtot=0.0d0
magtot(1)=Femag*dsin(theta0)*dcos(phi0)
magtot(2)=Femag*dsin(theta0)*dsin(phi0)
magtot(3)=Femag*dcos(theta0)
do ir=1,nr
   do ll=1,3
      magtot(ll)=magtot(ll)+dmult*stspin0(ll,ir)
   end do
end do
open(169,file=""//trim(adjustl(sys))//"_"//trim(adjustl(method))//"_"//trim(adjustl(wwoK1Fe))//"_&
     MsT_nJex"//trim(adjustl(nJexch))//"_rd"//trim(adjustl(rdch))//".txt") ! "# temperature (K), magnetization (per unit cell)
write(169,'(F10.2,3F15.5)') temp,magtot(1),magtot(2),magtot(3)

! angle of total magnetization
call vector_to_angle(magtot,thetatot,phitot,abstot)
do ir=1,nr
   call vector_to_angle(stspin0(1,ir),thetaj0(ir),phij0(ir),absj0(ir))
   call vector_to_angle(mspin(1,ir),thetaS0(ir),phiS0(ir),absS0(ir))
   call vector_to_angle(morbit(1,ir),thetaL0(ir),phiL0(ir),absL0(ir))
end do

! magnetic moments
open(171,file=""//trim(adjustl(sys))//"_"//trim(adjustl(method))//"_"//trim(adjustl(wwoK1Fe))//"_&
     Fe_theta-phi_nJex"//trim(adjustl(nJexch))//"_rd"//trim(adjustl(rdch))//".txt") !"# temperature (K), theta (degree), phi (degree), abs"
open(172,file=""//trim(adjustl(sys))//"_"//trim(adjustl(method))//"_"//trim(adjustl(wwoK1Fe))//"_&
     RE1-2_J_theta-phi_nJex"//trim(adjustl(nJexch))//"_rd"//trim(adjustl(rdch))//".txt") ! "# temperature (K), theta (degree), phi (degree), abs for RE 1-2"
open(173,file=""//trim(adjustl(sys))//"_"//trim(adjustl(method))//"_"//trim(adjustl(wwoK1Fe))//"_&
     RE1-2_S_theta-phi_nJex"//trim(adjustl(nJexch))//"_rd"//trim(adjustl(rdch))//".txt") ! "# temperature (K), theta (degree), phi (degree), abs for RE 1-2"
open(174,file=""//trim(adjustl(sys))//"_"//trim(adjustl(method))//"_"//trim(adjustl(wwoK1Fe))//"_&
     RE1-2_L_theta-phi_nJex"//trim(adjustl(nJexch))//"_rd"//trim(adjustl(rdch))//".txt") ! "# temperature (K), theta (degree), phi (degree), abs for RE 1-2"
open(175,file=""//trim(adjustl(sys))//"_"//trim(adjustl(method))//"_"//trim(adjustl(wwoK1Fe))//"_&
     Mtot_theta-phi_nJex"//trim(adjustl(nJexch))//"_rd"//trim(adjustl(rdch))//".txt")

write(171,"(4f12.5)") temp,theta0/pi*180.0d0,phi0/pi*180.0d0,Femag/(dmult*dble(nfe))
write(172,"(7f12.5)") temp,thetaj0(1)/pi*180.0d0,phij0(1)/pi*180.0d0,absj0(1),thetaj0(2)/pi*180.0d0,phij0(2)/pi*180.0d0,absj0(2)
write(173,"(7f12.5)") temp,thetaS0(1)/pi*180.0d0,phiS0(1)/pi*180.0d0,absS0(1),thetaS0(2)/pi*180.0d0,phiS0(2)/pi*180.0d0,absS0(2)
write(174,"(7f12.5)") temp,thetaL0(1)/pi*180.0d0,phiL0(1)/pi*180.0d0,absL0(1),thetaL0(2)/pi*180.0d0,phiL0(2)/pi*180.0d0,absL0(2)
if(nr.gt.2) then
   open(182,file=""//trim(adjustl(sys))//"_"//trim(adjustl(method))//"_"//trim(adjustl(wwoK1Fe))//"_&
        RE3-4_J_theta-phi_nJex"//trim(adjustl(nJexch))//"_rd"//trim(adjustl(rdch))//".txt") !"# temperature (K), theta (degree), phi (degree), abs for RE 3-4"
   open(183,file=""//trim(adjustl(sys))//"_"//trim(adjustl(method))//"_"//trim(adjustl(wwoK1Fe))//"_&
        RE3-4_S_theta-phi_nJex"//trim(adjustl(nJexch))//"_rd"//trim(adjustl(rdch))//".txt") !"# temperature (K), theta (degree), phi (degree), abs for RE 3-4"
   open(184,file=""//trim(adjustl(sys))//"_"//trim(adjustl(method))//"_"//trim(adjustl(wwoK1Fe))//"_&
        RE3-4_L_theta-phi_nJex"//trim(adjustl(nJexch))//"_rd"//trim(adjustl(rdch))//".txt") !"# temperature (K), theta (degree), phi (degree), abs for RE 3-4"
   write(182,"(7f12.5)") temp,thetaj0(3)/pi*180.0d0,phij0(3)/pi*180.0d0,absj0(3),thetaj0(4)/pi*180.0d0,phij0(4)/pi*180.0d0,absj0(4)
   write(183,"(7f12.5)") temp,thetaS0(3)/pi*180.0d0,phiS0(3)/pi*180.0d0,absS0(3),thetaS0(4)/pi*180.0d0,phiS0(4)/pi*180.0d0,absS0(4)
   write(184,"(7f12.5)") temp,thetaL0(3)/pi*180.0d0,phiL0(3)/pi*180.0d0,absL0(3),thetaL0(4)/pi*180.0d0,phiL0(4)/pi*180.0d0,absL0(4)
end if
write(175,"(4f12.5)") temp,thetatot/pi*180.0d0,phitot/pi*180.0d0,mu0*abstot*muB/volume

! history
do ir=1,nr
   write(6,"(a20,3f12.5)") "RE theta-phi-abs",thetaj0(ir)/pi*180.0d0,phij0(ir)/pi*180.0d0,absj0(ir)
end do
write(6,"(a20,3f12.5)") "Fe theta-phi-abs",theta0/pi*180.0d0,phi0/pi*180.0d0,Femag/(dmult*dble(nfe))
write(6,"(a20,3f12.5)") "total theta-phi-abs",thetatot/pi*180.0d0,phitot/pi*180.0d0,abstot
write(6,*) "magnetization [muB/unitcell, T]"
write(6,*) Femag0+Rmag0,(Femag0+Rmag0)*muB*mu0/volume,mu0*abstot*muB/volume

end subroutine moment_exp
