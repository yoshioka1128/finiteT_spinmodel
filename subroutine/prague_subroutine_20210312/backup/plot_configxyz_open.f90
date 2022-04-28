subroutine plot_configxyz_SmFe12(itime,nunit,nr,nfe,amp,vr,vfe,nunitx,nunity,nunitz,nrtot,nfetot,irsub,irunit,ifesub,ifeunit,&
     ax,cx,nplotx,nploty,nplotz,magR,magFe,magRT,magRST,magFeT,rrpt,rfept,pbcz,irpt,ifept)

implicit none
integer(4) :: nunit,nr,nfe,itime,iunit,i,nunitx,nunity,nunitz,ix,iy,iz,ir,ife,nrtot,nfetot
integer(4) :: ife0,ir0,ix0,iy0,iz0,irpt(nr,nunitx,nunity,nunitz),ifept(nfe,nunitx,nunity,nunitz)
real(8) :: rrpt(nr,0:nunitx+1,0:nunity+1,0:nunitz+1,3),rfept(nfe,0:nunitx+1,0:nunity+1,0:nunitz+1,3)
integer(4) :: nplotx(2),nploty(2),nplotz(2),idim,irsub(nr*nunit),ifesub(nfe*nunit),irunit(nr*nunit,3),ifeunit(nfe*nunit,3)
real(8) :: amp,ampS,ampL,vr(nr*nunit,3),vfe(nfe*nunit,3),ax,cx,posr(nr,3),posfe(nfe,3)
real(8) :: magR(nr*nunit),magFe(nr*nunit),magRT,magFeT,magRST,magRLT,signvr
character(5) :: time
character(10) :: pbcz
magRLT=magRT-magRST
!signvr=magRST/magRT*abs(magRT/magRST)
ampS=abs(magRST/magFeT)*magRT*magRST/abs(magRT*magRST)*amp
ampL=abs(magRLT/magFeT)*magRT*magRlT/abs(magRT*magRLT)*amp

write(time,'(I5)') itime
open(8,file="movie_source/TDRL"//trim(adjustl(time))//".txt") 
open(9,file="movie_source/TDRS"//trim(adjustl(time))//".txt") 
open(10,file="movie_source/TDFe"//trim(adjustl(time))//".txt") 

do iz=nplotz(1),nplotz(2)+1
   do iy=nploty(1),nploty(2)+1
      do ix=nplotx(1),nplotx(2)+1
         do i=1,nr
            if(rrpt(i,ix,iy,iz,1).le.dble(nplotx(2)*ax).and.&
                 rrpt(i,ix,iy,iz,2).le.dble(nploty(2)*ax).and.&
                 rrpt(i,ix,iy,iz,3).le.dble(nplotz(2)*cx)) then
               ir=irpt(i,ix,iy,iz)
               write(8,"(6f10.5)") rrpt(i,ix,iy,iz,1),rrpt(i,ix,iy,iz,2),rrpt(i,ix,iy,iz,3),&
                    ampL*vr(ir,1),ampL*vr(ir,2),ampL*vr(ir,3) ! fsite, Nd1
               write(9,"(6f10.5)") rrpt(i,ix,iy,iz,1),rrpt(i,ix,iy,iz,2),rrpt(i,ix,iy,iz,3),&
                    ampS*vr(ir,1),ampS*vr(ir,2),ampS*vr(ir,3) ! fsite, Nd1
            end if
         end do
         do i=1,nfe
            if(rfept(i,ix,iy,iz,1).le.dble(nplotx(2)*ax).and.&
                 rfept(i,ix,iy,iz,2).le.dble(nploty(2)*ax).and.&
                 rfept(i,ix,iy,iz,3).le.dble(nplotz(2)*cx)) then
               ife=ifept(i,ix,iy,iz)
               write(10,"(6f10.5)") rfept(i,ix,iy,iz,1),rfept(i,ix,iy,iz,2),rfept(i,ix,iy,iz,3),&
                    amp*magFe(ife)/magFeT*vfe(ife,1),amp*magFe(ife)/magFeT*vfe(ife,2),amp*magFe(ife)/magFeT*vfe(ife,3) ! fsite, Nd1
            end if
         end do

      end do
   end do
end do

end subroutine plot_configxyz_SmFe12




subroutine plot_configxyz_SmFe12_old(itime,nunit,nr,nfe,amp,vr,vfe,nunitx,nunity,nunitz,nrtot,nfetot,irsub,irunit,ifesub,ifeunit,&
     ax,cx,nplotx,nploty,nplotz,magR,magFe,magRT,magRST,magFeT,rrpt,rfept,pbcz)

implicit none
integer(4) :: nunit,nr,nfe,itime,iunit,i,nunitx,nunity,nunitz,ix,iy,iz,ir,ife,nrtot,nfetot
integer(4) :: ife0,ir0,ix0,iy0,iz0
real(8) :: rrpt(nr,0:nunitx+1,0:nunity+1,0:nunitz+1,3),rfept(nfe,0:nunitx+1,0:nunity+1,0:nunitz+1,3)
integer(4) :: nplotx(2),nploty(2),nplotz(2),idim,irsub(nr*nunit),ifesub(nfe*nunit),irunit(nr*nunit,3),ifeunit(nfe*nunit,3)
real(8) :: amp,ampS,ampL,vr(nr*nunit,3),vfe(nfe*nunit,3),ax,cx,posr(nr,3),posfe(nfe,3)
real(8) :: magR(nr*nunit),magFe(nr*nunit),magRT,magFeT,magRST,magRLT,signvr
character(5) :: time
character(10) :: pbcz
magRLT=magRT-magRST
signvr=magRST/magRT*abs(magRT/magRST)
ampS=signvr*abs(magRST/magFeT)*amp
ampL=-signvr*abs(magRLT/magFeT)*amp

write(time,'(I5)') itime
open(8,file="movie_source/TDRL"//trim(adjustl(time))//".txt") 
open(9,file="movie_source/TDRS"//trim(adjustl(time))//".txt") 
open(10,file="movie_source/TDFe"//trim(adjustl(time))//".txt") 
ir=1
do ir=1,nrtot
   i=irsub(ir)
   ix0=irunit(ir,1)
   iy0=irunit(ir,2)
   iz0=irunit(ir,3)
   if(ix0.ge.nplotx(1).and.ix0.le.nplotx(2).and.iy0.ge.nploty(1).and.iy0.le.nploty(2).and.iz0.ge.nplotz(1).and.iz0.le.nplotz(2)) then
      ix=ix0-nplotx(1)+1
      iy=iy0-nploty(1)+1
      iz=iz0-nplotz(1)+1
      if(pbcz.eq."off".and.nplotz(1).eq.1) then
         write(8,"(6f10.5)") rrpt(i,ix,iy,iz,1),rrpt(i,ix,iy,iz,2),rrpt(i,ix,iy,iz,3)-cx,&
              ampL*vr(ir,1),ampL*vr(ir,2),ampL*vr(ir,3) ! fsite, Nd1
      else
         write(6,"(6f10.5)") rrpt(i,ix,iy,iz,1),rrpt(i,ix,iy,iz,2),rrpt(i,ix,iy,iz,3),&
              ampL*vr(ir,1),ampL*vr(ir,2),ampL*vr(ir,3) ! fsite, Nd1
         write(8,"(6f10.5)") rrpt(i,ix,iy,iz,1),rrpt(i,ix,iy,iz,2),rrpt(i,ix,iy,iz,3),&
              ampL*vr(ir,1),ampL*vr(ir,2),ampL*vr(ir,3) ! fsite, Nd1
      end if
      if(pbcz.eq."off".and.nplotz(1).eq.1) then
         write(9,"(6f10.5)") rrpt(i,ix,iy,iz,1),rrpt(i,ix,iy,iz,2),rrpt(i,ix,iy,iz,3)-cx,&
              ampS*vr(ir,1),ampS*vr(ir,2),ampS*vr(ir,3) ! fsite, Nd1
      else
         write(9,"(6f10.5)") rrpt(i,ix,iy,iz,1),rrpt(i,ix,iy,iz,2),rrpt(i,ix,iy,iz,3),&
              ampS*vr(ir,1),ampS*vr(ir,2),ampS*vr(ir,3) ! fsite, Nd1
      end if
   end if
end do
stop
do ife=1,nfetot
   i=ifesub(ife)
   ix0=ifeunit(ife,1)
   iy0=ifeunit(ife,2)
   iz0=ifeunit(ife,3)
   if(ix0.ge.nplotx(1).and.ix0.le.nplotx(2).and.iy0.ge.nploty(1).and.iy0.le.nploty(2).and.iz0.ge.nplotz(1).and.iz0.le.nplotz(2)) then
      ix=ix0-nplotx(1)+1
      iy=iy0-nploty(1)+1
      iz=iz0-nplotz(1)+1
      if(pbcz.eq."off".and.nplotz(1).eq.1) then
         write(10,"(6f10.5)") rfept(i,ix,iy,iz,1),rfept(i,ix,iy,iz,2),rfept(i,ix,iy,iz,3)-cx,&
              amp*magFe(ife)/magFeT*vfe(ife,1),amp*magFe(ife)/magFeT*vfe(ife,2),amp*magFe(ife)/magFeT*vfe(ife,3) ! fsite, Nd1
      else
         write(10,"(6f10.5)") rfept(i,ix,iy,iz,1),rfept(i,ix,iy,iz,2),rfept(i,ix,iy,iz,3),&
              amp*magFe(ife)/magFeT*vfe(ife,1),amp*magFe(ife)/magFeT*vfe(ife,2),amp*magFe(ife)/magFeT*vfe(ife,3) ! fsite, Nd1
      end if
   end if
end do

end subroutine plot_configxyz_SmFe12_old


subroutine plot_configxyz_open(itime,nunit,nr,nfe,amp,vr,vfe,nunitx,nunity,nunitz,nrtot,nfetot,irsub,irunit,ifesub,ifeunit,&
     ax,cx,nplotx,nploty,nplotz,magR,magFe,magRT,magFeT,rrpt,rfept,pbcz)

implicit none
integer(4) :: nunit,nr,nfe,itime,iunit,i,nunitx,nunity,nunitz,ix,iy,iz,ir,ife,nrtot,nfetot
integer(4) :: ife0,ir0,ix0,iy0,iz0
real(8) :: rrpt(nr,0:nunitx+1,0:nunity+1,0:nunitz+1,3),rfept(nfe,0:nunitx+1,0:nunity+1,0:nunitz+1,3)
integer(4) :: nplotx(2),nploty(2),nplotz(2),idim,irsub(nr*nunit),ifesub(nfe*nunit),irunit(nr*nunit,3),ifeunit(nfe*nunit,3)
real(8) :: amp,vr(nr*nunit,3),vfe(nfe*nunit,3),ax,cx,posr(nr,3),posfe(nfe,3)
real(8) :: magR(nr*nunit),magFe(nr*nunit),magRT,magFeT
character(5) :: time
character(10) :: pbcz

write(time,'(I5)') itime
open(8,file="movie_source/TDf"//trim(adjustl(time))//".txt") 
open(9,file="movie_source/TDg"//trim(adjustl(time))//".txt") 
open(10,file="movie_source/TDFe"//trim(adjustl(time))//".txt") 

do ir=1,nrtot
   i=irsub(ir)
   ix0=irunit(ir,1)
   iy0=irunit(ir,2)
   iz0=irunit(ir,3)
   if(ix0.ge.nplotx(1).and.ix0.le.nplotx(2).and.iy0.ge.nploty(1).and.iy0.le.nploty(2).and.iz0.ge.nplotz(1).and.iz0.le.nplotz(2)) then
      ix=ix0-nplotx(1)+1
      iy=iy0-nploty(1)+1
      iz=iz0-nplotz(1)+1
      if(i.eq.1.or.i.eq.2.or.i.eq.5.or.i.eq.6) then
         if(pbcz.eq."off".and.nplotz(1).eq.1) then
            write(8,"(6f10.5)") rrpt(i,ix,iy,iz,1),rrpt(i,ix,iy,iz,2),rrpt(i,ix,iy,iz,3)-cx,&
                 amp*magR(ir)/magRT*vr(ir,1),amp*magR(ir)/magRT*vr(ir,2),amp*magR(ir)/magRT*vr(ir,3) ! fsite, Nd1
         else
            write(8,"(6f10.5)") rrpt(i,ix,iy,iz,1),rrpt(i,ix,iy,iz,2),rrpt(i,ix,iy,iz,3),&
                 amp*magR(ir)/magRT*vr(ir,1),amp*magR(ir)/magRT*vr(ir,2),amp*magR(ir)/magRT*vr(ir,3) ! fsite, Nd1
         end if
      end if
      if(i.eq.3.or.i.eq.4.or.i.eq.7.or.i.eq.8) then
         if(pbcz.eq."off".and.nplotz(1).eq.1) then
            write(9,"(6f10.5)") rrpt(i,ix,iy,iz,1),rrpt(i,ix,iy,iz,2),rrpt(i,ix,iy,iz,3)-cx,&
                 amp*magR(ir)/magRT*vr(ir,1),amp*magR(ir)/magRT*vr(ir,2),amp*magR(ir)/magRT*vr(ir,3) ! fsite, Nd1
         else
            write(9,"(6f10.5)") rrpt(i,ix,iy,iz,1),rrpt(i,ix,iy,iz,2),rrpt(i,ix,iy,iz,3),&
                 amp*magR(ir)/magRT*vr(ir,1),amp*magR(ir)/magRT*vr(ir,2),amp*magR(ir)/magRT*vr(ir,3) ! fsite, Nd1
         end if
      end if
   end if
end do

do ife=1,nfetot
   i=ifesub(ife)
   ix0=ifeunit(ife,1)
   iy0=ifeunit(ife,2)
   iz0=ifeunit(ife,3)
   if(ix0.ge.nplotx(1).and.ix0.le.nplotx(2).and.iy0.ge.nploty(1).and.iy0.le.nploty(2).and.iz0.ge.nplotz(1).and.iz0.le.nplotz(2)) then
      ix=ix0-nplotx(1)+1
      iy=iy0-nploty(1)+1
      iz=iz0-nplotz(1)+1
      if(pbcz.eq."off".and.nplotz(1).eq.1) then
         write(10,"(6f10.5)") rfept(i,ix,iy,iz,1),rfept(i,ix,iy,iz,2),rfept(i,ix,iy,iz,3)-cx,&
              amp*magFe(ife)/magFeT*vfe(ife,1),amp*magFe(ife)/magFeT*vfe(ife,2),amp*magFe(ife)/magFeT*vfe(ife,3) ! fsite, Nd1
      else
         write(10,"(6f10.5)") rfept(i,ix,iy,iz,1),rfept(i,ix,iy,iz,2),rfept(i,ix,iy,iz,3),&
              amp*magFe(ife)/magFeT*vfe(ife,1),amp*magFe(ife)/magFeT*vfe(ife,2),amp*magFe(ife)/magFeT*vfe(ife,3) ! fsite, Nd1
      end if
   end if
end do

end subroutine plot_configxyz_open


subroutine plot_configxyz_temp(itime,nunit,nr,nfe,amp,vr,vfe,nunitx,nunity,nunitz,nrtot,nfetot,irsub,irunit,ifesub,ifeunit,&
     ax,cx,nplotx,nploty,nplotz,magR,magFe,magRT,magFeT,rrpt,rfept,pbcz,Klcr,Tch,Hch,nunitxch,nunitych,nunitzch,badch)

implicit none
integer(4) :: nunit,nr,nfe,itime,iunit,i,nunitx,nunity,nunitz,ix,iy,iz,ir,ife,nrtot,nfetot
integer(4) :: ife0,ir0,ix0,iy0,iz0
real(8) :: rrpt(nr,0:nunitx+1,0:nunity+1,0:nunitz+1,3),rfept(nfe,0:nunitx+1,0:nunity+1,0:nunitz+1,3)
integer(4) :: nplotx(2),nploty(2),nplotz(2),idim,irsub(nr*nunit),ifesub(nfe*nunit),irunit(nr*nunit,3),ifeunit(nfe*nunit,3)
real(8) :: amp,vr(nr*nunit,3),vfe(nfe*nunit,3),ax,cx,posr(nr,3),posfe(nfe,3)
real(8) :: magR(nr*nunit),magFe(nr*nunit),magRT,magFeT,Klcr(9,nr*nunit)
character(5) :: time
character(10) :: pbcz
character(20) :: Tch,Hch,nunitxch,nunitych,nunitzch,badch

write(time,'(I5)') itime
open(86,file="movie_source/TDftemp_T"//trim(adjustl(Tch))//"_H"//trim(adjustl(Hch))//"_&
     "//trim(adjustl(nunitxch))//"x"//trim(adjustl(nunitych))//"x"//trim(adjustl(nunitzch))//"_&
     n"//trim(adjustl(badch))//"%.txt",status="replace") 
open(96,file="movie_source/TDgtemp_T"//trim(adjustl(Tch))//"_H"//trim(adjustl(Hch))//"_&
     "//trim(adjustl(nunitxch))//"x"//trim(adjustl(nunitych))//"x"//trim(adjustl(nunitzch))//"_&
     n"//trim(adjustl(badch))//"%.txt",status="replace")
open(106,file="movie_source/TDFetemp_T"//trim(adjustl(Tch))//"_H"//trim(adjustl(Hch))//"_&
     "//trim(adjustl(nunitxch))//"x"//trim(adjustl(nunitych))//"x"//trim(adjustl(nunitzch))//"_&
     n"//trim(adjustl(badch))//"%.txt",status="replace")
open(60,file="T"//trim(adjustl(Tch))//"_H"//trim(adjustl(Hch))//"_&
     n"//trim(adjustl(badch))//"_&
     "//trim(adjustl(nunitxch))//"x"//trim(adjustl(nunitych))//"x"//trim(adjustl(nunitzch))//"_&
     config.txt",status="replace")

do ir=1,nrtot
   i=irsub(ir)
   ix0=irunit(ir,1)
   iy0=irunit(ir,2)
   iz0=irunit(ir,3)
   if(ix0.ge.nplotx(1).and.ix0.le.nplotx(2).and.iy0.ge.nploty(1).and.iy0.le.nploty(2).and.iz0.ge.nplotz(1).and.iz0.le.nplotz(2)) then
      ix=ix0-nplotx(1)+1
      iy=iy0-nploty(1)+1
      iz=iz0-nplotz(1)+1
      if(Klcr(1,ir).lt.0.0d0) then
         if(i.eq.1.or.i.eq.2.or.i.eq.5.or.i.eq.6) then
            if(pbcz.eq."off".and.nplotz(1).eq.1) then
               write(86,"(6f10.5)") rrpt(i,ix,iy,iz,1),rrpt(i,ix,iy,iz,2),rrpt(i,ix,iy,iz,3)-cx,&
                    amp*magR(ir)/magRT*vr(ir,1),amp*magR(ir)/magRT*vr(ir,2),amp*magR(ir)/magRT*vr(ir,3) ! fsite, Nd1
            else
               write(86,"(6f10.5)") rrpt(i,ix,iy,iz,1),rrpt(i,ix,iy,iz,2),rrpt(i,ix,iy,iz,3),&
                    amp*magR(ir)/magRT*vr(ir,1),amp*magR(ir)/magRT*vr(ir,2),amp*magR(ir)/magRT*vr(ir,3) ! fsite, Nd1
            end if
         end if
      else
         if(i.eq.3.or.i.eq.4.or.i.eq.7.or.i.eq.8) then
            if(pbcz.eq."off".and.nplotz(1).eq.1) then
               write(96,"(6f10.5)") rrpt(i,ix,iy,iz,1),rrpt(i,ix,iy,iz,2),rrpt(i,ix,iy,iz,3)-cx,&
                    amp*magR(ir)/magRT*vr(ir,1),amp*magR(ir)/magRT*vr(ir,2),amp*magR(ir)/magRT*vr(ir,3) ! fsite, Nd1
            else
               write(96,"(6f10.5)") rrpt(i,ix,iy,iz,1),rrpt(i,ix,iy,iz,2),rrpt(i,ix,iy,iz,3),&
                    amp*magR(ir)/magRT*vr(ir,1),amp*magR(ir)/magRT*vr(ir,2),amp*magR(ir)/magRT*vr(ir,3) ! fsite, Nd1
            end if
         end if
      end if
   end if
end do

do ife=1,nfetot
   i=ifesub(ife)
   ix0=ifeunit(ife,1)
   iy0=ifeunit(ife,2)
   iz0=ifeunit(ife,3)
   if(ix0.ge.nplotx(1).and.ix0.le.nplotx(2).and.iy0.ge.nploty(1).and.iy0.le.nploty(2).and.iz0.ge.nplotz(1).and.iz0.le.nplotz(2)) then
      ix=ix0-nplotx(1)+1
      iy=iy0-nploty(1)+1
      iz=iz0-nplotz(1)+1
      if(pbcz.eq."off".and.nplotz(1).eq.1) then
         write(106,"(6f10.5)") rfept(i,ix,iy,iz,1),rfept(i,ix,iy,iz,2),rfept(i,ix,iy,iz,3)-cx,&
              amp*magFe(ife)/magFeT*vfe(ife,1),amp*magFe(ife)/magFeT*vfe(ife,2),amp*magFe(ife)/magFeT*vfe(ife,3) ! fsite, Nd1
      else
         write(106,"(6f10.5)") rfept(i,ix,iy,iz,1),rfept(i,ix,iy,iz,2),rfept(i,ix,iy,iz,3),&
              amp*magFe(ife)/magFeT*vfe(ife,1),amp*magFe(ife)/magFeT*vfe(ife,2),amp*magFe(ife)/magFeT*vfe(ife,3) ! fsite, Nd1
      end if
   end if
end do
close(86)
close(96)
close(106)

do ir=1,nrtot
   write(60,"(I5,3f25.17)") ir,vr(ir,1),vr(ir,2),vr(ir,3)
end do
do ife=1,nfetot
   write(60,"(I5,3f25.17)") ife,vfe(ife,1),vfe(ife,2),vfe(ife,3)
end do
close(60)

end subroutine plot_configxyz_temp


subroutine plot_configxyz_001(itime,nunit,nr,nfe,amp,vr,vfe,nunitx,nunity,nunitz,nrtot,nfetot,irsub,irunit,ifesub,ifeunit,&
     ax,cx,nplotx,nploty,nplotz,magR,magFe,magRT,magFeT,rrpt,rfept,pbcz,Klcr)

implicit none
integer(4) :: nunit,nr,nfe,itime,iunit,i,nunitx,nunity,nunitz,ix,iy,iz,ir,ife,nrtot,nfetot
integer(4) :: ife0,ir0,ix0,iy0,iz0
real(8) :: rrpt(nr,0:nunitx+1,0:nunity+1,0:nunitz+1,3),rfept(nfe,0:nunitx+1,0:nunity+1,0:nunitz+1,3)
integer(4) :: nplotx(2),nploty(2),nplotz(2),idim,irsub(nr*nunit),ifesub(nfe*nunit),irunit(nr*nunit,3),ifeunit(nfe*nunit,3)
real(8) :: amp,vr(nr*nunit,3),vfe(nfe*nunit,3),ax,cx,posr(nr,3),posfe(nfe,3)
real(8) :: magR(nr*nunit),magFe(nr*nunit),magRT,magFeT,Klcr(9,nr*nunit)
character(5) :: time
character(10) :: pbcz

write(time,'(I5)') itime
open(8,file="movie_source/TDf"//trim(adjustl(time))//".txt") 
open(9,file="movie_source/TDg"//trim(adjustl(time))//".txt") 
open(10,file="movie_source/TDFe"//trim(adjustl(time))//".txt") 

do ir=1,nrtot
   i=irsub(ir)
   ix0=irunit(ir,1)
   iy0=irunit(ir,2)
   iz0=irunit(ir,3)
   if(ix0.ge.nplotx(1).and.ix0.le.nplotx(2).and.iy0.ge.nploty(1).and.iy0.le.nploty(2).and.iz0.ge.nplotz(1).and.iz0.le.nplotz(2)) then
      ix=ix0-nplotx(1)+1
      iy=iy0-nploty(1)+1
      iz=iz0-nplotz(1)+1
      if(Klcr(1,ir).lt.0.0d0) then
         if(i.eq.1.or.i.eq.2.or.i.eq.5.or.i.eq.6) then
            if(pbcz.eq."off".and.nplotz(1).eq.1) then
               write(8,"(6f10.5)") rrpt(i,ix,iy,iz,1),rrpt(i,ix,iy,iz,2),rrpt(i,ix,iy,iz,3)-cx,&
                    amp*magR(ir)/magRT*vr(ir,1),amp*magR(ir)/magRT*vr(ir,2),amp*magR(ir)/magRT*vr(ir,3) ! fsite, Nd1
            else
               write(8,"(6f10.5)") rrpt(i,ix,iy,iz,1),rrpt(i,ix,iy,iz,2),rrpt(i,ix,iy,iz,3),&
                    amp*magR(ir)/magRT*vr(ir,1),amp*magR(ir)/magRT*vr(ir,2),amp*magR(ir)/magRT*vr(ir,3) ! fsite, Nd1
            end if
         end if
      else
         if(i.eq.3.or.i.eq.4.or.i.eq.7.or.i.eq.8) then
            if(pbcz.eq."off".and.nplotz(1).eq.1) then
               write(9,"(6f10.5)") rrpt(i,ix,iy,iz,1),rrpt(i,ix,iy,iz,2),rrpt(i,ix,iy,iz,3)-cx,&
                    amp*magR(ir)/magRT*vr(ir,1),amp*magR(ir)/magRT*vr(ir,2),amp*magR(ir)/magRT*vr(ir,3) ! fsite, Nd1
            else
               write(9,"(6f10.5)") rrpt(i,ix,iy,iz,1),rrpt(i,ix,iy,iz,2),rrpt(i,ix,iy,iz,3),&
                    amp*magR(ir)/magRT*vr(ir,1),amp*magR(ir)/magRT*vr(ir,2),amp*magR(ir)/magRT*vr(ir,3) ! fsite, Nd1
            end if
         end if
      end if
   end if
end do

do ife=1,nfetot
   i=ifesub(ife)
   ix0=ifeunit(ife,1)
   iy0=ifeunit(ife,2)
   iz0=ifeunit(ife,3)
   if(ix0.ge.nplotx(1).and.ix0.le.nplotx(2).and.iy0.ge.nploty(1).and.iy0.le.nploty(2).and.iz0.ge.nplotz(1).and.iz0.le.nplotz(2)) then
      ix=ix0-nplotx(1)+1
      iy=iy0-nploty(1)+1
      iz=iz0-nplotz(1)+1
      if(pbcz.eq."off".and.nplotz(1).eq.1) then
         write(10,"(6f10.5)") rfept(i,ix,iy,iz,1),rfept(i,ix,iy,iz,2),rfept(i,ix,iy,iz,3)-cx,&
              amp*magFe(ife)/magFeT*vfe(ife,1),amp*magFe(ife)/magFeT*vfe(ife,2),amp*magFe(ife)/magFeT*vfe(ife,3) ! fsite, Nd1
      else
         write(10,"(6f10.5)") rfept(i,ix,iy,iz,1),rfept(i,ix,iy,iz,2),rfept(i,ix,iy,iz,3),&
              amp*magFe(ife)/magFeT*vfe(ife,1),amp*magFe(ife)/magFeT*vfe(ife,2),amp*magFe(ife)/magFeT*vfe(ife,3) ! fsite, Nd1
      end if
   end if
end do

end subroutine plot_configxyz_001
