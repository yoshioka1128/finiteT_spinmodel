program qtl
implicit none
real(8) :: a,b,c,d,pdfnum(3,2)
real(8) :: dpdf(3,-3:3,2)
real(8) :: pi=dacos(-1.0d0),abc,are(2),aim
integer(4) :: i,j,k,l,m,lmax,iabc,inum,iatommax,iatnum,ia,ib,ic,id,ie,if,ist(100)
character(len=32) atom,case
character(32) rmtch,ach,bch,cch,dch
!integer(4) :: getcwd,status
!character(64) dirname

write(6,*) "directory name ?"
read(5,*) case
! write file
open(12,file=''//trim(adjustl(case))//'_population.txt')

! read file
open(10,file=''//trim(adjustl(case))//'.dmatup')
open(11,file=''//trim(adjustl(case))//'.dmatdn')

i=0
j=0
ist=0
do 
   i=i+1
   read(10,*,end=20) iatnum,ach,bch,cch
   read(11,*,end=20) iatnum,ach,bch,cch
   if(ach.eq."atom") then
      j=j+1
      ist(j)=i
   end if
end do

20 REWIND(10)
REWIND(11)
ist(j+1)=i

i=0
j=1
do j=1,99
   read(10,*,end=40) iatnum,ach,bch,cch
   read(11,*,end=40) iatnum,ach,bch,cch
   write(6,*) iatnum,trim(ach)," ",trim(bch)," ",trim(cch)," l"," m"
   write(12,*) trim(ach),iatnum
   dpdf=0.0d0
   do i=ist(j)+1,ist(j+1)-1
      read(10,*) are(1),aim,ia,ib,ic,id,ie,if
      read(11,*) are(2),aim,ia,ib,ic,id,ie,if
      if((ic.eq.1.or.ic.eq.2.or.ic.eq.3).and.(ie.eq.if)) then
         dpdf(ic,ie,1)=are(1)
         dpdf(ic,ie,2)=are(2)
!         write(6,*) dpdf(ic,ie,1),dpdf(ic,ie,2),ic,ie
!      else if((ic.eq.2).and.ie.eq.if) then
!         write(6,*) are(1),are(2),ic,ie
      end if
   end do

   pdfnum=0.0d0
   do i=1,3
      do k=-i,i
         pdfnum(i,1)=pdfnum(i,1)+dpdf(i,k,1)
         pdfnum(i,2)=pdfnum(i,2)+dpdf(i,k,2)
      end do
   end do
   
   write(6,*) "occupation p-orbital up, dn, tot"
   write(6,*) pdfnum(1,1),pdfnum(1,2),pdfnum(1,1)+pdfnum(1,2)
   write(6,*) "occupation d-orbital up, dn, tot"
   write(6,*) pdfnum(2,1),pdfnum(2,2),pdfnum(2,1)+pdfnum(2,2)
   write(6,*) "occupation f-orbital up, dn, tot"
   write(6,*) pdfnum(3,1),pdfnum(3,2),pdfnum(3,1)+pdfnum(3,2)

   write(12,*) "ptot=",0.5d0*(dpdf(1,-1,1)+dpdf(1,1,1))-dpdf(1,0,1)+0.5d0*(dpdf(1,-1,2)+dpdf(1,1,2))-dpdf(1,0,2)
   write(12,*) "dtot",dpdf(2,-2,1)+dpdf(2,2,1)-0.5d0*(dpdf(2,-1,1)+dpdf(2,1,1))-dpdf(2,0,1)+dpdf(2,-2,2)+dpdf(2,2,2)-0.5d0*(dpdf(2,-1,2)+dpdf(2,1,2))-dpdf(2,0,2)
   write(6,*) "ptot=",0.5d0*(dpdf(1,-1,1)+dpdf(1,1,1))-dpdf(1,0,1)+0.5d0*(dpdf(1,-1,2)+dpdf(1,1,2))-dpdf(1,0,2)
   write(6,*) "dtot",dpdf(2,-2,1)+dpdf(2,2,1)-0.5d0*(dpdf(2,-1,1)+dpdf(2,1,1))-dpdf(2,0,1)+dpdf(2,-2,2)+dpdf(2,2,2)-0.5d0*(dpdf(2,-1,2)+dpdf(2,1,2))-dpdf(2,0,2)
   write(6,*) "pup=",0.5d0*(dpdf(1,-1,1)+dpdf(1,1,1))-dpdf(1,0,1)
   write(6,*) "pdn=",0.5d0*(dpdf(1,-1,2)+dpdf(1,1,2))-dpdf(1,0,2)
   write(6,*) "dup=",dpdf(2,-2,1)+dpdf(2,2,1)-0.5d0*(dpdf(2,-1,1)+dpdf(2,1,1))-dpdf(2,0,1)
   write(6,*) "ddn=",dpdf(2,-2,2)+dpdf(2,2,2)-0.5d0*(dpdf(2,-1,2)+dpdf(2,1,2))-dpdf(2,0,2)
end do
40 stop

end program qtl

