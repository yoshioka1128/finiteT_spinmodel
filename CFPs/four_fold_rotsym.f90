program four_fold_rotsym
implicit none
integer(4) :: ir,i,l,m,ncfp,nr,nfe
real(8) :: Alm(6,-6:6),abc,rd,rl(6),Hex
character(20) :: ratom,tatom
open(33,file="~/research/CFPs/cfpSmFe11V_j_full_opencore.txt")
read(33,*) ratom,nr,tatom,nfe
read(33,*) 
read(33,*) 
read(33,*) 
read(33,*) 
read(33,*) 
Alm=0.0d0
read(33,*) ratom
read(33,*) Hex,rd ! Hex
read(33,*) rl(2),rl(4),rl(6),ncfp
do i=1,ncfp
   read(33,*) abc,l,m
   Alm(l,m)=abc
!   write(6,*) Alm(l,m),l,m
end do


write(6,*) ratom,2
write(6,"(2f10.5)") Hex,rd ! Hex
write(6,"(3f10.5,i5)") rl(2),rl(4),rl(6),ncfp
do l=2,6,2
   write(6,*) Alm(l,0),l,0
   do m=1,l
      if(m.eq.1) then
         if(Alm(l,m).ne.0.0d0) write(6,*) -Alm(l,m),l,-m
         if(Alm(l,-m).ne.0.0d0) write(6,*) Alm(l,-m),l,m
      else if(m.eq.2) then
         if(Alm(l,m).ne.0.0d0) write(6,*) -Alm(l,m),l,m
         if(Alm(l,-m).ne.0.0d0) write(6,*) -Alm(l,-m),l,-m
      else if(m.eq.3) then
         if(Alm(l,m).ne.0.0d0) write(6,*) Alm(l,m),l,-m
         if(Alm(l,-m).ne.0.0d0) write(6,*) -Alm(l,-m),l,m
      else if(m.eq.4) then
         if(Alm(l,m).ne.0.0d0) write(6,*) Alm(l,m),l,m
         if(Alm(l,-m).ne.0.0d0) write(6,*) Alm(l,-m),l,-m
      else if(m.eq.5) then
         if(Alm(l,m).ne.0.0d0) write(6,*) -Alm(l,m),l,-m
         if(Alm(l,-m).ne.0.0d0) write(6,*) Alm(l,-m),l,m
      else if(m.eq.6) then
         if(Alm(l,m).ne.0.0d0) write(6,*) -Alm(l,m),l,m
         if(Alm(l,-m).ne.0.0d0) write(6,*) -Alm(l,-m),l,-m
      end if
   end do
end do


write(6,*) ratom,3
write(6,"(2f10.5)") Hex,rd ! Hex
write(6,"(3f10.5,i5)") rl(2),rl(4),rl(6),ncfp
do l=2,6,2
   write(6,*) Alm(l,0),l,0
   do m=1,l
      if(m.eq.1) then
         if(Alm(l,m).ne.0.0d0) write(6,*) -Alm(l,m),l,m
         if(Alm(l,-m).ne.0.0d0) write(6,*) -Alm(l,-m),l,-m
      else if(m.eq.2) then
         if(Alm(l,m).ne.0.0d0) write(6,*) Alm(l,m),l,m
         if(Alm(l,-m).ne.0.0d0) write(6,*) Alm(l,-m),l,-m
      else if(m.eq.3) then
         if(Alm(l,m).ne.0.0d0) write(6,*) -Alm(l,m),l,m
         if(Alm(l,-m).ne.0.0d0) write(6,*) -Alm(l,-m),l,-m
      else if(m.eq.4) then
         if(Alm(l,m).ne.0.0d0) write(6,*) Alm(l,m),l,m
         if(Alm(l,-m).ne.0.0d0) write(6,*) Alm(l,-m),l,-m
      else if(m.eq.5) then
         if(Alm(l,m).ne.0.0d0) write(6,*) -Alm(l,m),l,m
         if(Alm(l,-m).ne.0.0d0) write(6,*) -Alm(l,-m),l,-m
      else if(m.eq.6) then
         if(Alm(l,m).ne.0.0d0) write(6,*) Alm(l,m),l,m
         if(Alm(l,-m).ne.0.0d0) write(6,*) Alm(l,-m),l,-m
      end if
   end do
end do


write(6,*) ratom,4
write(6,"(2f10.5)") Hex,rd ! Hex
write(6,"(3f10.5,i5)") rl(2),rl(4),rl(6),ncfp
do l=2,6,2
   write(6,*) Alm(l,0),l,0
   do m=1,l
      if(m.eq.1) then
         if(Alm(l,m).ne.0.0d0) write(6,*) Alm(l,m),l,-m
         if(Alm(l,-m).ne.0.0d0) write(6,*) -Alm(l,-m),l,m
      else if(m.eq.2) then
         if(Alm(l,m).ne.0.0d0) write(6,*) -Alm(l,m),l,m
         if(Alm(l,-m).ne.0.0d0) write(6,*) -Alm(l,-m),l,-m
      else if(m.eq.3) then
         if(Alm(l,m).ne.0.0d0) write(6,*) -Alm(l,m),l,-m
         if(Alm(l,-m).ne.0.0d0) write(6,*) Alm(l,-m),l,m
      else if(m.eq.4) then
         if(Alm(l,m).ne.0.0d0) write(6,*) Alm(l,m),l,m
         if(Alm(l,-m).ne.0.0d0) write(6,*) Alm(l,-m),l,-m
      else if(m.eq.5) then
         if(Alm(l,m).ne.0.0d0) write(6,*) Alm(l,m),l,-m
         if(Alm(l,-m).ne.0.0d0) write(6,*) -Alm(l,-m),l,m
      else if(m.eq.6) then
         if(Alm(l,m).ne.0.0d0) write(6,*) -Alm(l,m),l,m
         if(Alm(l,-m).ne.0.0d0) write(6,*) -Alm(l,-m),l,-m
      end if
   end do
end do


end program four_fold_rotsym
