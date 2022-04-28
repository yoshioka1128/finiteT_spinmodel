subroutine CkqJ(J1,CJst,JM)
implicit none
integer(4),parameter :: nJmax=21
real(8),parameter :: pi=dacos(-1.0d0)
complex(kind(0d0)),parameter :: im=(0.0d0,1.0d0)
integer(4) :: i,ii,j, ll, mm
real(8) :: J1, dl, abc,def,w3jsym,w6jsym, CJst(nJmax,-1:8,-8:8),JM(100,2)

! matrix element of spherical tensor operator
CJst=0.0d0
do ll=0,8
   do mm=-8,8
      dl=dble(ll)
      abc=1.0d0
      do i=int(2.0d0*J1+dl+1.0d0),int(2.0d0*J1-dl+1.0d0),-1
         abc=abc*dble(i)
      end do
      do ii=1,int(2.0d0*J1+1)
         CJst(ii,ll,mm)=(-1.0d0)**(int(J1-JM(ii,2)))*dsqrt(abc)/(2.0d0**ll)*w3jsym(J1,dl,J1,-JM(ii,2),dble(mm))
      end do
   end do
end do

do mm=-8,8
   do ii=1,int(2.0d0*J1+1)
      CJst(ii,-1,mm)=0.0d0
   end do
end do

return

end subroutine CkqJ
