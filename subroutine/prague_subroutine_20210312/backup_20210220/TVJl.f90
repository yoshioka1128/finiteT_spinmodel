subroutine TVJl(Dex,Dso,J1,L1,S1,Tmix1,Tmix2,HexT,extBJ,coffmix1,coffmix2)
integer(4) :: l,j,ll
real(8) :: Dex,Dso,coffmix1(6,-1:8),coffmix2(6,-1:8)
real(8) :: J1,L1,S1,extBJ(-1:8),Tmix1(6),Tmix2(6),HexT,eigen(100)

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

Tmix1=0.0d0
Tmix2=0.0d0
do l=1,6
   do ll=-1,8
      Tmix1(l)=Tmix1(l)+coffmix1(l,ll)*extBJ(ll)
      Tmix2(l)=Tmix2(l)+coffmix2(l,ll)*extBJ(ll)
   end do
end do

end subroutine TVJl
