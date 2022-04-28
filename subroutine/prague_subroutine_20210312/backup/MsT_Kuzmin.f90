function MsT_Kuzmin(temp,Ms,Tc,s,p)
implicit none
real(8) :: temp,s,p,Tc,Ms,MsT_Kuzmin

MsT_Kuzmin=Ms*(1.0d0-s*(temp/Tc)**(1.5d0)-(1.0d0-s)*(temp/Tc)**p)**(1.0d0/3.0d0)
if(MsT_Kuzmin.lt.0.0d0) MsT_Kuzmin=0.0d0

end function MsT_Kuzmin
