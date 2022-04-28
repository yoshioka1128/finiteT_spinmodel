if (exist("n")==0 || n<0) n=nini #変数の初期化
#title(n) = sprintf("time = %d",n)  #Gt*ipace Gt=0.01
title(n) = sprintf("time = %d (*0.4 ps)",n)  #Gt*ipace Gt=0.01

unset label 

set label 1 sprintf("time = %d",n)  font 'Times,20'  at 4 , -12 
set label 2 sprintf("                   (* %.1f ps)", dtime) font 'Times,20'  at 4 , -12 

splot \
sprintf("TDFe%d.txt",n) with vector arrowstyle 3 ,\
sprintf("TDf%d.txt",n) with vector arrowstyle 1 ,\
sprintf("TDg%d.txt",n) with vector arrowstyle 2 #,\

if (n<nend);n=n+dt; reread # loop