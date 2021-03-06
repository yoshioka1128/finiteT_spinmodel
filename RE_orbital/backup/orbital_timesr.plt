set terminal postscript eps enhanced
set output "orbital*r^2.eps"
set size 0.5,0.5
plot \
"Pr3.2_R4f.dat" using 1:($1**2*$2) w l lw 3 ,\
"Nd3.2_R4f.dat" using 1:($1**2*$2) w l lw 3,\
"Dy3.2_R4f.dat" using 1:($1**2*$2) w l lw 3,\
"Ho3.2_R4f.dat" using 1:($1**2*$2) w l lw 3
#"Tb3.2_R4f.dat" using 1:($1**2*$2) w l,\
reset