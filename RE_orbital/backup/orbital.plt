set terminal postscript eps enhanced color
set output "orbital.eps"
set size 0.5,0.5
unset key
plot \
"Pr3.2_R4f.dat" w l lw 3 ,\
"Nd3.2_R4f.dat" w l lw 3,\
"Dy3.2_R4f.dat" w l lw 3,\
"Ho3.2_R4f.dat" w l lw 3
reset