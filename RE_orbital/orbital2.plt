set terminal postscript eps enhanced color
set output "orbital2.eps"
set size 0.5,0.5
set xlabel "r (Bohr)"
set ylabel "r^2 R_{4f}(r)"
plot \
"Pr3.2_R4f.dat" w l lt 1 lw 3 lc rgb "red" title "Pr",\
"Nd3.2_R4f.dat" w l lt 1 lw 3 lc rgb "coral" title "Nd",\
"Sm3.2_R4f.dat" w l lt 1 lw 3 lc rgb "dark-green" title "Sm",\
"Dy3.2_R4f.dat" w l lt 1 lw 3 lc rgb "blue" title "Dy",\
"Ho3.2_R4f.dat" w l lt 1 lw 3 lc rgb "brown4" title "Ho"
reset