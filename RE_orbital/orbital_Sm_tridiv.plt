set terminal postscript eps enhanced color "Times-Roman"

set size 0.5,0.5
#set key font "Times-Roman"
set xlabel "r (Bohr)"
set ylabel "r^2R_{4f}(r)"

set output "orbital.eps"
plot \
"Sm3.2_R4f.dat" w l lw 2 title "trivalent" lc rgb "red",\
"Sm3.2_R4f.dat" w l lw 2 title "divalent" lc rgb "blue"


set title font "Times-Roman"
set output "orbital_Ce_chRmt.eps"
plot \
"Ce3.2_R4f.dat" w l lw 3 lt 1 lc rgb "red" title "rmt=3.2",\
"Ce2.8_R4f.dat" w l lw 2 lt 1 lc rgb "blue" title "rmt=2.8",\
"Ce2.5_R4f.dat" w lp lw 2 lt 2 lc rgb "black" title "rmt=2.5"

reset