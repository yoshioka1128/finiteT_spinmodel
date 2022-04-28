set terminal postscript eps enhanced color font "Times-Roman"

set output "magcurve_SmFe12_L5L3.eps"
set size 0.5,0.5
set xrange [0:10]
set yrange [0:2.2]
set ytics 0,0.5,2.5
set mxtics 2
set mytics 2
set xlabel "Applied Field [T]" 
set ylabel "Magnetization [T]" 2,0
set label "300 K" tc rgb "red" at 3.4,1.3 font ",12"
set label "0 K" tc rgb "blue" at 6,1.8 font ",12"
set key reverse Left
set key 9.5,0.8
plot \
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 2 lw 2 lc rgb "blue" title "L=5",\
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T300.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 2 lw 2 lc rgb "red" title "L=5",\
     "SmFe12_L3_opencore_wK1Fe_magcurve_T0.0_100_nJex5_rd1.0.txt" u 1:3 w l lt 1 lw 2 lc rgb "blue" title "L=3",\
     "SmFe12_L3_opencore_wK1Fe_magcurve_T300.0_100_nJex5_rd1.0.txt" u 1:3 w l lt 1 lw 2 lc rgb "red" title "L=3"
reset

