set terminal postscript eps enhanced color font "Times-Roman"

s=1.01
Tc=555
set ytics -2000,1000,7000
set yrange [-1500:6500]
set mytics 2
set mxtics 2
alpha(x)=(1-s*(x/Tc)**1.5-(1-s)*(x/Tc)**2.5)**(1.0/3.0)
set size 0.5,0.5
set output "englevel_per_alpha.eps"
set nokey
set xlabel "{/Times-Italic B}_{/Times-Roman ex}({/Times-Italic T})/\
{/Times-Italic B}_{/Times-Roman ex}(0)"
set ylabel "Energy Level [K]" 1,0
set label "{/Times-Italic J}=5/2" font ",11" at 0.05,-500
set label "{/Times-Italic J}=7/2" font ",11" at 0.05,1100
set label "{/Times-Italic J}=9/2" font ",11" at 0.05,3000
set label "{/Times-Italic J}=11/2" font ",11" at 0.05,5150
set lmargin 8
plot \
     "SmFe12_lambda411_opencore_englevel_J0_nJex5_rd1.0.txt" u (alpha($1)):2 w l lt 1 lc rgb "black",\
     "SmFe12_lambda411_opencore_englevel_J0_nJex5_rd1.0.txt" u (alpha($1)):3 w l lt 1 lc rgb "black",\
     "SmFe12_lambda411_opencore_englevel_J0_nJex5_rd1.0.txt" u (alpha($1)):4 w l lt 1 lc rgb "black",\
     "SmFe12_lambda411_opencore_englevel_J0_nJex5_rd1.0.txt" u (alpha($1)):5 w l lt 1 lc rgb "black",\
     "SmFe12_lambda411_opencore_englevel_J0_nJex5_rd1.0.txt" u (alpha($1)):6 w l lt 1 lc rgb "black",\
     "SmFe12_lambda411_opencore_englevel_J0_nJex5_rd1.0.txt" u (alpha($1)):7 w l lt 1 lc rgb "black",\
     "SmFe12_lambda411_opencore_englevel_J1_1_nJex5_rd1.0.txt" u (alpha($1)):2 w l lt 1 lc rgb "black",\
     "SmFe12_lambda411_opencore_englevel_J1_1_nJex5_rd1.0.txt" u (alpha($1)):3 w l lt 1 lc rgb "black",\
     "SmFe12_lambda411_opencore_englevel_J1_1_nJex5_rd1.0.txt" u (alpha($1)):4 w l lt 1 lc rgb "black",\
     "SmFe12_lambda411_opencore_englevel_J1_1_nJex5_rd1.0.txt" u (alpha($1)):5 w l lt 1 lc rgb "black",\
     "SmFe12_lambda411_opencore_englevel_J1_1_nJex5_rd1.0.txt" u (alpha($1)):6 w l lt 1 lc rgb "black",\
     "SmFe12_lambda411_opencore_englevel_J1_1_nJex5_rd1.0.txt" u (alpha($1)):7 w l lt 1 lc rgb "black",\
     "SmFe12_lambda411_opencore_englevel_J1_2_nJex5_rd1.0.txt" u (alpha($1)):2 w l lt 1 lc rgb "black",\
     "SmFe12_lambda411_opencore_englevel_J1_2_nJex5_rd1.0.txt" u (alpha($1)):3 w l lt 1 lc rgb "black",\
     "SmFe12_lambda411_opencore_englevel_J2_1_nJex5_rd1.0.txt" u (alpha($1)):2 w l lt 1 lc rgb "black",\
     "SmFe12_lambda411_opencore_englevel_J2_1_nJex5_rd1.0.txt" u (alpha($1)):3 w l lt 1 lc rgb "black",\
     "SmFe12_lambda411_opencore_englevel_J2_1_nJex5_rd1.0.txt" u (alpha($1)):4 w l lt 1 lc rgb "black",\
     "SmFe12_lambda411_opencore_englevel_J2_1_nJex5_rd1.0.txt" u (alpha($1)):5 w l lt 1 lc rgb "black",\
     "SmFe12_lambda411_opencore_englevel_J2_1_nJex5_rd1.0.txt" u (alpha($1)):6 w l lt 1 lc rgb "black",\
     "SmFe12_lambda411_opencore_englevel_J2_1_nJex5_rd1.0.txt" u (alpha($1)):7 w l lt 1 lc rgb "black",\
     "SmFe12_lambda411_opencore_englevel_J2_2_nJex5_rd1.0.txt" u (alpha($1)):2 w l lt 1 lc rgb "black",\
     "SmFe12_lambda411_opencore_englevel_J2_2_nJex5_rd1.0.txt" u (alpha($1)):3 w l lt 1 lc rgb "black",\
     "SmFe12_lambda411_opencore_englevel_J2_2_nJex5_rd1.0.txt" u (alpha($1)):4 w l lt 1 lc rgb "black",\
     "SmFe12_lambda411_opencore_englevel_J2_2_nJex5_rd1.0.txt" u (alpha($1)):5 w l lt 1 lc rgb "black",\
     "SmFe12_lambda411_opencore_englevel_J3_1_nJex5_rd1.0.txt" u (alpha($1)):2 w l lt 1 lc rgb "black",\
     "SmFe12_lambda411_opencore_englevel_J3_1_nJex5_rd1.0.txt" u (alpha($1)):3 w l lt 1 lc rgb "black",\
     "SmFe12_lambda411_opencore_englevel_J3_1_nJex5_rd1.0.txt" u (alpha($1)):4 w l lt 1 lc rgb "black",\
     "SmFe12_lambda411_opencore_englevel_J3_1_nJex5_rd1.0.txt" u (alpha($1)):5 w l lt 1 lc rgb "black",\
     "SmFe12_lambda411_opencore_englevel_J3_1_nJex5_rd1.0.txt" u (alpha($1)):6 w l lt 1 lc rgb "black",\
     "SmFe12_lambda411_opencore_englevel_J3_1_nJex5_rd1.0.txt" u (alpha($1)):7 w l lt 1 lc rgb "black",\
     "SmFe12_lambda411_opencore_englevel_J3_2_nJex5_rd1.0.txt" u (alpha($1)):2 w l lt 1 lc rgb "black",\
     "SmFe12_lambda411_opencore_englevel_J3_2_nJex5_rd1.0.txt" u (alpha($1)):3 w l lt 1 lc rgb "black",\
     "SmFe12_lambda411_opencore_englevel_J3_2_nJex5_rd1.0.txt" u (alpha($1)):4 w l lt 1 lc rgb "black",\
     "SmFe12_lambda411_opencore_englevel_J3_2_nJex5_rd1.0.txt" u (alpha($1)):5 w l lt 1 lc rgb "black",\
     "SmFe12_lambda411_opencore_englevel_J3_2_nJex5_rd1.0.txt" u (alpha($1)):6 w l lt 1 lc rgb "black",\
     "SmFe12_lambda411_opencore_englevel_J3_2_nJex5_rd1.0.txt" u (alpha($1)):7 w l lt 1 lc rgb "black"

reset

set size 0.5,0.5
set output "englevel.eps"
set nokey
plot \
     "SmFe12_lambda411_opencore_englevel_J0_nJex5_rd1.0.txt" u 1:2 w l lt 1 lc rgb "black",\
     "SmFe12_lambda411_opencore_englevel_J0_nJex5_rd1.0.txt" u 1:3 w l lt 1 lc rgb "black",\
     "SmFe12_lambda411_opencore_englevel_J0_nJex5_rd1.0.txt" u 1:4 w l lt 1 lc rgb "black",\
     "SmFe12_lambda411_opencore_englevel_J0_nJex5_rd1.0.txt" u 1:5 w l lt 1 lc rgb "black",\
     "SmFe12_lambda411_opencore_englevel_J0_nJex5_rd1.0.txt" u 1:6 w l lt 1 lc rgb "black",\
     "SmFe12_lambda411_opencore_englevel_J0_nJex5_rd1.0.txt" u 1:7 w l lt 1 lc rgb "black",\
     "SmFe12_lambda411_opencore_englevel_J1_1_nJex5_rd1.0.txt" u 1:2 w l lt 1 lc rgb "black",\
     "SmFe12_lambda411_opencore_englevel_J1_1_nJex5_rd1.0.txt" u 1:3 w l lt 1 lc rgb "black",\
     "SmFe12_lambda411_opencore_englevel_J1_1_nJex5_rd1.0.txt" u 1:4 w l lt 1 lc rgb "black",\
     "SmFe12_lambda411_opencore_englevel_J1_1_nJex5_rd1.0.txt" u 1:5 w l lt 1 lc rgb "black",\
     "SmFe12_lambda411_opencore_englevel_J1_1_nJex5_rd1.0.txt" u 1:6 w l lt 1 lc rgb "black",\
     "SmFe12_lambda411_opencore_englevel_J1_1_nJex5_rd1.0.txt" u 1:7 w l lt 1 lc rgb "black",\
     "SmFe12_lambda411_opencore_englevel_J1_2_nJex5_rd1.0.txt" u 1:2 w l lt 1 lc rgb "black",\
     "SmFe12_lambda411_opencore_englevel_J1_2_nJex5_rd1.0.txt" u 1:3 w l lt 1 lc rgb "black",\
     "SmFe12_lambda411_opencore_englevel_J2_1_nJex5_rd1.0.txt" u 1:2 w l lt 1 lc rgb "black",\
     "SmFe12_lambda411_opencore_englevel_J2_1_nJex5_rd1.0.txt" u 1:3 w l lt 1 lc rgb "black",\
     "SmFe12_lambda411_opencore_englevel_J2_1_nJex5_rd1.0.txt" u 1:4 w l lt 1 lc rgb "black",\
     "SmFe12_lambda411_opencore_englevel_J2_1_nJex5_rd1.0.txt" u 1:5 w l lt 1 lc rgb "black",\
     "SmFe12_lambda411_opencore_englevel_J2_1_nJex5_rd1.0.txt" u 1:6 w l lt 1 lc rgb "black",\
     "SmFe12_lambda411_opencore_englevel_J2_1_nJex5_rd1.0.txt" u 1:7 w l lt 1 lc rgb "black",\
     "SmFe12_lambda411_opencore_englevel_J2_2_nJex5_rd1.0.txt" u 1:2 w l lt 1 lc rgb "black",\
     "SmFe12_lambda411_opencore_englevel_J2_2_nJex5_rd1.0.txt" u 1:3 w l lt 1 lc rgb "black",\
     "SmFe12_lambda411_opencore_englevel_J2_2_nJex5_rd1.0.txt" u 1:4 w l lt 1 lc rgb "black",\
     "SmFe12_lambda411_opencore_englevel_J2_2_nJex5_rd1.0.txt" u 1:5 w l lt 1 lc rgb "black"
reset






set xrange [0:10]
set size 0.5,0.5
set nokey
eng0=0.0
set yrange [-2000:4000]
set ylabel "Energy Level [K]" 1,0
set noxtics
set mytics 2
set xtics ""
set lmargin 8
set label "A" at 2.8,-1600
set label "B" at 3.8,-1600
set label "C" at 4.8,-1600
set label "A" at 6.8,-1600
set label "B" at 7.8,-1600
set label "C" at 8.8,-1600

set output "eng_hierachy_SmFe12.eps"
plot \
     "SmFe12_lambda411_opencore_englevel_T0.0_diag_so_nJex2_rd1.0.txt"       u ($3+0.5):($2-eng0)   w l lw 1 lt 1 lc rgb "black",\
     "SmFe12_lambda411_opencore_englevel_T0.0_diag_soex_nJex2_rd1.0.txt"    u ($3*3+2.5):($2-eng0) w l lw 1 lt 1 lc rgb "black",\
     "SmFe12_lambda411_opencore_englevel_T0.0_diag_soexcf_nJex2_rd1.0.txt" u ($3*3+6.5):($2-eng0) w l lw 1 lt 1 lc rgb "black",\
     "SmFe12_lambda411_opencore_englevel_T0.0_lwstJso_rd1.0.txt" u ($3+0.5):($2-eng0) w l lw 1 lt 1 lc rgb "black",\
     "SmFe12_lambda411_opencore_englevel_T0.0_lwstJsoex_rd1.0.txt" u ($3+2.5):($2-eng0) w l lw 4 lt 1 lc rgb "dark-green",\
     "SmFe12_lambda411_opencore_englevel_T0.0_mixex1_rd1.0.txt" u ($3+3.5):($2-eng0)     w l lw 4 lt 1 lc rgb "blue",\
     "SmFe12_lambda411_opencore_englevel_T0.0_mixex2_rd1.0.txt" u ($3+4.5):($2-eng0)     w l lw 4 lt 1 lc rgb "red",\
     "SmFe12_lambda411_opencore_englevel_T0.0_lwstJsoexcf_rd1.0.txt" u ($3+6.5):($2-eng0) w l lw 4 lt 1 lc rgb "dark-green",\
     "SmFe12_lambda411_opencore_englevel_T0.0_mixexcf1_rd1.0.txt" u ($3+7.5):($2-eng0)     w l lw 4 lt 1 lc rgb "blue",\
     "SmFe12_lambda411_opencore_englevel_T0.0_mixexcf2_rd1.0.txt" u ($3+8.5):($2-eng0)     w l lw 4 lt 1 lc rgb "red"

reset


set terminal postscript eps enhanced color font "Times-Roman"

set size 0.5,0.5
set output "weight.eps"
set nokey
plot \
     "SmFe12_lambda411_opencore_weight_nJex5_rd1.0.txt" u 1:2 w l lt 1 lc rgb "red",\
     "SmFe12_lambda411_opencore_weight_nJex5_rd1.0.txt" u 1:3 w l lt 1 lc rgb "dark-green",\
     "SmFe12_lambda411_opencore_weight_nJex5_rd1.0.txt" u 1:4 w l lt 1 lc rgb "blue"

reset