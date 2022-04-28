set terminal postscript eps enhanced color font "Times-Roman,18"

kb=13.80649
VV=8.5205*8.5205*4.7693
VTi=8.54*8.54*4.78
KuFeVT0=47.65712/2.0
KuFeTiT0=47.65712/2.0
KuFeVT400=18.0218286652134/2.0
KuFeTiT400=14.6539768899952/2.0

set xrange [0:180]
set yrange [-40:220]
set mytics 2
ysize=0.8
set size 0.65,ysize
tmgn=0.98*ysize
bmgn=0.16*ysize
lmgn=0.105
rmgn=0.56
size=(tmgn-bmgn)/2.0
set output "SmFe11V_i_T0T400_eng_phi_theta90.eps"
set multiplot
set lmargin at screen lmgn
set rmargin at screen rmgn
set tmargin at screen tmgn
set bmargin at screen bmgn+size
set format x ""
set xtics 0,45,180
set ytics nomirror
set y2tics
set my2tics 2
set y2range [-40*2*kb/VV:220*2*kb/VV]
set y2tics offset -0.5,0
set label "(a) {/Times-Italic T}=0 K" at 8,195 font ",20"
set ylabel "MA Energy [K/f.u.]" 0.7,-5 font ",20"
set y2label "MA Energy [MJ/m^3]" -1.5,-5 font ",20"
#set arrow from 82.5,25 to 97.5,25 size graph 0.02,30,90 lw 1 lc rgb "dark-violet" 
set tics font ",20"
unset key
plot \
      0 w l lt 1 lw 1 lc rgb "black" notitle,\
     "SmFe11V_i_fullexp_lambda411_opencore_T0.0_eng_phi_theta90_nJex2_rd1.0.txt" u (180*$1/pi):($2+KuFeVT0) w l lt 1 lw 5 lc rgb "red" title "Sm_1",\
     "SmFe11V_i_fullexp_lambda411_opencore_T0.0_eng_phi_theta90_nJex2_rd1.0.txt" u (180*$1/pi):($3+KuFeVT0) w l lt 1 lw 5 lc rgb "blue" title "Sm_2",\
     "SmFe11V_i_fullexp_lambda411_opencore_T0.0_eng_phi_theta90_nJex2_rd1.0.txt" u (180*$1/pi):(($2+$3)/2.0+KuFeVT0) w l lt 1 lw 5 lc rgb "dark-violet" title "average",\
     "SmFe11V_i_fullexp_lambda411_opencore_T0.0_eng_an_phi_theta90_rd1.0.txt" u (180*$1/pi):($2+KuFeVT0) w l lt 2 lw 2.5 lc rgb "red" notitle "Sm 1",\
     "SmFe11V_i_fullexp_lambda411_opencore_T0.0_eng_an_phi_theta90_rd1.0.txt" u (180*$1/pi):($3+KuFeVT0) w l lt 2 lw 2.5 lc rgb "blue" notitle "Sm 2",\
     "SmFe11V_i_fullexp_lambda411_opencore_T0.0_eng_an_phi_theta90_rd1.0.txt" u (180*$1/pi):(($2+$3)/2.0+KuFeVT0) w l lt 2 lw 2.5 lc rgb "dark-violet" notitle "m0 and m4"
unset arrow
unset ylabel
unset y2label
unset label
set xtics 0,90,180
set xtics ("0"0,"{/Symbol p}/4"45,"{/Symbol p}/2"90,"3{/Symbol p}/4"135,"{/Symbol p}"180)
set xlabel "{/Symbol-Oblique f}^{TM} ({/Symbol-Oblique q}^{TM}={/Symbol p}/2)" 0.0,-0.32 font ",20"
set lmargin at screen lmgn
set rmargin at screen rmgn
set tmargin at screen bmgn+size
set bmargin at screen bmgn
set mxtics 2
set label "(b) {/Times-Italic T}=400 K" at 8,195 font ",20"
set key reverse Left
set key 185,200
#set arrow from 82.5,12 to 97.5,12 size graph 0.02,30,90 lw 1 lc rgb "dark-violet" 
plot \
      0 w l lt 1 lw 1 lc rgb "black" notitle,\
     "SmFe11V_i_fullexp_lambda411_opencore_T400.0_eng_phi_theta90_nJex2_rd1.0.txt" u (180*$1/pi):($2+KuFeVT400) w l lt 1 lw 5 lc rgb "red" title "Sm_1",\
     "SmFe11V_i_fullexp_lambda411_opencore_T400.0_eng_phi_theta90_nJex2_rd1.0.txt" u (180*$1/pi):($3+KuFeVT400) w l lt 1 lw 5 lc rgb "blue" title "Sm_2",\
     "SmFe11V_i_fullexp_lambda411_opencore_T400.0_eng_phi_theta90_nJex2_rd1.0.txt" u (180*$1/pi):(($2+$3)/2.0+KuFeVT400) w l lt 1 lw 5 lc rgb "dark-violet" title "average",\
     "SmFe11V_i_fullexp_lambda411_opencore_T400.0_eng_an_phi_theta90_rd1.0.txt" u (180*$1/pi):($2+KuFeVT400) w l lt 2 lw 2.5 lc rgb "red" notitle "Sm 1",\
     "SmFe11V_i_fullexp_lambda411_opencore_T400.0_eng_an_phi_theta90_rd1.0.txt" u (180*$1/pi):($3+KuFeVT400) w l lt 2 lw 2.5 lc rgb "blue" notitle "Sm 2",\
     "SmFe11V_i_fullexp_lambda411_opencore_T400.0_eng_an_phi_theta90_rd1.0.txt" u (180*$1/pi):(($2+$3)/2.0+KuFeVT400) w l lt 2 lw 2.5 lc rgb "dark-violet" notitle "m0 and m4"
#     "SmFe11V_i_fullexp_lambda411_opencore_T400.0_eng_an_phi_theta90_rd1.0.txt" u (180*$1/pi):($4*VV/kb/2.0) w l lt 2 lw 2.5 lc rgb "dark-violet" notitle "m0 and m4"

unset multiplot
reset







set xrange [0:180]
set yrange [-140:140]
set mytics 2
ysize=0.8
set size 0.65,ysize
tmgn=0.98*ysize
bmgn=0.16*ysize
lmgn=0.105
rmgn=0.56
size=(tmgn-bmgn)/2.0
set output "SmFe11Ti_i_T0T400_eng_phi_theta90.eps"
set multiplot
set lmargin at screen lmgn
set rmargin at screen rmgn
set tmargin at screen tmgn
set bmargin at screen bmgn+size
set format x ""
set xtics 0,45,180
set ytics nomirror
set y2tics
set my2tics 2
set y2range [-140*2*kb/VTi:140*2*kb/VTi]
set y2tics offset -0.5,0
set label "(a) {/Times-Italic T}=0 K" at 20,110 font ",20"
set ylabel "MA Energy [K/f.u.]" 0.7,-5 font ",20"
set y2label "MA Energy [MJ/m^3]" -1.5,-5 font ",20"
#set arrow from 82.5,25 to 97.5,25 size graph 0.02,30,90 lw 1 lc rgb "dark-violet" 
set tics font ",20"
unset key
plot \
      0 w l lt 1 lw 1 lc rgb "black" notitle,\
     "SmFe11Ti_i_fullexp_lambda411_opencore_T0.0_eng_phi_theta90_nJex2_rd1.0.txt" u (180*$1/pi):($2+KuFeTiT0) w l lt 1 lw 5 lc rgb "red" title "Sm_1",\
     "SmFe11Ti_i_fullexp_lambda411_opencore_T0.0_eng_phi_theta90_nJex2_rd1.0.txt" u (180*$1/pi):($3+KuFeTiT0) w l lt 1 lw 5 lc rgb "blue" title "Sm_2",\
     "SmFe11Ti_i_fullexp_lambda411_opencore_T0.0_eng_phi_theta90_nJex2_rd1.0.txt" u (180*$1/pi):($4*VTi/kb/2.0) w l lt 1 lw 5 lc rgb "dark-violet" title "average",\
     "SmFe11Ti_i_fullexp_lambda411_opencore_T0.0_eng_an_phi_theta90_rd1.0.txt" u (180*$1/pi):($2+KuFeTiT0) w l lt 2 lw 2.5 lc rgb "red" notitle "Sm 1",\
     "SmFe11Ti_i_fullexp_lambda411_opencore_T0.0_eng_an_phi_theta90_rd1.0.txt" u (180*$1/pi):($3+KuFeTiT0) w l lt 2 lw 2.5 lc rgb "blue" notitle "Sm 2",\
     "SmFe11Ti_i_fullexp_lambda411_opencore_T0.0_eng_an_phi_theta90_rd1.0.txt" u (180*$1/pi):($4*VTi/kb/2.0) w l lt 2 lw 2.5 lc rgb "dark-violet" notitle "m0 and m4"
unset arrow
unset ylabel
unset y2label
unset label
set xtics ("0"0,"{/Symbol p}/4"45,"{/Symbol p}/2"90,"3{/Symbol p}/4"135,"{/Symbol p}"180)
set xlabel "{/Symbol-Oblique f}^{TM} ({/Symbol-Oblique q}^{TM}={/Symbol p}/2)" 0.0,-0.32 font ",20"
set lmargin at screen lmgn
set rmargin at screen rmgn
set tmargin at screen bmgn+size
set bmargin at screen bmgn
set format x
set label "(b) {/Times-Italic T}=400 K" at 5,110 font ",20"
set key reverse Left
set key 185,-35
#set arrow from 82.5,12 to 97.5,12 size graph 0.02,30,90 lw 1 lc rgb "dark-violet" 
plot \
      0 w l lt 1 lw 1 lc rgb "black" notitle,\
     "SmFe11Ti_i_fullexp_lambda411_opencore_T400.0_eng_phi_theta90_nJex2_rd1.0.txt" u (180*$1/pi):($2+KuFeTiT400) w l lt 1 lw 5 lc rgb "red" title "Sm_1",\
     "SmFe11Ti_i_fullexp_lambda411_opencore_T400.0_eng_phi_theta90_nJex2_rd1.0.txt" u (180*$1/pi):($3+KuFeTiT400) w l lt 1 lw 5 lc rgb "blue" title "Sm_2",\
     "SmFe11Ti_i_fullexp_lambda411_opencore_T400.0_eng_phi_theta90_nJex2_rd1.0.txt" u (180*$1/pi):(($2+$3)/2.0+KuFeTiT400) w l lt 1 lw 5 lc rgb "dark-violet" title "average",\
     "SmFe11Ti_i_fullexp_lambda411_opencore_T400.0_eng_an_phi_theta90_rd1.0.txt" u (180*$1/pi):($2+KuFeTiT400) w l lt 2 lw 2.5 lc rgb "red" notitle "Sm 1",\
     "SmFe11Ti_i_fullexp_lambda411_opencore_T400.0_eng_an_phi_theta90_rd1.0.txt" u (180*$1/pi):($3+KuFeTiT400) w l lt 2 lw 2.5 lc rgb "blue" notitle "Sm 2",\
     "SmFe11Ti_i_fullexp_lambda411_opencore_T400.0_eng_an_phi_theta90_rd1.0.txt" u (180*$1/pi):(($2+$3)/2.0+KuFeTiT400) w l lt 2 lw 2.5 lc rgb "dark-violet" notitle "m0 and m4"
unset multiplot
reset









kb=13.80649
VV=8.5205*8.5205*4.7693
KuFe=47.65712/2.0
set xrange [0:180]
set yrange [-148:148]
set mxtics 2
set mytics 2
set output "SmFe11V_i_eng_dphi_phi_theta90.eps"
ysize=0.8
set size 0.5,ysize
tmgn=0.95*ysize
bmgn=0.12*ysize
lmgn=0.08
rmgn=0.44
size=(tmgn-bmgn)/2.0
set multiplot
set lmargin at screen lmgn
set rmargin at screen rmgn
set tmargin at screen tmgn
set bmargin at screen bmgn+size
set format x ""
set xtics 0,45,180
set ytics nomirror
set y2tics
set my2tics 2
set y2range [-148*2*kb/VV:148*2*kb/VV]
set y2tics offset -0.5,0
set label "(a)" at 5,135
set ylabel "MA Energy [K/atom]" 0.9,0
set y2label "MA Energy [MJ/m^3]" -1.7,0
set key reverse Left
set key 163,-69
plot \
      0 w l lt 1 lw 1 lc rgb "black" notitle,\
     "SmFe11V_i_fullexp_lambda411_opencore_T0.0_eng_phi_theta90_nJex2_rd1.0.txt" u (180*$1/pi):($2+KuFe) w l lt 1 lw 3 lc rgb "red" title "Sm_1",\
     "SmFe11V_i_fullexp_lambda411_opencore_T0.0_eng_phi_theta90_nJex2_rd1.0.txt" u (180*$1/pi):($3+KuFe) w l lt 1 lw 3 lc rgb "blue" title "Sm_2",\
     "SmFe11V_i_fullexp_lambda411_opencore_T0.0_eng_phi_theta90_nJex2_rd1.0.txt" u (180*$1/pi):(($2+$3)/2.0+KuFe) w l lt 1 lw 3 lc rgb "dark-violet" title "ave.",\
     "SmFe11V_i_fullexp_lambda411_opencore_T0.0_eng_an_phi_theta90_rd1.0.txt" u (180*$1/pi):($2+KuFe) w l lt 2 lw 1.5 lc rgb "red" notitle "Sm 1",\
     "SmFe11V_i_fullexp_lambda411_opencore_T0.0_eng_an_phi_theta90_rd1.0.txt" u (180*$1/pi):($3+KuFe) w l lt 2 lw 1.5 lc rgb "blue" notitle "Sm 2",\
     "SmFe11V_i_fullexp_lambda411_opencore_T0.0_eng_an_phi_theta90_rd1.0.txt" u (180*$1/pi):(($2+$3)/2.0+KuFe) w l lt 2 lw 1.5 lc rgb "dark-violet" notitle "m0 and m4"
unset key
unset y2tics
unset ylabel
unset y2label
set ytics mirror
set key reverse Left
set key 160,-75
unset label
set xlabel "{/Symbol-Oblique f}^{TM} [deg.] ({/Symbol-Oblique q}^{TM}=0)"
set lmargin at screen lmgn
set rmargin at screen rmgn
set tmargin at screen bmgn+size
set bmargin at screen bmgn
set format x
set yrange [-15:15]
set ylabel "Angular Difference [deg.]" 0,0
#set ylabel "{/Symbol D} {/Symbol-Oblique f}_{/Times-Italic S, L} [deg.]" 0,0
set label "(b)" at 5,13.5
set label "{/Symbol D}{/Symbol-Oblique f}_{{/Times-Italic S}, 1}" at 39,-6.3 tc rgb "red"
set label "{/Symbol D}{/Symbol-Oblique f}_{{/Times-Italic L}, 1}" at 47,-12.5 tc rgb "red"
set label "{/Symbol D}{/Symbol-Oblique f}_{{/Times-Italic S}, 2}" at 129,-6.3  tc rgb "blue"
set label "{/Symbol D}{/Symbol-Oblique f}_{{/Times-Italic L}, 2}" at 137,-12.5 tc rgb "blue"

plot \
      0 w l lt 1 lw 1 lc rgb "black" notitle,\
     "SmFe11V_i_fullexp_lambda411_opencore_T0.0_SLangle_phi_theta90_nJex2_rd1.0.txt" u (180*$1/pi):(180*($2-$1-pi)/pi) w l lw 3 lt 5 lc rgb "red" title "Sm 1 S",\
     "SmFe11V_i_fullexp_lambda411_opencore_T0.0_SLangle_phi_theta90_nJex2_rd1.0.txt" u (180*$1/pi):(180*($4-$1)/pi) w l lw 3 lt 1 lc rgb "red" title "Sm 1 L",\
     "SmFe11V_i_fullexp_lambda411_opencore_T0.0_SLangle_phi_theta90_nJex2_rd1.0.txt" u (180*$1/pi):(180*($3-$1-pi)/pi) w l lw 3 lt 5 lc rgb "blue" title "Sm 2 S",\
     "SmFe11V_i_fullexp_lambda411_opencore_T0.0_SLangle_phi_theta90_nJex2_rd1.0.txt" u (180*$1/pi):(180*($5-$1)/pi) w l lw 3 lt 1 lc rgb "blue" title "Sm 2 L"

unset multiplot
reset










K1num_Ti_i=6.9852653885
K2num_Ti_i=-5.9322741754
K3num_Ti_i=-2.5493092538
K1num_V_i=7.1267953406
K2num_V_i=-6.5338671789
K3num_V_i=-3.1416496487
K1num_Co_i=11.9234407953
K2num_Co_i=-0.6949271351
K3num_Co_i=-0.6517892882

K1num_Ti_j=14.5356312529
K2num_Ti_j=1.2097069144
K3num_Ti_j=0.2941871546
K1num_V_j=16.7603995326
K2num_V_j=2.4883705038
K3num_V_j=0.9033540414
K1num_Co_j=11.4994706740
K2num_Co_j=-1.3981713182
K3num_Co_j=-1.1350681540

K1an_Ti_i=5.80117897099174       
K2an_Ti_i=-4.82940477036506     
K3an_Ti_i=0.334017673822928
K1an_V_i=5.99278667407654      
K2an_V_i=-5.13294605774916     
K3an_V_i=0.356491816659877     
K1an_Co_i=12.4916392524232       
K2an_Co_i=-2.85691406853201     
K3an_Co_i=0.324937495879319     

K1an_Ti_j=17.1367454132188       
K2an_Ti_j=-2.27228431300856     
K3an_Ti_j=0.357896373593641     
K1an_V_j=21.4738256883728       
K2an_V_j=-1.87984624076897     
K3an_V_j=0.379631632406549     
K1an_Co_j=11.5697450957858       
K2an_Co_j=-3.49779596752900     
K3an_Co_j=0.315954847677776     

kb=13.80649
set nokey
set xrange [0:90]
set mxtics 2
set mytics 2
set xtics 0,45,180
set size 0.65,1.0

set output "SmFe11TiVCo_i_eng_theta_phi22.5.eps"
set multiplot
set yrange [0:11]

lmgn=0.12
rmgn=0.55
bmgn=0.12
tmgn=0.98
xsize=(rmgn-lmgn)
ysize=(tmgn-bmgn)/3.0
# set label position
lblx=5
lbly=9.5
set tmargin at screen bmgn+3*ysize
set bmargin at screen bmgn+2*ysize
set lmargin at screen lmgn
set rmargin at screen lmgn+xsize
set format x ""
set label "(a) SmFe_{11}Ti (8{/Times-Italic i} site)" at lblx,lbly
set nokey
plot \
      0 w l lt 1 lw 1 lc rgb "black" notitle,\
     "SmFe11Ti_i_fullexp_lambda411_opencore_T0.0_eng_theta_phi22.5_nJex2_rd1.0.txt" u (180*$1/pi):($4) w l lt 1 lw 3 lc rgb "black" notitle "ave.",\
     K1an_Ti_i*sin(pi*x/180)**2+K2an_Ti_i*sin(pi*x/180)**4+K3an_Ti_i*sin(pi*x/180)**6 w l lt 1 lw 3 lc rgb "red",\
     K1num_Ti_i*sin(pi*x/180)**2+K2num_Ti_i*sin(pi*x/180)**4+K3num_Ti_i*sin(pi*x/180)**6 w l lw 3 lc rgb "blue",\
     "SmFe11Ti_i_fullexp_lambda411_opencore_T0.0_eng_an_theta_phi22.5_rd1.0.txt" u (180*$1/pi):($4) w l lt 2 lw 3 lc rgb "violet" notitle ""
unset label
unset key
set ylabel "Magnetic Anisotropy Energy [MJ/m^3]" 0.5,0
set format y
set tmargin at screen bmgn+2*ysize
set bmargin at screen bmgn+ysize
set lmargin at screen lmgn
set rmargin at screen lmgn+xsize
set label "(b) SmFe_{11}V (8{/Times-Italic i} site)" at lblx,lbly
plot \
     0 w l lt 1 lw 1 lc rgb "black" notitle,\
     "SmFe11V_i_fullexp_lambda411_opencore_T0.0_eng_theta_phi22.5_nJex2_rd1.0.txt" u (180*$1/pi):($4) w l lt 1 lw 3 lc rgb "black" title "ave.",\
     K1an_V_i*sin(pi*x/180)**2+K2an_V_i*sin(pi*x/180)**4+K3an_V_i*sin(pi*x/180)**6 w l lt 1 lw 3 lc rgb "red",\
     K1num_V_i*sin(pi*x/180)**2+K2num_V_i*sin(pi*x/180)**4+K3num_V_i*sin(pi*x/180)**6 w l lw 3 lc rgb "blue",\
     "SmFe11V_i_fullexp_lambda411_opencore_T0.0_eng_an_theta_phi22.5_rd1.0.txt" u (180*$1/pi):($4) w l lt 2 lw 3 lc rgb "violet" notitle

unset label
unset ylabel
set format x
set format y
set xtics
set tmargin at screen bmgn+ysize
set bmargin at screen bmgn
set lmargin at screen lmgn
set rmargin at screen lmgn+xsize
set label "(c) SmFe_{11}Co (8{/Times-Italic i} site)" at lblx,lbly
set xlabel "{/Symbol-Oblique q}^{TM} [degree]" 0.0,0
set key reverse Left
set key right bottom
plot \
     0 w l lt 1 lw 1 lc rgb "black" notitle,\
     "SmFe11Co_i_fullexp_lambda411_opencore_T0.0_eng_theta_phi22.5_nJex2_rd1.0.txt" u (180*$1/pi):($4) w l lt 1 lw 3 lc rgb "black" notitle "ave.",\
     K1an_Co_i*sin(pi*x/180)**2+K2an_Co_i*sin(pi*x/180)**4+K3an_Co_i*sin(pi*x/180)**6 w l lt 1 lw 3 lc rgb "red" title "analytical",\
     K1num_Co_i*sin(pi*x/180)**2+K2num_Co_i*sin(pi*x/180)**4+K3num_Co_i*sin(pi*x/180)**6 w l lw 3 lc rgb "blue" title "numerical",\
     "SmFe11Co_i_fullexp_lambda411_opencore_T0.0_eng_an_theta_phi22.5_rd1.0.txt" u (180*$1/pi):($4) w l lt 2 lw 3 lc rgb "violet" title "analytical2"
unset multiplot
unset xlabel
unset yrange
unset label
reset






set terminal postscript eps enhanced color font "Times-Roman,12"
set size 0.5,0.5
set xrange [0:90]
set xtics 0,45,90
set mxtics 2
set mytics 2
set nokey
set output "SmFe12H_eng_theta_phi22.5.eps"

K1num_H=16.8194332782
K2num_H=0.5643342675
K3num_H=-1.0229105321
K1an_H=20.126122207709546
K2an_H=-6.891679198387643
K3an_H=0.232902698751507
set yrange [0:18]
set xlabel "{/Symbol-Oblique q}^{TM} [degree]" 0.0,0
set ylabel "MA Energy [MJ/m^3]" 0.5,0
set key reverse Left
set key left top
plot \
      0 w l lt 1 lw 1 lc rgb "black" notitle,\
     "SmFe12H_lambda411_opencore_T0.0_eng_theta_phi22.5_nJex2_rd1.0.txt" u (180*$1/pi):($4) w l lt 1 lw 3 lc rgb "black" title "statistical",\
     K1an_H*sin(pi*x/180)**2+K2an_H*sin(pi*x/180)**4+K3an_H*sin(pi*x/180)**6 w l lt 1 lw 3 lc rgb "red" title "Ki(analytical)",\
     K1num_H*sin(pi*x/180)**2+K2num_H*sin(pi*x/180)**4+K3num_H*sin(pi*x/180)**6 w l lw 3 lc rgb "blue" title "Ki(numerical)",\
     "SmFe12H_lambda411_opencore_T0.0_eng_an_theta_phi22.5_rd1.0.txt" u (180*$1/pi):4 w l lt 2 lw 3 lc rgb "violet" title "Ki(numerical2)"



set yrange [0:24]
set output "SmFe11TiVCo_j_eng_theta_phi22.5.eps"
set multiplot

lmgn=0.12
rmgn=0.55
bmgn=0.12
tmgn=0.98
xsize=(rmgn-lmgn)
ysize=(tmgn-bmgn)/3.0
# set label position
lblx=5
lbly=21
set tmargin at screen bmgn+3*ysize
set bmargin at screen bmgn+2*ysize
set lmargin at screen lmgn
set rmargin at screen lmgn+xsize
set format x ""
set label "(a) SmFe_{11}Ti (8{/Times-Italic i} site)" at lblx,lbly
set nokey
plot \
      0 w l lt 1 lw 1 lc rgb "black" notitle,\
     "SmFe11Ti_j_fullexp_lambda411_opencore_T0.0_eng_theta_phi22.5_nJex2_rd1.0.txt" u (180*$1/pi):($4) w l lt 1 lw 3 lc rgb "black" notitle "ave.",\
     K1an_Ti_j*sin(pi*x/180)**2+K2an_Ti_j*sin(pi*x/180)**4+K3an_Ti_j*sin(pi*x/180)**6 w l lw 3 lc rgb "red",\
     K1num_Ti_j*sin(pi*x/180)**2+K2num_Ti_j*sin(pi*x/180)**4+K3num_Ti_j*sin(pi*x/180)**6 w l lw 3 lc rgb "blue"
unset label
unset key
set ylabel "Magnetic Anisotropy Energy [MJ/m^3]" 0.5,0
set format y
set tmargin at screen bmgn+2*ysize
set bmargin at screen bmgn+ysize
set lmargin at screen lmgn
set rmargin at screen lmgn+xsize
set label "(b) SmFe_{11}V (8{/Times-Italic i} site)" at lblx,lbly
plot \
     0 w l lt 1 lw 1 lc rgb "black" notitle,\
     "SmFe11V_j_fullexp_lambda411_opencore_T0.0_eng_theta_phi22.5_nJex2_rd1.0.txt" u (180*$1/pi):($4) w l lt 1 lw 3 lc rgb "black" title "ave.",\
     K1an_V_j*sin(pi*x/180)**2+K2an_V_j*sin(pi*x/180)**4+K3an_V_j*sin(pi*x/180)**6 w l lw 3 lc rgb "red",\
     K1num_V_j*sin(pi*x/180)**2+K2num_V_j*sin(pi*x/180)**4+K3num_V_j*sin(pi*x/180)**6 w l lw 3 lc rgb "blue"
unset label
unset ylabel
set format x
set format y
set xtics
set tmargin at screen bmgn+ysize
set bmargin at screen bmgn
set lmargin at screen lmgn
set rmargin at screen lmgn+xsize
set label "(c) SmFe_{11}Co (8{/Times-Italic i} site)" at lblx,lbly
set xlabel "{/Symbol-Oblique q}^{TM} [degree]" 0.0,0
set key reverse Left
set key left center
plot \
     0 w l lt 1 lw 1 lc rgb "black" notitle,\
     "SmFe11Co_j_fullexp_lambda411_opencore_T0.0_eng_theta_phi22.5_nJex2_rd1.0.txt" u (180*$1/pi):($4) w l lt 1 lw 3 lc rgb "black" notitle "ave.",\
     K1an_Co_j*sin(pi*x/180)**2+K2an_Co_j*sin(pi*x/180)**4+K3an_Co_j*sin(pi*x/180)**6 w l lw 3 lc rgb "red" title "analytical",\
     K1num_Co_j*sin(pi*x/180)**2+K2num_Co_j*sin(pi*x/180)**4+K3num_Co_j*sin(pi*x/180)**6 w l lw 3 lc rgb "blue" title "numerical"
unset multiplot

reset





kb=13.80649
VTi=8.54*8.54*4.78
VV=8.5205*8.5205*4.7693
VCo=8.4*8.4*4.8
KuFe=47.65712/2.0

set nokey
set xrange [0:180]
set yrange [-175:175]
set mxtics 2
set mytics 2
set xtics 0,45,180
set ytics -200,100,200
set size 0.65,1.0

set output "SmFe11TiVCo_i_eng_phi_theta90.eps"
set multiplot
set y2tics
set my2tics 2
set y2range [-148*2*kb/VTi:148*2*kb/VTi]
set y2tics offset -0.5,0

lmgn=0.12
rmgn=0.55
bmgn=0.12
tmgn=0.98
xsize=(rmgn-lmgn)
ysize=(tmgn-bmgn)/3.0
# set label position
lblx=5
lbly=135
set tmargin at screen bmgn+3*ysize
set bmargin at screen bmgn+2*ysize
set lmargin at screen lmgn
set rmargin at screen lmgn+xsize
set format x ""
set ytics nomirror
set label "(a) SmFe_{11}Ti (8{/Times-Italic i} site)" at lblx,lbly
set nokey
plot \
      0 w l lt 1 lw 1 lc rgb "black" notitle,\
     "SmFe11Ti_i_fullexp_lambda411_opencore_T0.0_eng_phi_theta90_nJex2_rd1.0.txt" u (180*$1/pi):($2+KuFe) w l lt 1 lw 3 lc rgb "red" notitle "Sm1",\
     "SmFe11Ti_i_fullexp_lambda411_opencore_T0.0_eng_phi_theta90_nJex2_rd1.0.txt" u (180*$1/pi):($3+KuFe) w l lt 1 lw 3 lc rgb "blue" notitle "Sm2",\
     "SmFe11Ti_i_fullexp_lambda411_opencore_T0.0_eng_phi_theta90_nJex2_rd1.0.txt" u (180*$1/pi):(($2+$3)/2.0+KuFe) w l lt 1 lw 3 lc rgb "dark-violet" notitle "ave.",\
     "SmFe11Ti_i_fullexp_lambda411_opencore_T0.0_eng_an_phi_theta90_rd1.0.txt" u (180*$1/pi):($2+KuFe) w l lt 2 lw 3 lc rgb "red" notitle "Sm1",\
     "SmFe11Ti_i_fullexp_lambda411_opencore_T0.0_eng_an_phi_theta90_rd1.0.txt" u (180*$1/pi):($3+KuFe) w l lt 2 lw 3 lc rgb "blue" notitle "Sm2",\
     "SmFe11Ti_i_fullexp_lambda411_opencore_T0.0_eng_an_phi_theta90_rd1.0.txt" u (180*$1/pi):(($2+$3)/2.0+KuFe) w l lt 2 lw 3 lc rgb "black" notitle "ave."

unset label
unset key
set ylabel "Magnetic Anisotropy Energy [K]" 0.5,0
set y2label "Magnetic Anisotropy Energy [MJ/m^3]" -1.5,0
set format y
set tmargin at screen bmgn+2*ysize
set bmargin at screen bmgn+ysize
set lmargin at screen lmgn
set rmargin at screen lmgn+xsize
set label "(b) SmFe_{11}V (8{/Times-Italic i} site)" at lblx,lbly
set y2range [-148*2*kb/VV:148*2*kb/VV]
plot \
     0 w l lt 1 lw 1 lc rgb "black" notitle,\
     "SmFe11V_i_fullexp_lambda411_opencore_T0.0_eng_phi_theta90_nJex2_rd1.0.txt" u (180*$1/pi):($2+KuFe) w l lt 1 lw 3 lc rgb "red" title "Sm1",\
     "SmFe11V_i_fullexp_lambda411_opencore_T0.0_eng_phi_theta90_nJex2_rd1.0.txt" u (180*$1/pi):($3+KuFe) w l lt 1 lw 3 lc rgb "blue" title "Sm2",\
     "SmFe11V_i_fullexp_lambda411_opencore_T0.0_eng_phi_theta90_nJex2_rd1.0.txt" u (180*$1/pi):(($2+$3)/2.0+KuFe) w l lt 1 lw 3 lc rgb "dark-violet" title "ave.",\
     "SmFe11V_i_fullexp_lambda411_opencore_T0.0_eng_an_phi_theta90_rd1.0.txt" u (180*$1/pi):($2+KuFe) w l lt 2 lw 3 lc rgb "red" title "Sm1",\
     "SmFe11V_i_fullexp_lambda411_opencore_T0.0_eng_an_phi_theta90_rd1.0.txt" u (180*$1/pi):($3+KuFe) w l lt 2 lw 3 lc rgb "blue" title "Sm2",\
     "SmFe11V_i_fullexp_lambda411_opencore_T0.0_eng_an_phi_theta90_rd1.0.txt" u (180*$1/pi):(($2+$3)/2.0+KuFe) w l lt 2 lw 3 lc rgb "black" title "ave."

unset label
unset ylabel
unset y2label
set format x
set format y
set xtics
set tmargin at screen bmgn+ysize
set bmargin at screen bmgn
set lmargin at screen lmgn
set rmargin at screen lmgn+xsize
set label "(c) SmFe_{11}Co (8{/Times-Italic i} site)" at lblx,-lbly+0.2
set xlabel "{/Symbol-Oblique f}^{TM} [degree]" 0.0,0
set key reverse Left
set key right bottom
set y2range [-148*2*kb/VCo:148*2*kb/VCo]
plot \
     0 w l lt 1 lw 1 lc rgb "black" notitle,\
     "SmFe11Co_i_fullexp_lambda411_opencore_T0.0_eng_phi_theta90_nJex2_rd1.0.txt" u (180*$1/pi):($2+KuFe) w l lt 1 lw 3 lc rgb "red" title "Sm1",\
     "SmFe11Co_i_fullexp_lambda411_opencore_T0.0_eng_phi_theta90_nJex2_rd1.0.txt" u (180*$1/pi):($3+KuFe) w l lt 1 lw 3 lc rgb "blue" title "Sm2",\
     "SmFe11Co_i_fullexp_lambda411_opencore_T0.0_eng_phi_theta90_nJex2_rd1.0.txt" u (180*$1/pi):(($2+$3)/2.0+KuFe) w l lt 1 lw 3 lc rgb "dark-violet" title "ave.",\
     "SmFe11Co_i_fullexp_lambda411_opencore_T0.0_eng_an_phi_theta90_rd1.0.txt" u (180*$1/pi):($2+KuFe) w l lt 2 lw 3 lc rgb "red" title "Sm1",\
     "SmFe11Co_i_fullexp_lambda411_opencore_T0.0_eng_an_phi_theta90_rd1.0.txt" u (180*$1/pi):($3+KuFe) w l lt 2 lw 3 lc rgb "blue" title "Sm2",\
     "SmFe11Co_i_fullexp_lambda411_opencore_T0.0_eng_an_phi_theta90_rd1.0.txt" u (180*$1/pi):(($2+$3)/2.0+KuFe) w l lt 2 lw 3 lc rgb "black" title "ave."

unset multiplot
reset




















kb=13.80649
VTi=8.54*8.54*4.78
VV=8.5205*8.5205*4.7693
KuFe=47.65712/2.0

set xrange [0:180]
set yrange [-175:175]
set mxtics 2
set mytics 2
set output "SmFe11TiV_i_eng_phi_theta90.eps"
ysize=0.8
set size 0.5,ysize
tmgn=0.95*ysize
bmgn=0.12*ysize
lmgn=0.08
rmgn=0.44
size=(tmgn-bmgn)/2.0
set multiplot
set lmargin at screen lmgn
set rmargin at screen rmgn
set tmargin at screen tmgn
set bmargin at screen bmgn+size
set format x ""
set xtics 0,45,180
set ytics nomirror
set y2tics
set my2tics 2
set y2range [-148*2*kb/VTi:148*2*kb/VTi]
set y2tics offset -0.5,0
lbz=140
set label "(a) SmFe_{11}Ti (8{/Times-Italic i} site)" at 5,lbz
plot \
      0 w l lt 1 lw 1 lc rgb "black" notitle,\
     "SmFe11Ti_i_fullexp_lambda411_opencore_T0.0_eng_phi_theta90_nJex2_rd1.0.txt" u (180*$1/pi):($2+KuFe) w l lt 1 lw 3 lc rgb "red" notitle "Sm1",\
     "SmFe11Ti_i_fullexp_lambda411_opencore_T0.0_eng_phi_theta90_nJex2_rd1.0.txt" u (180*$1/pi):($3+KuFe) w l lt 1 lw 3 lc rgb "blue" notitle "Sm2",\
     "SmFe11Ti_i_fullexp_lambda411_opencore_T0.0_eng_phi_theta90_nJex2_rd1.0.txt" u (180*$1/pi):(($2+$3)/2.0+KuFe) w l lt 1 lw 3 lc rgb "dark-violet" notitle "ave.",\
     "SmFe11Ti_i_fullexp_lambda411_opencore_T0.0_eng_an_phi_theta90_rd1.0.txt" u (180*$1/pi):($2+KuFe) w l lt 2 lw 1.5 lc rgb "red" notitle "Sm 1",\
     "SmFe11Ti_i_fullexp_lambda411_opencore_T0.0_eng_an_phi_theta90_rd1.0.txt" u (180*$1/pi):($3+KuFe) w l lt 2 lw 1.5 lc rgb "blue" notitle "Sm 2",\
     "SmFe11Ti_i_fullexp_lambda411_opencore_T0.0_eng_an_phi_theta90_rd1.0.txt" u (180*$1/pi):(($2+$3)/2.0+KuFe) w l lt 2 lw 1.5 lc rgb "dark-violet" notitle "m0 and m4"
unset key
set key reverse Left
set key 160,-75
unset label
set xlabel "{/Symbol-Oblique f}^{TM} [deg.]"
set ylabel "Magnetic Anisotropy Energy [K/atom]" 0.9,6
set y2label "Magnetic Anisotropy Energy [MJ/m^3]" -2.0,6
set lmargin at screen lmgn
set rmargin at screen rmgn
set tmargin at screen bmgn+size
set bmargin at screen bmgn
set format x
set y2range [-148*2*kb/VV:148*2*kb/VV]
set label "(b) SmFe_{11}V (8{/Times-Italic i} site)" at 5,lbz
plot \
     0 w l lt 1 lw 1 lc rgb "black" notitle,\
     "SmFe11V_i_fullexp_lambda411_opencore_T0.0_eng_phi_theta90_nJex2_rd1.0.txt" u (180*$1/pi):($2+KuFe) w l lt 1 lw 3 lc rgb "red" title "Sm1",\
     "SmFe11V_i_fullexp_lambda411_opencore_T0.0_eng_phi_theta90_nJex2_rd1.0.txt" u (180*$1/pi):($3+KuFe) w l lt 1 lw 3 lc rgb "blue" title "Sm2",\
     "SmFe11V_i_fullexp_lambda411_opencore_T0.0_eng_phi_theta90_nJex2_rd1.0.txt" u (180*$1/pi):(($2+$3)/2.0+KuFe) w l lt 1 lw 3 lc rgb "dark-violet" title "ave.",\
     "SmFe11V_i_fullexp_lambda411_opencore_T0.0_eng_an_phi_theta90_rd1.0.txt" u (180*$1/pi):($2+KuFe) w l lt 2 lw 1.5 lc rgb "red" notitle "Sm 1",\
     "SmFe11V_i_fullexp_lambda411_opencore_T0.0_eng_an_phi_theta90_rd1.0.txt" u (180*$1/pi):($3+KuFe) w l lt 2 lw 1.5 lc rgb "blue" notitle "Sm 2",\
     "SmFe11V_i_fullexp_lambda411_opencore_T0.0_eng_an_phi_theta90_rd1.0.txt" u (180*$1/pi):(($2+$3)/2.0+KuFe) w l lt 2 lw 1.5 lc rgb "dark-violet" notitle "total \n m0 and m4"
unset multiplot
reset







set size 0.5,0.5
set key reverse Left
set xlabel "phi [deg.]"
set ylabel "energy [K]"
set output "SmFe11Ti_i_eng_phi_theta90.eps"
plot \
     "SmFe11Ti_i_fullexp_lambda411_opencore_T0.0_eng_phi_theta90_nJex2_rd1.0.txt" u (180*$1/pi):2 w l lt 1 lc rgb "red" title "Sm 1",\
     "SmFe11Ti_i_fullexp_lambda411_opencore_T0.0_eng_phi_theta90_nJex2_rd1.0.txt" u (180*$1/pi):3 w l lt 1 lc rgb "blue" title "Sm 2",\
     "SmFe11Ti_i_fullexp_lambda411_opencore_T0.0_eng_phi_theta90_nJex2_rd1.0.txt" u (180*$1/pi):(($2+$3)/2.0) w l lt 1 lc rgb "dark-violet" title "full",\
     "SmFe11Ti_i_fullexp_lambda411_opencore_T0.0_eng_phi_theta90_nJex2_rd1.0.txt" u (180*$1/pi):2 w l lt 2 lc rgb "red" notitle "Sm 1",\
     "SmFe11Ti_i_fullexp_lambda411_opencore_T0.0_eng_phi_theta90_nJex2_rd1.0.txt" u (180*$1/pi):2 w l lt 2 lc rgb "blue" notitle "Sm 2",\
     "SmFe11Ti_i_fullexp_lambda411_opencore_T0.0_eng_phi_theta90_nJex2_rd1.0.txt" u (180*$1/pi):2 w l lt 2 lc rgb "dark-violet" title "m0 and m4"
reset


set size 0.5,0.5
set key reverse Left
set xlabel "phi [deg.]"
set ylabel "energy [K]"
set output "SmFe11V_i_eng_phi_theta90.eps"
plot \
     "SmFe11V_i_fullexp_lambda411_opencore_T0.0_eng_phi_theta90_nJex2_rd1.0.txt" u (180*$1/pi):2 w l lt 1 lc rgb "red" title "Sm 1",\
     "SmFe11V_i_fullexp_lambda411_opencore_T0.0_eng_phi_theta90_nJex2_rd1.0.txt" u (180*$1/pi):3 w l lt 1 lc rgb "blue" title "Sm 2",\
     "SmFe11V_i_fullexp_lambda411_opencore_T0.0_eng_phi_theta90_nJex2_rd1.0.txt" u (180*$1/pi):(($2+$3)/2.0) w l lt 1 lc rgb "dark-violet" title "full",\
     "SmFe11V_i_fullexp_lambda411_opencore_T0.0_eng_phi_theta90_nJex2_rd1.0.txt" u (180*$1/pi):2 w l lt 2 lc rgb "red" notitle "Sm 1",\
     "SmFe11V_i_fullexp_lambda411_opencore_T0.0_eng_phi_theta90_nJex2_rd1.0.txt" u (180*$1/pi):2 w l lt 2 lc rgb "blue" notitle "Sm 2",\
     "SmFe11V_i_fullexp_lambda411_opencore_T0.0_eng_phi_theta90_nJex2_rd1.0.txt" u (180*$1/pi):2 w l lt 2 lc rgb "dark-violet" title "m0 and m4"
reset




set size 0.5,0.5
set key reverse Left
set xlabel "phi"
set ylabel "energy"
set output "SmFe11Ti_i_eng_phi_theta0.eps"
plot \
     "SmFe11Ti_i_fullexp_lambda411_opencore_T0.0_eng_phi_theta0_nJex2_rd1.0.txt" u (180*$1/pi):2 w l lt 1 lc rgb "red" title "Sm 1",\
     "SmFe11Ti_i_fullexp_lambda411_opencore_T0.0_eng_phi_theta0_nJex2_rd1.0.txt" u (180*$1/pi):3 w l lt 1 lc rgb "blue" title "Sm 2",\
     "SmFe11Ti_i_fullexp_lambda411_opencore_T0.0_eng_phi_theta0_nJex2_rd1.0.txt" u (180*$1/pi):(($2+$3)/2.0) w l lt 1 lc rgb "dark-violet" title "full",\
     "SmFe11Ti_i_fullexp_lambda411_opencore_T0.0_eng_phi_theta0_nJex2_rd1.0.txt" u (180*$1/pi):2 w l lt 2 lc rgb "red" notitle "Sm 1",\
     "SmFe11Ti_i_fullexp_lambda411_opencore_T0.0_eng_phi_theta0_nJex2_rd1.0.txt" u (180*$1/pi):3 w l lt 2 lc rgb "blue" notitle "Sm 2",\
     "SmFe11Ti_i_fullexp_lambda411_opencore_T0.0_eng_phi_theta0_nJex2_rd1.0.txt" u (180*$1/pi):2 w l lt 2 lc rgb "dark-violet" title "m0 and m4"
reset

set size 0.5,0.5
set key reverse Left
set xlabel "{/Symbol f} [deg.]"
set ylabel "{/Symbol Df} [deg.]"
set yrange [-80:80]
set output "SmFe11Ti_i_delta_phi_SL.eps"
set yrange [-10:10]
plot \
     "SmFe11Ti_i_fullexp_lambda411_opencore_T0.0_SLangle_phi_theta0_nJex2_rd1.0.txt" u (180*$1/pi):(180*($2-$1-pi)/pi) w l lt 2 lc rgb "red" title "Sm 1 S",\
     "SmFe11Ti_i_fullexp_lambda411_opencore_T0.0_SLangle_phi_theta0_nJex2_rd1.0.txt" u (180*$1/pi):(180*($3-$1-pi)/pi) w l lt 2 lc rgb "blue" title "Sm 2 S",\
     "SmFe11Ti_i_fullexp_lambda411_opencore_T0.0_SLangle_phi_theta0_nJex2_rd1.0.txt" u (180*$1/pi):(180*($4-$1)/pi) w l lt 1 lc rgb "red" title "Sm 1 L",\
     "SmFe11Ti_i_fullexp_lambda411_opencore_T0.0_SLangle_phi_theta0_nJex2_rd1.0.txt" u (180*$1/pi):(180*($5-$1)/pi) w l lt 1 lc rgb "blue" title "Sm 2 L"

set size 0.5,0.5
set key reverse Left
set xlabel "{/Symbol f} [deg.]"
set ylabel "{/Symbol Df} [deg.]"
set output "SmFe11Ti_i_delta_phi_SL_theta90.eps"
set yrange [-10:10]
plot \
     "SmFe11Ti_i_fullexp_lambda411_opencore_T0.0_SLangle_phi_theta90_nJex2_rd1.0.txt" u (180*$1/pi):(180*($2+$1-pi)/pi) w l lt 2 lc rgb "red" title "Sm 1 S",\
     "SmFe11Ti_i_fullexp_lambda411_opencore_T0.0_SLangle_phi_theta90_nJex2_rd1.0.txt" u (180*$1/pi):(180*($3+$1-pi)/pi) w l lt 2 lc rgb "blue" title "Sm 2 S",\
     "SmFe11Ti_i_fullexp_lambda411_opencore_T0.0_SLangle_phi_theta90_nJex2_rd1.0.txt" u (180*$1/pi):(180*($4-$1)/pi) w l lt 1 lc rgb "red" title "Sm 1 L",\
     "SmFe11Ti_i_fullexp_lambda411_opencore_T0.0_SLangle_phi_theta90_nJex2_rd1.0.txt" u (180*$1/pi):(180*($5-$1)/pi) w l lt 1 lc rgb "blue" title "Sm 2 L"


reset














set size 0.5,0.5
set key reverse Left
set xlabel "phi"
set ylabel "energy"
set output "SmFe11V_i_eng_phi_theta0.eps"
plot \
     "SmFe11V_i_fullexp_lambda411_opencore_T0.0_eng_phi_theta0_nJex2_rd1.0.txt" u (180*$1/pi):2 w l lt 1 lc rgb "red" title "Sm 1",\
     "SmFe11V_i_fullexp_lambda411_opencore_T0.0_eng_phi_theta0_nJex2_rd1.0.txt" u (180*$1/pi):3 w l lt 1 lc rgb "blue" title "Sm 2",\
     "SmFe11V_i_fullexp_lambda411_opencore_T0.0_eng_phi_theta0_nJex2_rd1.0.txt" u (180*$1/pi):(($2+$3)/2.0) w l lt 1 lc rgb "dark-violet" title "full",\
     "SmFe11V_i_fullexp_lambda411_opencore_T0.0_eng_phi_theta0_nJex2_rd1.0.txt" u (180*$1/pi):2 w l lt 2 lc rgb "red" notitle "Sm 1",\
     "SmFe11V_i_fullexp_lambda411_opencore_T0.0_eng_phi_theta0_nJex2_rd1.0.txt" u (180*$1/pi):3 w l lt 2 lc rgb "blue" notitle "Sm 2",\
     "SmFe11V_i_fullexp_lambda411_opencore_T0.0_eng_phi_theta0_nJex2_rd1.0.txt" u (180*$1/pi):2 w l lt 2 lc rgb "dark-violet" title "m0 and m4"
reset

set size 0.5,0.5
set key reverse Left
set xlabel "{/Symbol f} [deg.]"
set ylabel "{/Symbol Df} [deg.]"
set yrange [-80:80]
set output "SmFe11V_i_delta_phi_SL.eps"
set yrange [-10:10]
plot \
     "SmFe11V_i_fullexp_lambda411_opencore_T0.0_SLangle_phi_theta0_nJex2_rd1.0.txt" u (180*$1/pi):(180*($2-$1-pi)/pi) w l lt 2 lc rgb "red" title "Sm 1 S",\
     "SmFe11V_i_fullexp_lambda411_opencore_T0.0_SLangle_phi_theta0_nJex2_rd1.0.txt" u (180*$1/pi):(180*($3-$1-pi)/pi) w l lt 2 lc rgb "blue" title "Sm 2 S",\
     "SmFe11V_i_fullexp_lambda411_opencore_T0.0_SLangle_phi_theta0_nJex2_rd1.0.txt" u (180*$1/pi):(180*($4-$1)/pi) w l lt 1 lc rgb "red" title "Sm 1 L",\
     "SmFe11V_i_fullexp_lambda411_opencore_T0.0_SLangle_phi_theta0_nJex2_rd1.0.txt" u (180*$1/pi):(180*($5-$1)/pi) w l lt 1 lc rgb "blue" title "Sm 2 L"

set size 0.5,0.5
set key reverse Left
set xlabel "{/Symbol f} [deg.]"
set ylabel "{/Symbol Df} [deg.]"
set yrange [-80:80]
set output "SmFe11V_i_delta_phi_SL_theta90.eps"
set yrange [-10:10]
plot \
     "SmFe11V_i_fullexp_lambda411_opencore_T0.0_SLangle_phi_theta90_nJex2_rd1.0.txt" u (180*$1/pi):(180*($2-$1-pi)/pi) w l lt 2 lc rgb "red" title "Sm 1 S",\
     "SmFe11V_i_fullexp_lambda411_opencore_T0.0_SLangle_phi_theta90_nJex2_rd1.0.txt" u (180*$1/pi):(180*($3-$1-pi)/pi) w l lt 2 lc rgb "blue" title "Sm 2 S",\
     "SmFe11V_i_fullexp_lambda411_opencore_T0.0_SLangle_phi_theta90_nJex2_rd1.0.txt" u (180*$1/pi):(180*($4-$1)/pi) w l lt 1 lc rgb "red" title "Sm 1 L",\
     "SmFe11V_i_fullexp_lambda411_opencore_T0.0_SLangle_phi_theta90_nJex2_rd1.0.txt" u (180*$1/pi):(180*($5-$1)/pi) w l lt 1 lc rgb "blue" title "Sm 2 L"


reset





set terminal postscript eps enhanced color font "Times-Roman,13"

set size 0.5,0.6
set output "Helmholtz_eng_T0_T400_coexist.eps"
set multiplot
set key reverse Left
# set left and right margin
lmgn=0.08
rmgn=0.48
bmgn=0.085
tmgn=0.84
xsize=(rmgn-lmgn)
ysize=(tmgn-bmgn)/3.0
# set label position
lblx=250
lbly=16.5
set tmargin at screen bmgn+2*ysize
set bmargin at screen bmgn+ysize
set lmargin at screen lmgn
set rmargin at screen lmgn+xsize
set label "(a) {/Times-Italic T}=0 K" at 0.05,7
set format x ""
set ylabel "{/Times-Italic F}({/Times-BoldItalic e}^{/Times-Italic=8 TM}, {/Times-Italic T}) [MJ/m^3]" 0.8,-5
set xrange [0:1]
set yrange [-1:8]
set ytics -2,2,8
set mxtics 2
set mytics 2
K1TM0=1.966
set key 0.6,6
plot \
     0 notitle lt 1 lw 0.5 lc rgb "black",\
     "SmFe12_lambda411_opencore_T0.0_eng_theta_phi0_nJex2_rd1.0.txt"\
     u (sin($1)):4 w l lw 3 lt 2 lc rgb "blue" title"{/Times-Roman=12 statistics}",\
     "SmFe12_lambda411_opencore_T0.0_eng_an_theta_phi0_rd1.0.txt" \
     u (sin($1)):4 w l lw 3 lt 1 lc rgb "blue" title "{/Times-Roman=12 model C}",\
     K1TM0*x**2 lw 2 lt 5 lc rgb "black" title "{/Times-Roman=12 Fe sublattice}"
unset ylabel
unset label
unset arrow
set format x
set tmargin at screen bmgn+ysize
set bmargin at screen bmgn
set lmargin at screen lmgn
set rmargin at screen lmgn+xsize
set label "(b) {/Times-Italic T}=400 K" at 0.05,7
set xlabel "{/Times-BoldItalic e^{/Times-Italic=8 TM}}{/Symbol \327} {/Times-BoldItalic n}_{/Times-Italic a}" 0,0.4
K1TM400=0.387
plot \
     0 notitle lt 1 lw 0.5 lc rgb "black",\
     "SmFe12_lambda411_opencore_T400.0_eng_theta_phi0_nJex2_rd1.0.txt"\
     u (sin($1)):4 w l lw 3 lt 2 lc rgb "red" title"{/Times-Roman=12 statistics}",\
     "SmFe12_lambda411_opencore_T400.0_eng_an_theta_phi0_rd1.0.txt" \
     u (sin($1)):4 w l lw 3 lt 1 lc rgb "red" title "{/Times-Roman=12 model C}",\
     K1TM400*x**2 lw 2 lt 5 lc rgb "black" title "{/Times-Roman=12 Fe sublattice}"

unset multiplot
reset

