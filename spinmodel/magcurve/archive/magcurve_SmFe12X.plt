set terminal postscript eps enhanced color font "Times-Roman,12.5"


set xrange [0:29]
set yrange [0:2.5]
set nokey
set ytics 0,1,3
set mxtics 2
set mytics 2
keyx=28
keyy=1.4
set size 0.5,0.55
set output "magcurve_SmFe12X_annum_100_001.eps"
set multiplot
# set left and right margin
lmgn=0.05
rmgn=0.49
bmgn=0.08
tmgn=0.54
xsize=(rmgn-lmgn)/2.0
ysize=(tmgn-bmgn)/2.0
# set label position

lblx=19
lbly=2.2
set tmargin at screen bmgn+2*ysize
set bmargin at screen bmgn+ysize
set lmargin at screen lmgn
set rmargin at screen lmgn+xsize-0.002
set format x ""
set label "(a)" at 1.0,2.3
set label "{/Times-Italic T}=0 K" at lblx,lbly
set xtics scale 3.0/4.0,3.0/8.0
set ytics scale 3.0/4.0,3.0/8.0
unset key
plot \
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T0.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 1 lw 3 lc rgb "black" notitle "SmFe_{12}",\
     "SmFe12H_lambda411_opencore_wK1Fe_magcurve_T0.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 1 lw 3 lc rgb "red" notitle "SmFeH_{12}",\
     "SmFe12B_lambda411_opencore_wK1Fe_magcurve_T0.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 1 lw 3 lc rgb "dark-green" notitle "0 K",\
     "SmFe12C_lambda411_opencore_wK1Fe_magcurve_T0.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 1 lw 3 lc rgb "blue" notitle "0 K",\
     "SmFe12N_lambda411_opencore_wK1Fe_magcurve_T0.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 1 lw 3 lc rgb "dark-violet" notitle "0 K",\
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 2 lw 1.5 lc rgb "black" title "0 K",\
     "SmFe12H_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 2 lw 1.5 lc rgb "red" title "0 K",\
     "SmFe12B_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 2 lw 1.5 lc rgb "dark-green" title "0 K",\
     "SmFe12C_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 2 lw 1.5 lc rgb "blue" title "0 K",\
     "SmFe12N_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 2 lw 1.5 lc rgb "dark-violet" title "0 K"
unset label
set format y ""
set tmargin at screen bmgn+2*ysize
set bmargin at screen bmgn+ysize
set lmargin at screen lmgn+xsize+0.002
set rmargin at screen rmgn
set format x ""
set label "{/Times-Italic T}=0 K" at lblx,lbly
set label "(b)" at 1.0,2.3
unset key 
plot \
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T0.0_001_mixexcf2_rd1.0.txt" u 1:2 w l lt 1 lw 3 lc rgb "black" title "SmFe_{12}",\
     "SmFe12H_lambda411_opencore_wK1Fe_magcurve_T0.0_001_mixexcf2_rd1.0.txt" u 1:2 w l lt 1 lw 3 lc rgb "red" title "SmFeH_{12}",\
     "SmFe12B_lambda411_opencore_wK1Fe_magcurve_T0.0_001_mixexcf2_rd1.0.txt" u 1:2 w l lt 1 lw 3 lc rgb "dark-green" notitle "0 K",\
     "SmFe12C_lambda411_opencore_wK1Fe_magcurve_T0.0_001_mixexcf2_rd1.0.txt" u 1:2 w l lt 1 lw 3 lc rgb "blue" notitle "0 K",\
     "SmFe12N_lambda411_opencore_wK1Fe_magcurve_T0.0_001_mixexcf2_rd1.0.txt" u 1:2 w l lt 1 lw 3 lc rgb "dark-violet" notitle "0 K",\
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T0.0_001_nJex2_rd1.0.txt" u 1:3 w l lt 2 lw 1.5 lc rgb "black" title "0 K",\
     "SmFe12H_lambda411_opencore_wK1Fe_magcurve_T0.0_001_nJex2_rd1.0.txt" u 1:3 w l lt 2 lw 1.5 lc rgb "red" title "0 K",\
     "SmFe12B_lambda411_opencore_wK1Fe_magcurve_T0.0_001_nJex2_rd1.0.txt" u 1:3 w l lt 2 lw 1.5 lc rgb "dark-green" title "0 K",\
     "SmFe12C_lambda411_opencore_wK1Fe_magcurve_T0.0_001_nJex2_rd1.0.txt" u 1:3 w l lt 2 lw 1.5 lc rgb "blue" title "0 K",\
     "SmFe12N_lambda411_opencore_wK1Fe_magcurve_T0.0_001_nJex2_rd1.0.txt" u 1:3 w l lt 2 lw 1.5 lc rgb "dark-violet" title "0 K"

set tmargin at screen bmgn+ysize
set bmargin at screen bmgn
set lmargin at screen lmgn
set rmargin at screen lmgn+xsize-0.002
set format x
set format y
unset label
set xlabel "{/Times-Italic B} [T] (|| {/Times-Italic a}-axis)" 0,0.2
set ylabel "{/Symbol-Oblique m}_0{/Times-BoldItalic M}_{s}({/Times-Italic T}){/Symbol \327}({/Times-BoldItalic B}/{/Times-Italic B}) [T]" 1.0,4.5
#set ylabel "Magnetization [T]" 1.0,5.5
#set ylabel "{/Symbol-Oblique m}_0{/Bold-Times-Italic M}_s({/Times-Italic T}){/Bold-Times-Italic n}_{a,c} [T]" 1.0,5.5
set label "{/Times-Italic T}/{/Times-Italic T}_C=0.721" at lblx-2,lbly
set key reverse Left
set key at 26,1.8
plot \
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T400.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 1 lw 3 lc rgb "black" title "SmFe_{12} [24]",\
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 2 lw 1.5 lc rgb "black" notitle "SmFe_{12} [99]",\
     "SmFe12H_lambda411_opencore_wK1Fe_magcurve_T400.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 1 lw 3 lc rgb "red" notitle "SmFeH_{12}",\
     "SmFe12B_lambda411_opencore_wK1Fe_magcurve_T400.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 1 lw 3 lc rgb "dark-green" notitle "0 K",\
     "SmFe12C_lambda411_opencore_wK1Fe_magcurve_T400.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 1 lw 3 lc rgb "blue" notitle "0 K",\
     "SmFe12N_lambda411_opencore_wK1Fe_magcurve_T400.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 1 lw 3 lc rgb "dark-violet" notitle "0 K",\
     "SmFe12H_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 2 lw 1.5 lc rgb "red" notitle "SmFeH_{12}",\
     "SmFe12B_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 2 lw 1.5 lc rgb "dark-green" notitle "SmFe_{12}B",\
     "SmFe12C_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 2 lw 1.5 lc rgb "blue" notitle "SmFe_{12}C",\
     "SmFe12N_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 2 lw 1.5 lc rgb "dark-violet" notitle "SmFe_{12}N"
unset ylabel

set tmargin at screen bmgn+ysize
set bmargin at screen bmgn
set lmargin at screen lmgn+xsize+0.002
set rmargin at screen rmgn
unset label
set label "{/Times-Italic T}/{/Times-Italic T}_C=0.721" at lblx-2,lbly
set xlabel "{/Times-Italic B} [T] (|| {/Times-Italic c}-axis)" 0,0.2
set format y ""
set key reverse Left
set key at -3.2,1.2
plot \
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T400.0_001_mixexcf2_rd1.0.txt" u 1:2 w l lt 1 lw 3 lc rgb "black" notitle "SmFe_{12}",\
     "SmFe12H_lambda411_opencore_wK1Fe_magcurve_T400.0_001_mixexcf2_rd1.0.txt" u 1:2 w l lt 1 lw 3 lc rgb "red" title "SmFe_{12}H",\
     "SmFe12B_lambda411_opencore_wK1Fe_magcurve_T400.0_001_mixexcf2_rd1.0.txt" u 1:2 w l lt 1 lw 3 lc rgb "dark-green" title "SmFe_{12}B",\
     "SmFe12C_lambda411_opencore_wK1Fe_magcurve_T400.0_001_mixexcf2_rd1.0.txt" u 1:2 w l lt 1 lw 3 lc rgb "blue" title "SmFe_{12}C",\
     "SmFe12N_lambda411_opencore_wK1Fe_magcurve_T400.0_001_mixexcf2_rd1.0.txt" u 1:2 w l lt 1 lw 3 lc rgb "dark-violet" title "SmFe_{12}N",\
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T400.0_001_nJex2_rd1.0.txt" u 1:3 w l lt 2 lw 1.5 lc rgb "black" notitle "SmFe_{12}",\
     "SmFe12H_lambda411_opencore_wK1Fe_magcurve_T400.0_001_nJex2_rd1.0.txt" u 1:3 w l lt 2 lw 1.5 lc rgb "red" notitle "SmFe_{12}H",\
     "SmFe12B_lambda411_opencore_wK1Fe_magcurve_T400.0_001_nJex2_rd1.0.txt" u 1:3 w l lt 2 lw 1.5 lc rgb "dark-green" notitle "SmFe_{12}B",\
     "SmFe12C_lambda411_opencore_wK1Fe_magcurve_T400.0_001_nJex2_rd1.0.txt" u 1:3 w l lt 2 lw 1.5 lc rgb "blue" notitle "SmFe_{12}C",\
     "SmFe12N_lambda411_opencore_wK1Fe_magcurve_T400.0_001_nJex2_rd1.0.txt" u 1:3 w l lt 2 lw 1.5 lc rgb "dark-violet" notitle "SmFe_{12}N"
unset multiplot
reset




set xrange [0:29]
set yrange [0:2.5]
set nokey
set ytics 0,1,3
set mxtics 2
set mytics 2
keyx=28
keyy=1.4
set size 0.5,0.55
set output "magcurve_SmFe12X_annum2_100_001.eps"
set multiplot
# set left and right margin
lmgn=0.05
rmgn=0.49
bmgn=0.08
tmgn=0.54
xsize=(rmgn-lmgn)/2.0
ysize=(tmgn-bmgn)/2.0
# set label position

lblx=19
lbly=2.2
set tmargin at screen bmgn+2*ysize
set bmargin at screen bmgn+ysize
set lmargin at screen lmgn
set rmargin at screen lmgn+xsize-0.002
set format x ""
set label "(a)" at 1.0,2.3
set label "{/Times-Italic T}=0 K" at lblx,lbly
set xtics scale 3.0/4.0,3.0/8.0
set ytics scale 3.0/4.0,3.0/8.0
unset key
plot \
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T0.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 2 lw 1.5 lc rgb "black" notitle "SmFe_{12}",\
     "SmFe12H_lambda411_opencore_wK1Fe_magcurve_T0.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 2 lw 1.5 lc rgb "red" notitle "SmFeH_{12}",\
     "SmFe12B_lambda411_opencore_wK1Fe_magcurve_T0.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 2 lw 1.5 lc rgb "dark-green" notitle "0 K",\
     "SmFe12C_lambda411_opencore_wK1Fe_magcurve_T0.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 2 lw 1.5 lc rgb "blue" notitle "0 K",\
     "SmFe12N_lambda411_opencore_wK1Fe_magcurve_T0.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 2 lw 1.5 lc rgb "dark-violet" notitle "0 K",\
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "black" title "0 K",\
     "SmFe12H_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "red" title "0 K",\
     "SmFe12B_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "dark-green" title "0 K",\
     "SmFe12C_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "blue" title "0 K",\
     "SmFe12N_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "dark-violet" title "0 K"
unset label
set format y ""
set tmargin at screen bmgn+2*ysize
set bmargin at screen bmgn+ysize
set lmargin at screen lmgn+xsize+0.002
set rmargin at screen rmgn
set format x ""
set label "{/Times-Italic T}=0 K" at lblx,lbly
set label "(b)" at 1.0,2.3
unset key 
plot \
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T0.0_001_mixexcf2_rd1.0.txt" u 1:2 w l lt 2 lw 1.5 lc rgb "black" title "SmFe_{12}",\
     "SmFe12H_lambda411_opencore_wK1Fe_magcurve_T0.0_001_mixexcf2_rd1.0.txt" u 1:2 w l lt 2 lw 1.5 lc rgb "red" title "SmFeH_{12}",\
     "SmFe12B_lambda411_opencore_wK1Fe_magcurve_T0.0_001_mixexcf2_rd1.0.txt" u 1:2 w l lt 2 lw 1.5 lc rgb "dark-green" notitle "0 K",\
     "SmFe12C_lambda411_opencore_wK1Fe_magcurve_T0.0_001_mixexcf2_rd1.0.txt" u 1:2 w l lt 2 lw 1.5 lc rgb "blue" notitle "0 K",\
     "SmFe12N_lambda411_opencore_wK1Fe_magcurve_T0.0_001_mixexcf2_rd1.0.txt" u 1:2 w l lt 2 lw 1.5 lc rgb "dark-violet" notitle "0 K",\
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T0.0_001_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "black" title "0 K",\
     "SmFe12H_lambda411_opencore_wK1Fe_magcurve_T0.0_001_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "red" title "0 K",\
     "SmFe12B_lambda411_opencore_wK1Fe_magcurve_T0.0_001_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "dark-green" title "0 K",\
     "SmFe12C_lambda411_opencore_wK1Fe_magcurve_T0.0_001_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "blue" title "0 K",\
     "SmFe12N_lambda411_opencore_wK1Fe_magcurve_T0.0_001_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "dark-violet" title "0 K"

set tmargin at screen bmgn+ysize
set bmargin at screen bmgn
set lmargin at screen lmgn
set rmargin at screen lmgn+xsize-0.002
set format x
set format y
unset label
set xlabel "{/Times-Italic B} [T] (|| {/Times-Italic a}-axis)" 0,0.2
set ylabel "{/Symbol-Oblique m}_0{/Times-BoldItalic M}_{s}({/Times-Italic T}){/Symbol \327}({/Times-BoldItalic B}/{/Times-Italic B}) [T]" 1.0,4.5
#set ylabel "Magnetization [T]" 1.0,5.5
#set ylabel "{/Symbol-Oblique m}_0{/Bold-Times-Italic M}_s({/Times-Italic T}){/Bold-Times-Italic n}_{a,c} [T]" 1.0,5.5
set label "{/Times-Italic T}/{/Times-Italic T}_C=0.721" at lblx-2,lbly
set key reverse Left
set key at 26,1.8
plot \
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T400.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 2 lw 1.5 lc rgb "black" notitle "SmFe_{12}",\
     "SmFe12H_lambda411_opencore_wK1Fe_magcurve_T400.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 2 lw 1.5 lc rgb "red" notitle "SmFeH_{12}",\
     "SmFe12B_lambda411_opencore_wK1Fe_magcurve_T400.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 2 lw 1.5 lc rgb "dark-green" notitle "0 K",\
     "SmFe12C_lambda411_opencore_wK1Fe_magcurve_T400.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 2 lw 1.5 lc rgb "blue" notitle "0 K",\
     "SmFe12N_lambda411_opencore_wK1Fe_magcurve_T400.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 2 lw 1.5 lc rgb "dark-violet" notitle "0 K",\
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "black" title "SmFe_{12} [99]",\
     "SmFe12H_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "red" notitle "SmFeH_{12}",\
     "SmFe12B_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "dark-green" notitle "SmFe_{12}B",\
     "SmFe12C_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "blue" notitle "SmFe_{12}C",\
     "SmFe12N_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "dark-violet" notitle "SmFe_{12}N"
unset ylabel

set tmargin at screen bmgn+ysize
set bmargin at screen bmgn
set lmargin at screen lmgn+xsize+0.002
set rmargin at screen rmgn
unset label
set label "{/Times-Italic T}/{/Times-Italic T}_C=0.721" at lblx-2,lbly
set xlabel "{/Times-Italic B} [T] (|| {/Times-Italic c}-axis)" 0,0.2
set format y ""
set key reverse Left
set key at -3.2,1.2
plot \
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T400.0_001_mixexcf2_rd1.0.txt" u 1:2 w l lt 2 lw 1.5 lc rgb "black" notitle "SmFe_{12}",\
     "SmFe12H_lambda411_opencore_wK1Fe_magcurve_T400.0_001_mixexcf2_rd1.0.txt" u 1:2 w l lt 2 lw 1.5 lc rgb "red" notitle "SmFe_{12}H",\
     "SmFe12B_lambda411_opencore_wK1Fe_magcurve_T400.0_001_mixexcf2_rd1.0.txt" u 1:2 w l lt 2 lw 1.5 lc rgb "dark-green" notitle "SmFe_{12}B",\
     "SmFe12C_lambda411_opencore_wK1Fe_magcurve_T400.0_001_mixexcf2_rd1.0.txt" u 1:2 w l lt 2 lw 1.5 lc rgb "blue" notitle "SmFe_{12}C",\
     "SmFe12N_lambda411_opencore_wK1Fe_magcurve_T400.0_001_mixexcf2_rd1.0.txt" u 1:2 w l lt 2 lw 1.5 lc rgb "dark-violet" notitle "SmFe_{12}N",\
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T400.0_001_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "black" notitle "SmFe_{12}",\
     "SmFe12H_lambda411_opencore_wK1Fe_magcurve_T400.0_001_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "red" title "SmFe_{12}H",\
     "SmFe12B_lambda411_opencore_wK1Fe_magcurve_T400.0_001_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "dark-green" title "SmFe_{12}B",\
     "SmFe12C_lambda411_opencore_wK1Fe_magcurve_T400.0_001_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "blue" title "SmFe_{12}C",\
     "SmFe12N_lambda411_opencore_wK1Fe_magcurve_T400.0_001_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "dark-violet" title "SmFe_{12}N"
unset multiplot
reset





set xrange [0:29]
set yrange [0:2.5]
set nokey
set ytics 0,1,3
set mxtics 2
set mytics 2
keyx=28
keyy=1.4
set size 0.5,0.65
set output "magcurve_SmFe12X_annum.eps"
set multiplot
# set left and right margin
lmgn=0.06
#rmgn=0.4
rmgn=0.48
bmgn=0.08
tmgn=0.6
xsize=(rmgn-lmgn)
ysize=(tmgn-bmgn)/2.0
# set label position
lblx=0.5
lbly=2.2
set tmargin at screen bmgn+2*ysize
set bmargin at screen bmgn+ysize
set lmargin at screen lmgn
set rmargin at screen lmgn+xsize
set format x ""
set label "(a) SmFe_{12}" at lblx,lbly
set key 16,1.1
set key reverse Left
unset key
set xtics scale 0,0
plot \
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T0.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 2 lw 2 lc rgb "black" title "0 K",\
     "SmFe12H_lambda411_opencore_wK1Fe_magcurve_T0.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 2 lw 2 lc rgb "red" title "0 K",\
     "SmFe12B_lambda411_opencore_wK1Fe_magcurve_T0.0_001_mixexcf2_rd1.0.txt" u 1:2 w l lt 2 lw 2 lc rgb "dark-green" title "0 K",\
     "SmFe12C_lambda411_opencore_wK1Fe_magcurve_T0.0_001_mixexcf2_rd1.0.txt" u 1:2 w l lt 2 lw 2 lc rgb "blue" title "0 K",\
     "SmFe12N_lambda411_opencore_wK1Fe_magcurve_T0.0_001_mixexcf2_rd1.0.txt" u 1:2 w l lt 2 lw 2 lc rgb "dark-violet" title "0 K",\
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 2 lc rgb "black" title "0 K",\
     "SmFe12H_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 2 lc rgb "red" title "0 K",\
     "SmFe12B_lambda411_opencore_wK1Fe_magcurve_T0.0_001_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 2 lc rgb "dark-green" title "0 K",\
     "SmFe12C_lambda411_opencore_wK1Fe_magcurve_T0.0_001_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 2 lc rgb "blue" title "0 K",\
     "SmFe12N_lambda411_opencore_wK1Fe_magcurve_T0.0_001_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 2 lc rgb "dark-violet" title "0 K"
set format y
set format x
set tmargin at screen bmgn+ysize
set bmargin at screen bmgn
set lmargin at screen lmgn
set rmargin at screen lmgn+xsize
unset label
set label "(b) SmFe_{12}H" at lblx,lbly
set xlabel "{/Times-Italic B} [T]" 0,0.2
unset key
set label "400 K" tc rgb "red" at graph 0.37,0.5 font ",12"
set label "300 K" tc rgb "orange" at graph 0.51,0.56 font ",12"
set label "200 K" tc rgb "dark-yellow" at graph 0.64,0.65 font ",12"
set label "100 K" tc rgb "dark-green" at graph 0.75,0.74 font ",12"
set label "0 K" tc rgb "blue" at graph 0.88,0.77 font ",12"
plot \
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T400.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 2 lw 2 lc rgb "black" title "400 K",\
     "SmFe12H_lambda411_opencore_wK1Fe_magcurve_T400.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 2 lw 2 lc rgb "red" title "400 K",\
     "SmFe12B_lambda411_opencore_wK1Fe_magcurve_T400.0_001_mixexcf2_rd1.0.txt" u 1:2 w l lt 2 lw 2 lc rgb "dark-green" title "400 K",\
     "SmFe12C_lambda411_opencore_wK1Fe_magcurve_T400.0_001_mixexcf2_rd1.0.txt" u 1:2 w l lt 2 lw 2 lc rgb "dark-blue" title "400 K",\
     "SmFe12N_lambda411_opencore_wK1Fe_magcurve_T400.0_001_mixexcf2_rd1.0.txt" u 1:2 w l lt 2 lw 2 lc rgb "dark-violet" title "400 K",\
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 2 lc rgb "black" title "400 K",\
     "SmFe12H_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 2 lc rgb "red" title "400 K",\
     "SmFe12B_lambda411_opencore_wK1Fe_magcurve_T400.0_001_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 2 lc rgb "dark-green" title "0 K",\
     "SmFe12C_lambda411_opencore_wK1Fe_magcurve_T400.0_001_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 2 lc rgb "blue" title "0 K",\
     "SmFe12N_lambda411_opencore_wK1Fe_magcurve_T400.0_001_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 2 lc rgb "dark-violet" title "0 K"
unset multiplot
reset












set terminal postscript eps enhanced color font "Times-Roman,12.5"

set xrange [0:29]
set yrange [0:2.5]
set nokey
set ytics 0,1,3
set mxtics 2
set mytics 2
keyx=28
keyy=1.4
set size 0.5,0.65
set output "magcurve_SmFe12X_num_100_001.eps"
set multiplot
# set left and right margin
lmgn=0.06
#rmgn=0.4
rmgn=0.48
bmgn=0.08
tmgn=0.6
xsize=(rmgn-lmgn)/2.0
ysize=(tmgn-bmgn)/2.0
# set label position
lblx=0.6
lbly=2.2

set tmargin at screen bmgn+2*ysize
set bmargin at screen bmgn+ysize
set lmargin at screen lmgn
set rmargin at screen lmgn+xsize
set format x ""
set label "(a) {/Times-Italic T}=0 K" at lblx,lbly
unset key
plot \
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "black" title "0 K",\
     "SmFe12H_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "red" title "0 K",\
     "SmFe12B_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "dark-green" title "0 K",\
     "SmFe12C_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "blue" title "0 K",\
     "SmFe12N_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "dark-violet" title "0 K"
unset label
set format y ""
set tmargin at screen bmgn+2*ysize
set bmargin at screen bmgn+ysize
set lmargin at screen lmgn+xsize
set rmargin at screen rmgn
set format x ""
set label "(b) {/Times-Italic T}=0 K" at lblx,lbly
unset key 
plot \
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T0.0_001_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "black" title "0 K",\
     "SmFe12H_lambda411_opencore_wK1Fe_magcurve_T0.0_001_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "red" title "0 K",\
     "SmFe12B_lambda411_opencore_wK1Fe_magcurve_T0.0_001_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "dark-green" title "0 K",\
     "SmFe12C_lambda411_opencore_wK1Fe_magcurve_T0.0_001_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "blue" title "0 K",\
     "SmFe12N_lambda411_opencore_wK1Fe_magcurve_T0.0_001_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "dark-violet" title "0 K"

set tmargin at screen bmgn+ysize
set bmargin at screen bmgn
set lmargin at screen lmgn
set rmargin at screen lmgn+xsize
set format x
set format y
unset label
set xlabel "{/Times-Italic B} [T] (|| {/Times-Italic a}-axis)" 0,0.2
set ylabel "{/Symbol-Oblique m}_0{/Times-Italic M} [T]" 1.0,4.8
set label "(c) {/Times-Italic T}=400 K" at lblx,lbly
set key reverse Left
set key at 27,0.65
plot \
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "black" title "SmFe_{12}",\
     "SmFe12H_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "red" title "SmFeH_{12}",\
     "SmFe12B_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "dark-green" title "SmFe_{12}B",\
     "SmFe12C_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "blue" title "SmFe_{12}C",\
     "SmFe12N_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "dark-violet" title "SmFe_{12}N"
unset ylabel

set tmargin at screen bmgn+ysize
set bmargin at screen bmgn
set lmargin at screen lmgn+xsize
set rmargin at screen rmgn
unset label
set label "(d) {/Times-Italic T}=400 K" at lblx,lbly
set xlabel "{/Times-Italic B} [T] (|| {/Times-Italic c}-axis)" 0,0.2
set format y ""
set key reverse Left
set key at 27,0.85
plot \
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T400.0_001_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "black" notitle "SmFe_{12}",\
     "SmFe12H_lambda411_opencore_wK1Fe_magcurve_T400.0_001_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "red" notitle "SmFe_{12}H",\
     "SmFe12B_lambda411_opencore_wK1Fe_magcurve_T400.0_001_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "dark-green" title "SmFe_{12}B",\
     "SmFe12C_lambda411_opencore_wK1Fe_magcurve_T400.0_001_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "blue" title "SmFe_{12}C",\
     "SmFe12N_lambda411_opencore_wK1Fe_magcurve_T400.0_001_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "dark-violet" title "SmFe_{12}N"

unset multiplot
reset






















set xrange [0:29]
set yrange [0:2.5]
set nokey
set ytics 0,1,3
set mxtics 2
set mytics 2
keyx=28
keyy=1.4
set size 0.5,0.65
set output "magcurve_SmFe12_SmFe12H_annum3.eps"
set multiplot
# set left and right margin
lmgn=0.06
#rmgn=0.4
rmgn=0.48
bmgn=0.08
tmgn=0.6
xsize=(rmgn-lmgn)
ysize=(tmgn-bmgn)/2.0
# set label position
lblx=0.5
lbly=2.2
set tmargin at screen bmgn+2*ysize
set bmargin at screen bmgn+ysize
set lmargin at screen lmgn
set rmargin at screen lmgn+xsize
set format x ""
set label "(a) SmFe_{12}" at lblx,lbly
set key 16,1.1
set key reverse Left
unset key
plot \
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T0.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 2 lw 2 lc rgb "blue" title "0 K",\
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T200.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 2 lw 2 lc rgb "dark-green" title "200 K",\
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T400.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 2 lw 2 lc rgb "red" title "400 K",\
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 2 lc rgb "blue" title "0 K",\
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T200.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 2 lc rgb "dark-green" title "200 K",\
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 2 lc rgb "red" title "400 K"
set ylabel "{/Symbol-Oblique m}_0{/Times-Italic M} [T]" 1.2,4.8
set format y
set format x
set tmargin at screen bmgn+ysize
set bmargin at screen bmgn
set lmargin at screen lmgn
set rmargin at screen lmgn+xsize
unset label
set label "(b) SmFe_{12}H" at lblx,lbly
set xlabel "{/Times-Italic B} [T]" 0,0.2
unset key
set label "400 K" tc rgb "red" at graph 0.37,0.5 font ",12"
set label "300 K" tc rgb "orange" at graph 0.51,0.56 font ",12"
set label "200 K" tc rgb "dark-yellow" at graph 0.64,0.65 font ",12"
set label "100 K" tc rgb "dark-green" at graph 0.75,0.74 font ",12"
set label "0 K" tc rgb "blue" at graph 0.88,0.77 font ",12"
plot \
     "SmFe12H_lambda411_opencore_wK1Fe_magcurve_T0.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 2 lw 2 lc rgb "blue" title "0 K",\
     "SmFe12H_lambda411_opencore_wK1Fe_magcurve_T200.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 2 lw 2 lc rgb "dark-green" title "200 K",\
     "SmFe12H_lambda411_opencore_wK1Fe_magcurve_T400.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 2 lw 2 lc rgb "red" title "400 K",\
     "SmFe12H_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 2 lc rgb "blue" title "0 K",\
     "SmFe12H_lambda411_opencore_wK1Fe_magcurve_T200.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 2 lc rgb "dark-green" title "200 K",\
     "SmFe12H_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 2 lc rgb "red" title "400 K"
unset multiplot
reset





set xrange [0:29]
set yrange [0:2.5]
set nokey
set ytics 0,1,3
set mxtics 2
set mytics 2
keyx=28
keyy=1.4
set size 0.5,0.65
set output "magcurve_SmFe12_SmFe12H_annum.eps"
set multiplot
# set left and right margin
lmgn=0.06
#rmgn=0.4
rmgn=0.48
bmgn=0.08
tmgn=0.6
xsize=(rmgn-lmgn)
ysize=(tmgn-bmgn)/2.0
# set label position
lblx=0.5
lbly=2.2
set tmargin at screen bmgn+2*ysize
set bmargin at screen bmgn+ysize
set lmargin at screen lmgn
set rmargin at screen lmgn+xsize
set format x ""
set label "(a) SmFe_{12}" at lblx,lbly
set key 16,1.1
set key reverse Left
unset key
plot \
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T0.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 2 lw 2 lc rgb "blue" title "0 K",\
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 2 lc rgb "blue" title "0 K",\
     "SmFe12H_lambda411_opencore_wK1Fe_magcurve_T0.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 2 lw 2 lc rgb "red" title "0 K",\
     "SmFe12H_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 2 lc rgb "red" title "0 K"#,\
#     "SmFe12B_lambda411_opencore_wK1Fe_magcurve_T0.0_001_mixexcf2_rd1.0.txt" u 1:2 w l lt 2 lw 2 lc rgb "dark-green" title "0 K",\
#     "SmFe12B_lambda411_opencore_wK1Fe_magcurve_T0.0_001_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 2 lc rgb "dark-green" title "0 K"

set ylabel "{/Symbol-Oblique m}_0{/Times-Italic M} [T]" 1.2,4.8
set format y
set format x
set tmargin at screen bmgn+ysize
set bmargin at screen bmgn
set lmargin at screen lmgn
set rmargin at screen lmgn+xsize
unset label
set label "(b) SmFe_{12}H" at lblx,lbly
set xlabel "{/Times-Italic B} [T]" 0,0.2
unset key
plot \
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T400.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 2 lw 2 lc rgb "blue" title "400 K",\
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 2 lc rgb "bluea" title "400 K",\
     "SmFe12H_lambda411_opencore_wK1Fe_magcurve_T400.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 2 lw 2 lc rgb "red" title "400 K",\
     "SmFe12H_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 2 lc rgb "red" title "400 K"#,\
#     "SmFe12B_lambda411_opencore_wK1Fe_magcurve_T400.0_001_mixexcf2_rd1.0.txt" u 1:2 w l lt 2 lw 2 lc rgb "dark-green" title "400 K",\
#     "SmFe12B_lambda411_opencore_wK1Fe_magcurve_T400.0_001_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 2 lc rgb "dark-green" title "400 K"

unset multiplot
reset





set xrange [0:20]
set yrange [0:2.5]
set nokey
set ytics 0,1,3
set mxtics 2
set mytics 2
keyx=28
keyy=1.4
set size 0.5,0.65
set output "magcurve_SmFe12_SmFe12H_annum2.eps"
set multiplot
# set left and right margin
lmgn=0.06
#rmgn=0.4
rmgn=0.48
bmgn=0.08
tmgn=0.6
xsize=(rmgn-lmgn)
ysize=(tmgn-bmgn)/2.0
# set label position
lblx=0.5
lbly=2.2
set tmargin at screen bmgn+2*ysize
set bmargin at screen bmgn+ysize
set lmargin at screen lmgn
set rmargin at screen lmgn+xsize
set format x ""
set label "(a) SmFe_{12}" at lblx,lbly
set key 16,1.1
set key reverse Left
unset key
plot \
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T0.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 2 lw 2 lc rgb "dark-violet" title "0 K",\
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T100.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 2 lw 2 lc rgb "dark-green" title "100 K",\
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T200.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 2 lw 2 lc rgb "dark-yellow" title "200 K",\
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T300.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 2 lw 2 lc rgb "orange" title "300 K",\
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T400.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 2 lw 2 lc rgb "red" title "400 K",\
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 2 lc rgb "dark-violet" title "0 K",\
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T100.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 2 lc rgb "dark-green" title "100 K",\
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T200.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 2 lc rgb "dark-yellow" title "200 K",\
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T300.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 2 lc rgb "orange" title "300 K",\
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 2 lc rgb "red" title "400 K"
set ylabel "{/Symbol-Oblique m}_0{/Times-Italic M} [T]" 1.2,4.8
set format y
set format x
set tmargin at screen bmgn+ysize
set bmargin at screen bmgn
set lmargin at screen lmgn
set rmargin at screen lmgn+xsize
unset label
set label "(b) SmFe_{12}H" at lblx,lbly
set xlabel "{/Times-Italic B} [T]" 0,0.2
unset key
set label "400 K" tc rgb "red" at graph 0.37,0.5 font ",12"
set label "300 K" tc rgb "orange" at graph 0.51,0.56 font ",12"
set label "200 K" tc rgb "dark-yellow" at graph 0.64,0.65 font ",12"
set label "100 K" tc rgb "dark-green" at graph 0.75,0.74 font ",12"
set label "0 K" tc rgb "blue" at graph 0.88,0.77 font ",12"
plot \
     "SmFe12H_lambda411_opencore_wK1Fe_magcurve_T0.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 2 lw 2 lc rgb "dark-violet" title "0 K",\
     "SmFe12H_lambda411_opencore_wK1Fe_magcurve_T100.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 2 lw 2 lc rgb "dark-green" title "100 K",\
     "SmFe12H_lambda411_opencore_wK1Fe_magcurve_T200.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 2 lw 2 lc rgb "dark-yellow" title "200 K",\
     "SmFe12H_lambda411_opencore_wK1Fe_magcurve_T300.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 2 lw 2 lc rgb "orange" title "300 K",\
     "SmFe12H_lambda411_opencore_wK1Fe_magcurve_T400.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 2 lw 2 lc rgb "red" title "400 K",\
     "SmFe12H_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 2 lc rgb "dark-violet" title "0 K",\
     "SmFe12H_lambda411_opencore_wK1Fe_magcurve_T100.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 2 lc rgb "dark-green" title "100 K",\
     "SmFe12H_lambda411_opencore_wK1Fe_magcurve_T200.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 2 lc rgb "dark-yellow" title "200 K",\
     "SmFe12H_lambda411_opencore_wK1Fe_magcurve_T300.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 2 lc rgb "orange" title "300 K",\
     "SmFe12H_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 2 lc rgb "red" title "400 K"
unset multiplot
reset





set output "magcurve_SmFe12_lambda411_wK1Fe_an.eps"
set size 0.5,0.5
set xrange [0:10]
set yrange [0:2.2]
set ytics 0,0.5,2.5
set mxtics 2
set mytics 2
set xlabel "Applied Field [T]" 
set ylabel "Magnetization [T]" 2,0
set label "400 K" tc rgb "red" at 3.5,1.37 font ",12"
set label "300 K" tc rgb "orange" at 4.5,1.58 font ",12"
set label "200 K" tc rgb "dark-yellow" at 5,1.73 font ",12"
set label "100 K" tc rgb "dark-green" at 5.2,1.84 font ",12"
set label "0 K" tc rgb "blue" at 6.7,1.89 font ",12"
set nokey
plot \
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T0.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 1 lw 2 lc rgb "blue" title "0 K",\
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T100.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 1 lw 2 lc rgb "dark-green" title "100 K",\
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T200.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 1 lw 2 lc rgb "dark-yellow" title "200 K",\
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T300.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 1 lw 2 lc rgb "orange" title "300 K",\
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T400.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 1 lw 2 lc rgb "red" title "400 K"
reset








set output "magcurve_SmFe12_lambda411_wK1Fe_num.eps"
set size 0.5,0.5
set xrange [0:10]
set yrange [0:2.2]
set ytics 0,0.5,2.5
set mxtics 2
set mytics 2
set xlabel "Applied Field [T]" 
set ylabel "Magnetization [T]" 2,0
set label "400 K" tc rgb "red" at 3.5,1.37 font ",12"
set label "300 K" tc rgb "orange" at 4.5,1.58 font ",12"
set label "200 K" tc rgb "dark-yellow" at 5,1.73 font ",12"
set label "100 K" tc rgb "dark-green" at 5.5,1.84 font ",12"
set label "0 K" tc rgb "blue" at 6.7,1.89 font ",12"
set nokey
plot \
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 2 lc rgb "blue" title "0 K",\
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T100.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 2 lc rgb "dark-green" title "100 K",\
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T200.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 2 lc rgb "dark-yellow" title "200 K",\
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T300.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 2 lc rgb "orange" title "300 K",\
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 2 lc rgb "red" title "400 K"
reset


set xrange [0:17]
set yrange [0:2.5]
set nokey
set ytics 0,1,3
set mxtics 2
set mytics 2
keyx=28
keyy=1.4
set size 0.5,0.65
set output "magcurve_SmFe12H_lambda411_wK1Fe.eps"
set multiplot
# set left and right margin
lmgn=0.06
#rmgn=0.4
rmgn=0.48
bmgn=0.08
tmgn=0.6
xsize=(rmgn-lmgn)
ysize=(tmgn-bmgn)/2.0
# set label position
lblx=0.5
lbly=2.2
set tmargin at screen bmgn+2*ysize
set bmargin at screen bmgn+ysize
set lmargin at screen lmgn
set rmargin at screen lmgn+xsize
set format x ""
set label "(a) SmFe_{12}" at lblx,lbly
set key 16,1.1
set key reverse Left
unset key
plot \
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 2 lc rgb "blue" title "0 K",\
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T100.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 2 lc rgb "dark-green" title "100 K",\
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T200.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 2 lc rgb "dark-yellow" title "200 K",\
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T300.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 2 lc rgb "orange" title "300 K",\
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 2 lc rgb "red" title "400 K"
set ylabel "{/Symbol-Oblique m}_0{/Times-Italic M} [T]" 1.2,4.8
set format y
set format x
set tmargin at screen bmgn+ysize
set bmargin at screen bmgn
set lmargin at screen lmgn
set rmargin at screen lmgn+xsize
unset label
set label "(b) SmFe_{12}H" at lblx,lbly
set xlabel "{/Times-Italic B} [T]" 0,0.2
unset key
set label "400 K" tc rgb "red" at graph 0.37,0.5 font ",12"
set label "300 K" tc rgb "orange" at graph 0.51,0.56 font ",12"
set label "200 K" tc rgb "dark-yellow" at graph 0.64,0.65 font ",12"
set label "100 K" tc rgb "dark-green" at graph 0.75,0.74 font ",12"
set label "0 K" tc rgb "blue" at graph 0.88,0.77 font ",12"
plot \
     "SmFe12H_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 2 lc rgb "blue" title "0 K",\
     "SmFe12H_lambda411_opencore_wK1Fe_magcurve_T100.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 2 lc rgb "dark-green" title "100 K",\
     "SmFe12H_lambda411_opencore_wK1Fe_magcurve_T200.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 2 lc rgb "dark-yellow" title "200 K",\
     "SmFe12H_lambda411_opencore_wK1Fe_magcurve_T300.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 2 lc rgb "orange" title "300 K",\
     "SmFe12H_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 2 lc rgb "red" title "400 K"
unset multiplot
reset


