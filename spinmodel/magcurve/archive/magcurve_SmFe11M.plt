set terminal postscript eps enhanced color font "Times-Roman,12.5"




set xrange [0:29]
set yrange [0:2.5]
set nokey
set ytics 0,1,2
set mxtics 2
set mytics 2
set size 0.5,0.75
set output "magcurve_SmFe11M_annum.eps"

VDiop=8.568*8.568*4.797 # 1.0d-30
mu0=4.0*pi # 1.0d-7
muB=9.2740100783 # 1.0d-24
coffDiop=muB*mu0*0.1/(VDiop/2.0)

set multiplot
# set left and right margin
lmgn=0.057
rmgn=0.49
bmgn=0.08
tmgn=0.74
xsize=(rmgn-lmgn)/2.0
ysize=(tmgn-bmgn)/3.0
# set label position
lblx=1
lbly=2.15
set tmargin at screen bmgn+3*ysize
set bmargin at screen bmgn+2*ysize
set lmargin at screen lmgn
set rmargin at screen lmgn+xsize
set format x ""
set label "(a) SmFe_{11}Ti   {/Times-Italic T}=0 K" at lblx,lbly
set xtics scale 3.0/4.0,3.0/8.0
set ytics scale 3.0/4.0,3.0/8.0
#set label "[100]"  at 8.5,0.35
plot \
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T0.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 1 lw 1 lc rgb "black" title "SmFe_{12}",\
     "SmFe11Ti_f_fullexp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 2 lw 1.5 lc rgb "red" title "f-site",\
     "SmFe11Ti_i_fullexp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 2 lw 1.5 lc rgb "dark-green" title "i-site",\
     "SmFe11Ti_j_fullexp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 2 lw 1.5 lc rgb "blue" title "j-site",\
     "SmFe11Ti_f_fullexp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 1 lw 3 lc rgb "red" title "f-site",\
     "SmFe11Ti_i_fullexp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 1 lw 3 lc rgb "dark-green" title "i-site",\
     "SmFe11Ti_j_fullexp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 1 lw 3 lc rgb "blue" title "j-site"#,\
#     "magcurve_SmFe11Ti_Diop.txt" index 0 u 1:($2*coffDiop) w lp pt 2 ps 0.5 lc rgb "black"
unset label
set tmargin at screen bmgn+3*ysize
set bmargin at screen bmgn+2*ysize
set lmargin at screen lmgn+xsize
set rmargin at screen lmgn+2*xsize
set format x ""
set format y ""
set label "(b) SmFe_{11}Ti   {/Times-Italic T}=400 K"  at lblx,lbly
#set label "[100]"  at 6.8,0.35
plot \
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T400.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 1 lw 1 lc rgb "black" title "SmFe_{12}",\
     "SmFe11Ti_f_fullexp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 2 lw 1.5 lc rgb "red" title "f-site",\
     "SmFe11Ti_i_fullexp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 2 lw 1.5 lc rgb "dark-green" title "i-site",\
     "SmFe11Ti_j_fullexp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 2 lw 1.5 lc rgb "blue" title "j-site",\
     "SmFe11Ti_f_fullexp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 1 lw 3 lc rgb "red" title "f-site",\
     "SmFe11Ti_i_fullexp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 1 lw 3 lc rgb "dark-green" title "i-site",\
     "SmFe11Ti_j_fullexp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 1 lw 3 lc rgb "blue" title "j-site"#,\
#     "magcurve_SmFe11Ti_Diop.txt" index 1 u 1:($2*coffDiop) w lp pt 2 ps 0.5 lc rgb "black"
unset label
#set ylabel "Magnetization [T]" 0.5
set ylabel "{/Symbol-Oblique m}_0{/Times-BoldItalic M}_{s}({/Times-Italic T}){/Symbol \327}({/Times-BoldItalic B}/{/Times-Italic B}) [T]" 0.5
set format y
set tmargin at screen bmgn+2*ysize
set bmargin at screen bmgn+ysize
set lmargin at screen lmgn
set rmargin at screen lmgn+xsize
set label "(c) SmFe_{11}V    {/Times-Italic T}=0 K"  at lblx,lbly
#set label "[100]"  at 10.5,0.35
plot \
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T0.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 1 lw 1 lc rgb "black" title "SmFe_{12}",\
     "SmFe11V_f_fullexp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 2 lw 1.5 lc rgb "red" title "f-site",\
     "SmFe11V_i_fullexp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 2 lw 1.5 lc rgb "dark-green" title "i-site",\
     "SmFe11V_j_fullexp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 2 lw 1.5 lc rgb "blue" title "j-site",\
     "SmFe11V_f_fullexp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 1 lw 3 lc rgb "red" title "f-site",\
     "SmFe11V_i_fullexp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 1 lw 3 lc rgb "dark-green" title "i-site",\
     "SmFe11V_j_fullexp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 1 lw 3 lc rgb "blue" title "j-site"
unset label
unset ylabel
set tmargin at screen bmgn+2*ysize
set bmargin at screen bmgn+ysize
set lmargin at screen lmgn+xsize
set rmargin at screen lmgn+2*xsize
set format x ""
set format y ""
set label "(d) SmFe_{11}V   {/Times-Italic T}=400 K"  at lblx,lbly
#set label "[100]"  at 9,0.35
plot \
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T400.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 1 lw 1 lc rgb "black" title "SmFe_{12}",\
     "SmFe11V_f_fullexp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 2 lw 1.5 lc rgb "red" title "f-site",\
     "SmFe11V_i_fullexp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 2 lw 1.5 lc rgb "dark-green" title "i-site",\
     "SmFe11V_j_fullexp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 2 lw 1.5 lc rgb "blue" title "j-site",\
     "SmFe11V_f_fullexp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 1 lw 3 lc rgb "red" title "f-site",\
     "SmFe11V_i_fullexp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 1 lw 3 lc rgb "dark-green" title "i-site",\
     "SmFe11V_j_fullexp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 1 lw 3 lc rgb "blue" title "j-site"
unset label
set format x
set format y
set xtics 
set tmargin at screen bmgn+ysize
set bmargin at screen bmgn
set lmargin at screen lmgn
set rmargin at screen lmgn+xsize
set key reverse Left
keyx=30
keyy=1.3
set key keyx,keyy
set label "(e) SmFe_{11}Co    {/Times-Italic T}=0 K"  at lblx,lbly
#set label "[100]"  at 5,0.35
plot \
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T0.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 1 lw 1 lc rgb "black" title "SmFe_{12}",\
     "SmFe11Co_f_fullexp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 2 lw 1.5 lc rgb "red" notitle "8{/Times-Italic f}-site",\
     "SmFe11Co_i_fullexp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 2 lw 1.5 lc rgb "dark-green" notitle "8{/Times-Italic i}-site",\
     "SmFe11Co_j_fullexp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 2 lw 1.5 lc rgb "blue" notitle "8{/Times-Italic j}-site",\
     "SmFe11Co_f_fullexp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 1 lw 3 lc rgb "red" title "8{/Times-Italic f}-site",\
     "SmFe11Co_i_fullexp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 1 lw 3 lc rgb "dark-green" title "8{/Times-Italic i}-site",\
     "SmFe11Co_j_fullexp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 1 lw 3 lc rgb "blue" title "8{/Times-Italic j}-site"
unset label
unset key
set tmargin at screen bmgn+ysize
set bmargin at screen bmgn
set lmargin at screen lmgn+xsize
set rmargin at screen lmgn+2*xsize
set format y ""
set label "(f) SmFe_{11}Co   {/Times-Italic T}=400 K"  at lblx,lbly
#set label "[100]"  at 3.2,0.35
set xlabel "{/Times-Italic B} [T] (|| {/Times-Italic a}-axis)"  offset screen -0.11,0
plot \
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T400.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 1 lw 1 lc rgb "black" \
     title "{/Times-Roman=12 SmFe_{12}}",\
     "SmFe11Co_f_fullexp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 2 lw 1.5 lc rgb "red"\
     title "{/Times-Roman=15 8{/Times-Italic f}-site}",\
     "SmFe11Co_i_fullexp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 2 lw 1.5 lc rgb "dark-green" \
     title "{/Times-Roman=15 8{/Times-Italic i}-site}",\
     "SmFe11Co_j_fullexp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 2 lw 1.5 lc rgb "blue" \
     title "{/Times-Roman=15 8{/Times-Italic j}-site}",\
     "SmFe11Co_f_fullexp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 1 lw 3 lc rgb "red"\
     title "{/Times-Roman=15 8{/Times-Italic f}-site}",\
     "SmFe11Co_i_fullexp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 1 lw 3 lc rgb "dark-green" \
     title "{/Times-Roman=15 8{/Times-Italic i}-site}",\
     "SmFe11Co_j_fullexp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 1 lw 3 lc rgb "blue" \
     title "{/Times-Roman=15 8{/Times-Italic j}-site}"
unset multiplot
reset








set xrange [0:29]
set yrange [0:2.5]
set nokey
set ytics 0,1,2
set mxtics 2
set mytics 2
set size 0.5,0.75
set output "magcurve_SmFe11M_annum_wexp.eps"

VDiop=8.568*8.568*4.797 # 1.0d-30
mu0=4.0*pi # 1.0d-7
muB=9.2740100783 # 1.0d-24
coffDiop=muB*mu0*0.1/(VDiop/2.0)

set multiplot
# set left and right margin
lmgn=0.057
rmgn=0.49
bmgn=0.08
tmgn=0.74
xsize=(rmgn-lmgn)/2.0
ysize=(tmgn-bmgn)/3.0
# set label position
lblx=1
lbly=2.15
set tmargin at screen bmgn+3*ysize
set bmargin at screen bmgn+2*ysize
set lmargin at screen lmgn
set rmargin at screen lmgn+xsize
set format x ""
set label "(a) SmFe_{11}Ti   {/Times-Italic T}=0 K" at lblx,lbly
set xtics scale 3.0/4.0,3.0/8.0
set ytics scale 3.0/4.0,3.0/8.0
#set label "[100]"  at 8.5,0.35
plot \
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T0.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 1 lw 1 lc rgb "black" title "SmFe_{12}",\
     "SmFe11Ti_f_fullexp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 2 lw 1.5 lc rgb "red" title "f-site",\
     "SmFe11Ti_i_fullexp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 2 lw 1.5 lc rgb "dark-green" title "i-site",\
     "SmFe11Ti_j_fullexp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 2 lw 1.5 lc rgb "blue" title "j-site",\
     "SmFe11Ti_f_fullexp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 1 lw 3 lc rgb "red" title "f-site",\
     "SmFe11Ti_i_fullexp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 1 lw 3 lc rgb "dark-green" title "i-site",\
     "SmFe11Ti_j_fullexp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 1 lw 3 lc rgb "blue" title "j-site",\
     "magcurve_SmFe11Ti_Diop.txt" index 0 u 1:($2*coffDiop) w lp pt 2 ps 0.5 lc rgb "black"
unset label
set tmargin at screen bmgn+3*ysize
set bmargin at screen bmgn+2*ysize
set lmargin at screen lmgn+xsize
set rmargin at screen lmgn+2*xsize
set format x ""
set format y ""
set label "(b) SmFe_{11}Ti   {/Times-Italic T}=300 K"  at lblx,lbly
#set label "[100]"  at 6.8,0.35
plot \
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T300.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 1 lw 1 lc rgb "black" title "SmFe_{12}",\
     "SmFe11Ti_f_fullexp_lambda411_opencore_wK1Fe_magcurve_T300.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 2 lw 1.5 lc rgb "red" title "f-site",\
     "SmFe11Ti_i_fullexp_lambda411_opencore_wK1Fe_magcurve_T300.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 2 lw 1.5 lc rgb "dark-green" title "i-site",\
     "SmFe11Ti_j_fullexp_lambda411_opencore_wK1Fe_magcurve_T300.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 2 lw 1.5 lc rgb "blue" title "j-site",\
     "SmFe11Ti_f_fullexp_lambda411_opencore_wK1Fe_magcurve_T300.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 1 lw 3 lc rgb "red" title "f-site",\
     "SmFe11Ti_i_fullexp_lambda411_opencore_wK1Fe_magcurve_T300.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 1 lw 3 lc rgb "dark-green" title "i-site",\
     "SmFe11Ti_j_fullexp_lambda411_opencore_wK1Fe_magcurve_T300.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 1 lw 3 lc rgb "blue" title "j-site",\
     "magcurve_SmFe11Ti_Diop.txt" index 1 u 1:($2*coffDiop) w lp pt 2 ps 0.5 lc rgb "black"
unset label
#set ylabel "Magnetization [T]" 0.5
set ylabel "{/Symbol-Oblique m}_0{/Times-BoldItalic M}_{s}({/Times-Italic T}){/Symbol \327}({/Times-BoldItalic B}/{/Times-Italic B}) [T]" 0.5
set format y
set tmargin at screen bmgn+2*ysize
set bmargin at screen bmgn+ysize
set lmargin at screen lmgn
set rmargin at screen lmgn+xsize
set label "(c) SmFe_{11}V    {/Times-Italic T}=0 K"  at lblx,lbly
#set label "[100]"  at 10.5,0.35
plot \
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T0.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 1 lw 1 lc rgb "black" title "SmFe_{12}",\
     "SmFe11V_f_fullexp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 2 lw 1.5 lc rgb "red" title "f-site",\
     "SmFe11V_i_fullexp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 2 lw 1.5 lc rgb "dark-green" title "i-site",\
     "SmFe11V_j_fullexp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 2 lw 1.5 lc rgb "blue" title "j-site",\
     "SmFe11V_f_fullexp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 1 lw 3 lc rgb "red" title "f-site",\
     "SmFe11V_i_fullexp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 1 lw 3 lc rgb "dark-green" title "i-site",\
     "SmFe11V_j_fullexp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 1 lw 3 lc rgb "blue" title "j-site"
unset label
unset ylabel
set tmargin at screen bmgn+2*ysize
set bmargin at screen bmgn+ysize
set lmargin at screen lmgn+xsize
set rmargin at screen lmgn+2*xsize
set format x ""
set format y ""
set label "(d) SmFe_{11}V   {/Times-Italic T}=300 K"  at lblx,lbly
#set label "[100]"  at 9,0.35
plot \
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T300.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 1 lw 1 lc rgb "black" title "SmFe_{12}",\
     "SmFe11V_f_fullexp_lambda411_opencore_wK1Fe_magcurve_T300.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 2 lw 1.5 lc rgb "red" title "f-site",\
     "SmFe11V_i_fullexp_lambda411_opencore_wK1Fe_magcurve_T300.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 2 lw 1.5 lc rgb "dark-green" title "i-site",\
     "SmFe11V_j_fullexp_lambda411_opencore_wK1Fe_magcurve_T300.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 2 lw 1.5 lc rgb "blue" title "j-site",\
     "SmFe11V_f_fullexp_lambda411_opencore_wK1Fe_magcurve_T300.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 1 lw 3 lc rgb "red" title "f-site",\
     "SmFe11V_i_fullexp_lambda411_opencore_wK1Fe_magcurve_T300.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 1 lw 3 lc rgb "dark-green" title "i-site",\
     "SmFe11V_j_fullexp_lambda411_opencore_wK1Fe_magcurve_T300.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 1 lw 3 lc rgb "blue" title "j-site"
unset label
set format x
set format y
set xtics 
set tmargin at screen bmgn+ysize
set bmargin at screen bmgn
set lmargin at screen lmgn
set rmargin at screen lmgn+xsize
set key reverse Left
keyx=30
keyy=1.3
set key keyx,keyy
set label "(e) SmFe_{11}Co    {/Times-Italic T}=0 K"  at lblx,lbly
#set label "[100]"  at 5,0.35
plot \
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T0.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 1 lw 1 lc rgb "black" title "SmFe_{12}",\
     "SmFe11Co_f_fullexp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 2 lw 1.5 lc rgb "red" notitle "8{/Times-Italic f}-site",\
     "SmFe11Co_i_fullexp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 2 lw 1.5 lc rgb "dark-green" notitle "8{/Times-Italic i}-site",\
     "SmFe11Co_j_fullexp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 2 lw 1.5 lc rgb "blue" notitle "8{/Times-Italic j}-site",\
     "SmFe11Co_f_fullexp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 1 lw 3 lc rgb "red" title "8{/Times-Italic f}-site",\
     "SmFe11Co_i_fullexp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 1 lw 3 lc rgb "dark-green" title "8{/Times-Italic i}-site",\
     "SmFe11Co_j_fullexp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 1 lw 3 lc rgb "blue" title "8{/Times-Italic j}-site"
unset label
unset key
set tmargin at screen bmgn+ysize
set bmargin at screen bmgn
set lmargin at screen lmgn+xsize
set rmargin at screen lmgn+2*xsize
set format y ""
set label "(f) SmFe_{11}Co   {/Times-Italic T}=300 K"  at lblx,lbly
#set label "[100]"  at 3.2,0.35
set xlabel "{/Times-Italic B} [T] (|| {/Times-Italic a}-axis)"  offset screen -0.11,0
plot \
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T300.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 1 lw 1 lc rgb "black" \
     title "{/Times-Roman=12 SmFe_{12}}",\
     "SmFe11Co_f_fullexp_lambda411_opencore_wK1Fe_magcurve_T300.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 2 lw 1.5 lc rgb "red"\
     title "{/Times-Roman=15 8{/Times-Italic f}-site}",\
     "SmFe11Co_i_fullexp_lambda411_opencore_wK1Fe_magcurve_T300.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 2 lw 1.5 lc rgb "dark-green" \
     title "{/Times-Roman=15 8{/Times-Italic i}-site}",\
     "SmFe11Co_j_fullexp_lambda411_opencore_wK1Fe_magcurve_T300.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 2 lw 1.5 lc rgb "blue" \
     title "{/Times-Roman=15 8{/Times-Italic j}-site}",\
     "SmFe11Co_f_fullexp_lambda411_opencore_wK1Fe_magcurve_T300.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 1 lw 3 lc rgb "red"\
     title "{/Times-Roman=15 8{/Times-Italic f}-site}",\
     "SmFe11Co_i_fullexp_lambda411_opencore_wK1Fe_magcurve_T300.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 1 lw 3 lc rgb "dark-green" \
     title "{/Times-Roman=15 8{/Times-Italic i}-site}",\
     "SmFe11Co_j_fullexp_lambda411_opencore_wK1Fe_magcurve_T300.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 1 lw 3 lc rgb "blue" \
     title "{/Times-Roman=15 8{/Times-Italic j}-site}"
unset multiplot
reset





set xrange [0:29]
set yrange [0:2.5]
set nokey
set ytics 0,1,2
set mxtics 2
set mytics 2
set size 0.5,0.75
set output "magcurve_SmFe11M_annum2.eps"

VDiop=8.568*8.568*4.797 # 1.0d-30
mu0=4.0*pi # 1.0d-7
muB=9.2740100783 # 1.0d-24
coffDiop=muB*mu0*0.1/(VDiop/2.0)

set multiplot
# set left and right margin
lmgn=0.057
rmgn=0.49
bmgn=0.08
tmgn=0.74
xsize=(rmgn-lmgn)/2.0
ysize=(tmgn-bmgn)/3.0
# set label position
lblx=1
lbly=2.15
set tmargin at screen bmgn+3*ysize
set bmargin at screen bmgn+2*ysize
set lmargin at screen lmgn
set rmargin at screen lmgn+xsize
set format x ""
set label "(a) SmFe_{11}Ti   {/Times-Italic T}=0 K" at lblx,lbly
set xtics scale 3.0/4.0,3.0/8.0
set ytics scale 3.0/4.0,3.0/8.0
#set label "[100]"  at 8.5,0.35
plot \
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 1 lc rgb "black" title "SmFe_{12}",\
     "SmFe11Ti_f_fullexp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 2 lw 1.5 lc rgb "red" title "f-site",\
     "SmFe11Ti_i_fullexp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 2 lw 1.5 lc rgb "dark-green" title "i-site",\
     "SmFe11Ti_j_fullexp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 2 lw 1.5 lc rgb "blue" title "j-site",\
     "SmFe11Ti_f_fullexp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "red" title "f-site",\
     "SmFe11Ti_i_fullexp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "dark-green" title "i-site",\
     "SmFe11Ti_j_fullexp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "blue" title "j-site",\
     "magcurve_SmFe11Ti_Diop.txt" u 1:($2*coffDiop) w lp
unset label
set tmargin at screen bmgn+3*ysize
set bmargin at screen bmgn+2*ysize
set lmargin at screen lmgn+xsize
set rmargin at screen lmgn+2*xsize
set format x ""
set format y ""
set label "(b) SmFe_{11}Ti   {/Times-Italic T}=400 K"  at lblx,lbly
#set label "[100]"  at 6.8,0.35
plot \
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 1 lc rgb "black" title "SmFe_{12}",\
     "SmFe11Ti_f_fullexp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 2 lw 1.5 lc rgb "red" title "f-site",\
     "SmFe11Ti_i_fullexp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 2 lw 1.5 lc rgb "dark-green" title "i-site",\
     "SmFe11Ti_j_fullexp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 2 lw 1.5 lc rgb "blue" title "j-site",\
     "SmFe11Ti_f_fullexp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "red" title "f-site",\
     "SmFe11Ti_i_fullexp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "dark-green" title "i-site",\
     "SmFe11Ti_j_fullexp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "blue" title "j-site"
unset label
#set ylabel "Magnetization [T]" 0.5
set ylabel "{/Symbol-Oblique m}_0{/Times-BoldItalic M}_{s}({/Times-Italic T}){/Symbol \327}({/Times-BoldItalic B}/{/Times-Italic B}) [T]" 0.5
set format y
set tmargin at screen bmgn+2*ysize
set bmargin at screen bmgn+ysize
set lmargin at screen lmgn
set rmargin at screen lmgn+xsize
set label "(c) SmFe_{11}V    {/Times-Italic T}=0 K"  at lblx,lbly
#set label "[100]"  at 10.5,0.35
plot \
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 1 lc rgb "black" title "SmFe_{12}",\
     "SmFe11V_f_fullexp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 2 lw 1.5 lc rgb "red" title "f-site",\
     "SmFe11V_i_fullexp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 2 lw 1.5 lc rgb "dark-green" title "i-site",\
     "SmFe11V_j_fullexp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 2 lw 1.5 lc rgb "blue" title "j-site",\
     "SmFe11V_f_fullexp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "red" title "f-site",\
     "SmFe11V_i_fullexp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "dark-green" title "i-site",\
     "SmFe11V_j_fullexp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "blue" title "j-site"
unset label
unset ylabel
set tmargin at screen bmgn+2*ysize
set bmargin at screen bmgn+ysize
set lmargin at screen lmgn+xsize
set rmargin at screen lmgn+2*xsize
set format x ""
set format y ""
set label "(d) SmFe_{11}V   {/Times-Italic T}=400 K"  at lblx,lbly
#set label "[100]"  at 9,0.35
plot \
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 1 lc rgb "black" title "SmFe_{12}",\
     "SmFe11V_f_fullexp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 2 lw 1.5 lc rgb "red" title "f-site",\
     "SmFe11V_i_fullexp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 2 lw 1.5 lc rgb "dark-green" title "i-site",\
     "SmFe11V_j_fullexp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 2 lw 1.5 lc rgb "blue" title "j-site",\
     "SmFe11V_f_fullexp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "red" title "f-site",\
     "SmFe11V_i_fullexp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "dark-green" title "i-site",\
     "SmFe11V_j_fullexp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "blue" title "j-site"
unset label
set format x
set format y
set xtics 
set tmargin at screen bmgn+ysize
set bmargin at screen bmgn
set lmargin at screen lmgn
set rmargin at screen lmgn+xsize
set key reverse Left
keyx=30
keyy=1.3
set key keyx,keyy
set label "(e) SmFe_{11}Co    {/Times-Italic T}=0 K"  at lblx,lbly
#set label "[100]"  at 5,0.35
plot \
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 1 lc rgb "black" title "SmFe_{12}",\
     "SmFe11Co_f_fullexp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 2 lw 1.5 lc rgb "red" notitle "8{/Times-Italic f}-site",\
     "SmFe11Co_i_fullexp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 2 lw 1.5 lc rgb "dark-green" notitle "8{/Times-Italic i}-site",\
     "SmFe11Co_j_fullexp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 2 lw 1.5 lc rgb "blue" notitle "8{/Times-Italic j}-site",\
     "SmFe11Co_f_fullexp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "red" title "8{/Times-Italic f}-site",\
     "SmFe11Co_i_fullexp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "dark-green" title "8{/Times-Italic i}-site",\
     "SmFe11Co_j_fullexp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "blue" title "8{/Times-Italic j}-site"
unset label
unset key
set tmargin at screen bmgn+ysize
set bmargin at screen bmgn
set lmargin at screen lmgn+xsize
set rmargin at screen lmgn+2*xsize
set format y ""
set label "(f) SmFe_{11}Co   {/Times-Italic T}=400 K"  at lblx,lbly
#set label "[100]"  at 3.2,0.35
set xlabel "{/Times-Italic B} [T] (|| {/Times-Italic a}-axis)"  offset screen -0.11,0
plot \
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 1 lc rgb "black" \
     title "{/Times-Roman=12 SmFe_{12}}",\
     "SmFe11Co_f_fullexp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 2 lw 1.5 lc rgb "red"\
     title "{/Times-Roman=15 8{/Times-Italic f}-site}",\
     "SmFe11Co_i_fullexp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 2 lw 1.5 lc rgb "dark-green" \
     title "{/Times-Roman=15 8{/Times-Italic i}-site}",\
     "SmFe11Co_j_fullexp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 2 lw 1.5 lc rgb "blue" \
     title "{/Times-Roman=15 8{/Times-Italic j}-site}",\
     "SmFe11Co_f_fullexp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "red"\
     title "{/Times-Roman=15 8{/Times-Italic f}-site}",\
     "SmFe11Co_i_fullexp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "dark-green" \
     title "{/Times-Roman=15 8{/Times-Italic i}-site}",\
     "SmFe11Co_j_fullexp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "blue" \
     title "{/Times-Roman=15 8{/Times-Italic j}-site}"
unset multiplot
reset










set xrange [0:29]
set yrange [0:2.5]
set nokey
set ytics 0,1,2
set mxtics 2
set mytics 2
set size 0.5,0.75
set output "magcurve_SmFe11M_annum_m0m4.eps"
set multiplot
# set left and right margin
lmgn=0.057
rmgn=0.49
bmgn=0.08
tmgn=0.74
xsize=(rmgn-lmgn)/2.0
ysize=(tmgn-bmgn)/3.0
# set label position
lblx=1
lbly=2.15
set tmargin at screen bmgn+3*ysize
set bmargin at screen bmgn+2*ysize
set lmargin at screen lmgn
set rmargin at screen lmgn+xsize
set format x ""
set label "(a) SmFe_{11}Ti   {/Times-Italic T}=0 K" at lblx,lbly
set xtics scale 3.0/4.0,3.0/8.0
set ytics scale 3.0/4.0,3.0/8.0
#set label "[100]"  at 8.5,0.35
plot \
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 1 lc rgb "black" title "SmFe_{12}",\
     "SmFe11Ti_f_m0m4exp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 2 lw 1.5 lc rgb "red" title "f-site",\
     "SmFe11Ti_i_m0m4exp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 2 lw 1.5 lc rgb "dark-green" title "i-site",\
     "SmFe11Ti_j_m0m4exp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 2 lw 1.5 lc rgb "blue" title "j-site",\
     "SmFe11Ti_f_fullexp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "red" title "f-site",\
     "SmFe11Ti_i_fullexp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "dark-green" title "i-site",\
     "SmFe11Ti_j_fullexp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "blue" title "j-site"
unset label
set tmargin at screen bmgn+3*ysize
set bmargin at screen bmgn+2*ysize
set lmargin at screen lmgn+xsize
set rmargin at screen lmgn+2*xsize
set format x ""
set format y ""
set label "(b) SmFe_{11}Ti   {/Times-Italic T}=400 K"  at lblx,lbly
#set label "[100]"  at 6.8,0.35
plot \
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 1 lc rgb "black" title "SmFe_{12}",\
     "SmFe11Ti_f_m0m4exp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 2 lw 1.5 lc rgb "red" title "f-site",\
     "SmFe11Ti_i_m0m4exp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 2 lw 1.5 lc rgb "dark-green" title "i-site",\
     "SmFe11Ti_j_m0m4exp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 2 lw 1.5 lc rgb "blue" title "j-site",\
     "SmFe11Ti_f_fullexp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "red" title "f-site",\
     "SmFe11Ti_i_fullexp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "dark-green" title "i-site",\
     "SmFe11Ti_j_fullexp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "blue" title "j-site"
unset label
set ylabel "Magnetization [T]" 0.5
set format y
set tmargin at screen bmgn+2*ysize
set bmargin at screen bmgn+ysize
set lmargin at screen lmgn
set rmargin at screen lmgn+xsize
set label "(c) SmFe_{11}V    {/Times-Italic T}=0 K"  at lblx,lbly
#set label "[100]"  at 10.5,0.35
plot \
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 1 lc rgb "black" title "SmFe_{12}",\
     "SmFe11V_f_m0m4exp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 2 lw 1.5 lc rgb "red" title "f-site",\
     "SmFe11V_i_m0m4exp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 2 lw 1.5 lc rgb "dark-green" title "i-site",\
     "SmFe11V_j_m0m4exp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 2 lw 1.5 lc rgb "blue" title "j-site",\
     "SmFe11V_f_fullexp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "red" title "f-site",\
     "SmFe11V_i_fullexp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "dark-green" title "i-site",\
     "SmFe11V_j_fullexp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "blue" title "j-site"
unset label
unset ylabel
set tmargin at screen bmgn+2*ysize
set bmargin at screen bmgn+ysize
set lmargin at screen lmgn+xsize
set rmargin at screen lmgn+2*xsize
set format x ""
set format y ""
set label "(d) SmFe_{11}V   {/Times-Italic T}=400 K"  at lblx,lbly
#set label "[100]"  at 9,0.35
plot \
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 1 lc rgb "black" title "SmFe_{12}",\
     "SmFe11V_f_m0m4exp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 2 lw 1.5 lc rgb "red" title "f-site",\
     "SmFe11V_i_m0m4exp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 2 lw 1.5 lc rgb "dark-green" title "i-site",\
     "SmFe11V_j_m0m4exp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 2 lw 1.5 lc rgb "blue" title "j-site",\
     "SmFe11V_f_fullexp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "red" title "f-site",\
     "SmFe11V_i_fullexp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "dark-green" title "i-site",\
     "SmFe11V_j_fullexp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "blue" title "j-site"
unset label
set format x
set format y
set xtics 
set tmargin at screen bmgn+ysize
set bmargin at screen bmgn
set lmargin at screen lmgn
set rmargin at screen lmgn+xsize
set key reverse Left
keyx=30
keyy=1.3
set key keyx,keyy
set label "(e) SmFe_{11}Co    {/Times-Italic T}=0 K"  at lblx,lbly
#set label "[100]"  at 5,0.35
plot \
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 1 lc rgb "black" title "SmFe_{12}",\
     "SmFe11Co_f_m0m4exp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 2 lw 1.5 lc rgb "red" notitle "8{/Times-Italic f}-site",\
     "SmFe11Co_i_m0m4exp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 2 lw 1.5 lc rgb "dark-green" notitle "8{/Times-Italic i}-site",\
     "SmFe11Co_j_m0m4exp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 2 lw 1.5 lc rgb "blue" notitle "8{/Times-Italic j}-site",\
     "SmFe11Co_f_fullexp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "red" title "8{/Times-Italic f}-site",\
     "SmFe11Co_i_fullexp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "dark-green" title "8{/Times-Italic i}-site",\
     "SmFe11Co_j_fullexp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "blue" title "8{/Times-Italic j}-site"
unset label
unset key
set tmargin at screen bmgn+ysize
set bmargin at screen bmgn
set lmargin at screen lmgn+xsize
set rmargin at screen lmgn+2*xsize
set format y ""
set label "(f) SmFe_{11}Co   {/Times-Italic T}=400 K"  at lblx,lbly
#set label "[100]"  at 3.2,0.35
set xlabel "{/Times-Italic B} [T] (|| {/Times-Italic a}-axis)"  offset screen -0.11,0
plot \
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 1 lc rgb "black" \
     title "{/Times-Roman=12 SmFe_{12}}",\
     "SmFe11Co_f_m0m4exp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 2 lw 1.5 lc rgb "red"\
     title "{/Times-Roman=15 8{/Times-Italic f}-site}",\
     "SmFe11Co_i_m0m4exp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 2 lw 1.5 lc rgb "dark-green" \
     title "{/Times-Roman=15 8{/Times-Italic i}-site}",\
     "SmFe11Co_j_m0m4exp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 2 lw 1.5 lc rgb "blue" \
     title "{/Times-Roman=15 8{/Times-Italic j}-site}",\
     "SmFe11Co_f_fullexp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "red"\
     title "{/Times-Roman=15 8{/Times-Italic f}-site}",\
     "SmFe11Co_i_fullexp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "dark-green" \
     title "{/Times-Roman=15 8{/Times-Italic i}-site}",\
     "SmFe11Co_j_fullexp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "blue" \
     title "{/Times-Roman=15 8{/Times-Italic j}-site}"
unset multiplot
reset











set xrange [0:29]
set yrange [0:2.5]
set nokey
set ytics 0,1,2
set ytics font ",15"
set mxtics 2
set mytics 2
set size 0.75,0.75
set output "magcurve_SmFe11M_comp_lambda411_wK1Fe_100.eps"
set multiplot
# set left and right margin
lmgn=0.06
rmgn=0.73
bmgn=0.09
tmgn=0.74
xsize=(rmgn-lmgn)/2.0
ysize=(tmgn-bmgn)/3.0
# set label position
lblx=1
lbly=2.15
set tmargin at screen bmgn+3*ysize
set bmargin at screen bmgn+2*ysize
set lmargin at screen lmgn
set rmargin at screen lmgn+xsize
set format x ""
set label "(a) SmFe_{11}Ti   {/Times-Italic T}=0 K" font ",15" at lblx,lbly
#set label "[100]" font ",15" at 8.5,0.35
plot \
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 1 lc rgb "black" title "SmFe_{12}",\
     "SmFe11Ti_f_m0m4exp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 2 lw 3 lc rgb "red" title "f-site",\
     "SmFe11Ti_i_m0m4exp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 2 lw 3 lc rgb "dark-green" title "i-site",\
     "SmFe11Ti_j_m0m4exp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 2 lw 3 lc rgb "blue" title "j-site",\
     "SmFe11Ti_f_fullexp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "red" title "f-site",\
     "SmFe11Ti_i_fullexp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "dark-green" title "i-site",\
     "SmFe11Ti_j_fullexp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "blue" title "j-site"
unset label
set tmargin at screen bmgn+3*ysize
set bmargin at screen bmgn+2*ysize
set lmargin at screen lmgn+xsize
set rmargin at screen lmgn+2*xsize
set format x ""
set format y ""
set label "(b) SmFe_{11}Ti   {/Times-Italic T}=400 K" font ",15" at lblx,lbly
#set label "[100]" font ",15" at 6.8,0.35
plot \
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 1 lc rgb "black" title "SmFe_{12}",\
     "SmFe11Ti_f_m0m4exp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 2 lw 3 lc rgb "red" title "f-site",\
     "SmFe11Ti_i_m0m4exp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 2 lw 3 lc rgb "dark-green" title "i-site",\
     "SmFe11Ti_j_m0m4exp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 2 lw 3 lc rgb "blue" title "j-site",\
     "SmFe11Ti_f_fullexp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "red" title "f-site",\
     "SmFe11Ti_i_fullexp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "dark-green" title "i-site",\
     "SmFe11Ti_j_fullexp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "blue" title "j-site"
unset label
set ylabel "{/Symbol-Oblique m}_0{/Times-Italic M} [T]" font ",15" offset screen 0.01,0
set format y
set tmargin at screen bmgn+2*ysize
set bmargin at screen bmgn+ysize
set lmargin at screen lmgn
set rmargin at screen lmgn+xsize
set label "(c) SmFe_{11}V    {/Times-Italic T}=0 K" font ",15" at lblx,lbly
#set label "[100]" font ",15" at 10.5,0.35
plot \
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 1 lc rgb "black" title "SmFe_{12}",\
     "SmFe11V_f_m0m4exp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 2 lw 3 lc rgb "red" title "f-site",\
     "SmFe11V_i_m0m4exp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 2 lw 3 lc rgb "dark-green" title "i-site",\
     "SmFe11V_j_m0m4exp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 2 lw 3 lc rgb "blue" title "j-site",\
     "SmFe11V_f_fullexp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "red" title "f-site",\
     "SmFe11V_i_fullexp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "dark-green" title "i-site",\
     "SmFe11V_j_fullexp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "blue" title "j-site"
unset label
unset ylabel
set tmargin at screen bmgn+2*ysize
set bmargin at screen bmgn+ysize
set lmargin at screen lmgn+xsize
set rmargin at screen lmgn+2*xsize
set format x ""
set format y ""
set label "(d) SmFe_{11}V   {/Times-Italic T}=400 K" font ",15" at lblx,lbly
#set label "[100]" font ",15" at 9,0.35
plot \
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 1 lc rgb "black" title "SmFe_{12}",\
     "SmFe11V_f_m0m4exp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 2 lw 3 lc rgb "red" title "f-site",\
     "SmFe11V_i_m0m4exp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 2 lw 3 lc rgb "dark-green" title "i-site",\
     "SmFe11V_j_m0m4exp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 2 lw 3 lc rgb "blue" title "j-site",\
     "SmFe11V_f_fullexp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "red" title "f-site",\
     "SmFe11V_i_fullexp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "dark-green" title "i-site",\
     "SmFe11V_j_fullexp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "blue" title "j-site"
unset label
set format x
set format y
set xtics font ",15"
set tmargin at screen bmgn+ysize
set bmargin at screen bmgn
set lmargin at screen lmgn
set rmargin at screen lmgn+xsize
set key reverse Left
keyx=28
keyy=1.45
set key keyx,keyy
set label "(e) SmFe_{11}Co    {/Times-Italic T}=0 K" font ",15" at lblx,lbly
#set label "[100]" font ",15" at 5,0.35
plot \
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 1 lc rgb "black" title "SmFe_{12}",\
     "SmFe11Co_f_m0m4exp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 2 lw 3 lc rgb "red" notitle "8{/Times-Italic f}-site",\
     "SmFe11Co_i_m0m4exp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 2 lw 3 lc rgb "dark-green" notitle "8{/Times-Italic i}-site",\
     "SmFe11Co_j_m0m4exp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 2 lw 3 lc rgb "blue" notitle "8{/Times-Italic j}-site",\
     "SmFe11Co_f_fullexp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "red" title "8{/Times-Italic f}-site",\
     "SmFe11Co_i_fullexp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "dark-green" title "8{/Times-Italic i}-site",\
     "SmFe11Co_j_fullexp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "blue" title "8{/Times-Italic j}-site"
unset label
unset key
set tmargin at screen bmgn+ysize
set bmargin at screen bmgn
set lmargin at screen lmgn+xsize
set rmargin at screen lmgn+2*xsize
set format y ""
set label "(f) SmFe_{11}Co   {/Times-Italic T}=400 K" font ",15" at lblx,lbly
#set label "[100]" font ",15" at 3.2,0.35
set xlabel "{/Times-Italic B} [T] (|| {/Times-Italic a}-axis)" font ",15" offset screen -0.16,0
plot \
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 1 lc rgb "black" \
     title "{/Times-Roman=12 SmFe_{12}}",\
     "SmFe11Co_f_m0m4exp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 2 lw 3 lc rgb "red"\
     title "{/Times-Roman=15 8{/Times-Italic f}-site}",\
     "SmFe11Co_i_m0m4exp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 2 lw 3 lc rgb "dark-green" \
     title "{/Times-Roman=15 8{/Times-Italic i}-site}",\
     "SmFe11Co_j_m0m4exp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_mixexcf2_rd1.0.txt" u 1:2 w l lt 2 lw 3 lc rgb "blue" \
     title "{/Times-Roman=15 8{/Times-Italic j}-site}",\
     "SmFe11Co_f_fullexp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "red"\
     title "{/Times-Roman=15 8{/Times-Italic f}-site}",\
     "SmFe11Co_i_fullexp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "dark-green" \
     title "{/Times-Roman=15 8{/Times-Italic i}-site}",\
     "SmFe11Co_j_fullexp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "blue" \
     title "{/Times-Roman=15 8{/Times-Italic j}-site}"
unset multiplot
reset







set xrange [0:29]
set yrange [0:2.5]
set nokey
set ytics 0,1,2
set ytics font ",15"
set mxtics 2
set mytics 2
set size 0.75,0.75
set output "magcurve_SmFe11M_m0m4exp_lambda411_wK1Fe_100.eps"
set multiplot
# set left and right margin
lmgn=0.06
rmgn=0.73
bmgn=0.09
tmgn=0.74
xsize=(rmgn-lmgn)/2.0
ysize=(tmgn-bmgn)/3.0
# set label position
lblx=1
lbly=2.15
set tmargin at screen bmgn+3*ysize
set bmargin at screen bmgn+2*ysize
set lmargin at screen lmgn
set rmargin at screen lmgn+xsize
set format x ""
set label "(a) SmFe_{11}Ti   {/Times-Italic T}=0 K" font ",15" at lblx,lbly
#set label "[100]" font ",15" at 8.5,0.35
plot \
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 1 lc rgb "black" title "SmFe_{12}",\
     "SmFe11Ti_f_m0m4exp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 2 lw 3 lc rgb "red" title "f-site",\
     "SmFe11Ti_i_m0m4exp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 2 lw 3 lc rgb "dark-green" title "i-site",\
     "SmFe11Ti_j_m0m4exp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 2 lw 3 lc rgb "blue" title "j-site",\
     "SmFe11Ti_f_fullexp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "red" title "f-site",\
     "SmFe11Ti_i_fullexp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "dark-green" title "i-site",\
     "SmFe11Ti_j_fullexp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "blue" title "j-site"
unset label
set tmargin at screen bmgn+3*ysize
set bmargin at screen bmgn+2*ysize
set lmargin at screen lmgn+xsize
set rmargin at screen lmgn+2*xsize
set format x ""
set format y ""
set label "(b) SmFe_{11}Ti   {/Times-Italic T}=400 K" font ",15" at lblx,lbly
#set label "[100]" font ",15" at 6.8,0.35
plot \
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 1 lc rgb "black" title "SmFe_{12}",\
     "SmFe11Ti_f_m0m4exp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 2 lw 3 lc rgb "red" title "f-site",\
     "SmFe11Ti_i_m0m4exp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 2 lw 3 lc rgb "dark-green" title "i-site",\
     "SmFe11Ti_j_m0m4exp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 2 lw 3 lc rgb "blue" title "j-site",\
     "SmFe11Ti_f_fullexp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "red" title "f-site",\
     "SmFe11Ti_i_fullexp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "dark-green" title "i-site",\
     "SmFe11Ti_j_fullexp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "blue" title "j-site"
unset label
set ylabel "{/Symbol-Oblique m}_0{/Times-Italic M} [T]" font ",15" offset screen 0.01,0
set format y
set tmargin at screen bmgn+2*ysize
set bmargin at screen bmgn+ysize
set lmargin at screen lmgn
set rmargin at screen lmgn+xsize
set label "(c) SmFe_{11}V    {/Times-Italic T}=0 K" font ",15" at lblx,lbly
#set label "[100]" font ",15" at 10.5,0.35
plot \
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 1 lc rgb "black" title "SmFe_{12}",\
     "SmFe11V_f_m0m4exp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 2 lw 3 lc rgb "red" title "f-site",\
     "SmFe11V_i_m0m4exp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 2 lw 3 lc rgb "dark-green" title "i-site",\
     "SmFe11V_j_m0m4exp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 2 lw 3 lc rgb "blue" title "j-site",\
     "SmFe11V_f_fullexp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "red" title "f-site",\
     "SmFe11V_i_fullexp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "dark-green" title "i-site",\
     "SmFe11V_j_fullexp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "blue" title "j-site"
unset label
unset ylabel
set tmargin at screen bmgn+2*ysize
set bmargin at screen bmgn+ysize
set lmargin at screen lmgn+xsize
set rmargin at screen lmgn+2*xsize
set format x ""
set format y ""
set label "(d) SmFe_{11}V   {/Times-Italic T}=400 K" font ",15" at lblx,lbly
#set label "[100]" font ",15" at 9,0.35
plot \
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 1 lc rgb "black" title "SmFe_{12}",\
     "SmFe11V_f_m0m4exp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 2 lw 3 lc rgb "red" title "f-site",\
     "SmFe11V_i_m0m4exp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 2 lw 3 lc rgb "dark-green" title "i-site",\
     "SmFe11V_j_m0m4exp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 2 lw 3 lc rgb "blue" title "j-site",\
     "SmFe11V_f_fullexp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "red" title "f-site",\
     "SmFe11V_i_fullexp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "dark-green" title "i-site",\
     "SmFe11V_j_fullexp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "blue" title "j-site"
unset label
set format x
set format y
set xtics font ",15"
set tmargin at screen bmgn+ysize
set bmargin at screen bmgn
set lmargin at screen lmgn
set rmargin at screen lmgn+xsize
set key reverse Left
keyx=28
keyy=1.45
set key keyx,keyy
set label "(e) SmFe_{11}Co    {/Times-Italic T}=0 K" font ",15" at lblx,lbly
#set label "[100]" font ",15" at 5,0.35
plot \
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 1 lc rgb "black" title "SmFe_{12}",\
     "SmFe11Co_f_m0m4exp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 2 lw 3 lc rgb "red" notitle "8{/Times-Italic f}-site",\
     "SmFe11Co_i_m0m4exp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 2 lw 3 lc rgb "dark-green" notitle "8{/Times-Italic i}-site",\
     "SmFe11Co_j_m0m4exp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 2 lw 3 lc rgb "blue" notitle "8{/Times-Italic j}-site",\
     "SmFe11Co_f_fullexp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "red" title "8{/Times-Italic f}-site",\
     "SmFe11Co_i_fullexp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "dark-green" title "8{/Times-Italic i}-site",\
     "SmFe11Co_j_fullexp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "blue" title "8{/Times-Italic j}-site"
unset label
unset key
set tmargin at screen bmgn+ysize
set bmargin at screen bmgn
set lmargin at screen lmgn+xsize
set rmargin at screen lmgn+2*xsize
set format y ""
set label "(f) SmFe_{11}Co   {/Times-Italic T}=400 K" font ",15" at lblx,lbly
#set label "[100]" font ",15" at 3.2,0.35
set xlabel "{/Times-Italic B} [T]" font ",15" offset screen -0.16,0
plot \
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 1 lc rgb "black" \
     title "{/Times-Roman=12 SmFe_{12}}",\
     "SmFe11Co_f_m0m4exp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 2 lw 3 lc rgb "red"\
     title "{/Times-Roman=15 8{/Times-Italic f}-site}",\
     "SmFe11Co_i_m0m4exp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 2 lw 3 lc rgb "dark-green" \
     title "{/Times-Roman=15 8{/Times-Italic i}-site}",\
     "SmFe11Co_j_m0m4exp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 2 lw 3 lc rgb "blue" \
     title "{/Times-Roman=15 8{/Times-Italic j}-site}",\
     "SmFe11Co_f_fullexp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "red"\
     title "{/Times-Roman=15 8{/Times-Italic f}-site}",\
     "SmFe11Co_i_fullexp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "dark-green" \
     title "{/Times-Roman=15 8{/Times-Italic i}-site}",\
     "SmFe11Co_j_fullexp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "blue" \
     title "{/Times-Roman=15 8{/Times-Italic j}-site}"
unset multiplot
reset










set terminal postscript eps enhanced color font "Times-Roman"

set xrange [0:29]
set yrange [0:2.5]
set nokey
set ytics 0,1,2
set ytics font ",15"
set mxtics 2
set mytics 2
set size 0.75,0.75
set output "magcurve_SmFe11M_m0m4onlyexp_lambda411_wK1Fe_100.eps"
set multiplot
# set left and right margin
lmgn=0.06
rmgn=0.73
bmgn=0.09
tmgn=0.74
xsize=(rmgn-lmgn)/2.0
ysize=(tmgn-bmgn)/3.0
# set label position
lblx=1
lbly=2.15
set tmargin at screen bmgn+3*ysize
set bmargin at screen bmgn+2*ysize
set lmargin at screen lmgn
set rmargin at screen lmgn+xsize
set format x ""
set label "(a) SmFe_{11}Ti   {/Times-Italic T}=0 K" font ",15" at lblx,lbly
#set label "[100]" font ",15" at 8.5,0.35
plot \
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 1 lc rgb "black" title "SmFe_{12}",\
     "SmFe11Ti_f_m0m4exp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "red" title "f-site",\
     "SmFe11Ti_i_m0m4exp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "dark-green" title "i-site",\
     "SmFe11Ti_j_m0m4exp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "blue" title "j-site"
unset label
set tmargin at screen bmgn+3*ysize
set bmargin at screen bmgn+2*ysize
set lmargin at screen lmgn+xsize
set rmargin at screen lmgn+2*xsize
set format x ""
set format y ""
set label "(b) SmFe_{11}Ti   {/Times-Italic T}=400 K" font ",15" at lblx,lbly
#set label "[100]" font ",15" at 6.8,0.35
plot \
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 1 lc rgb "black" title "SmFe_{12}",\
     "SmFe11Ti_f_m0m4exp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "red" title "f-site",\
     "SmFe11Ti_i_m0m4exp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "dark-green" title "i-site",\
     "SmFe11Ti_j_m0m4exp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "blue" title "j-site"
unset label
set ylabel "{/Symbol-Oblique m}_0{/Times-Italic M} [T]" font ",15" offset screen 0.01,0
set format y
set tmargin at screen bmgn+2*ysize
set bmargin at screen bmgn+ysize
set lmargin at screen lmgn
set rmargin at screen lmgn+xsize
set label "(c) SmFe_{11}V    {/Times-Italic T}=0 K" font ",15" at lblx,lbly
#set label "[100]" font ",15" at 10.5,0.35
plot \
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 1 lc rgb "black" title "SmFe_{12}",\
     "SmFe11V_f_m0m4exp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "red" title "f-site",\
     "SmFe11V_i_m0m4exp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "dark-green" title "i-site",\
     "SmFe11V_j_m0m4exp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "blue" title "j-site"
unset label
unset ylabel
set tmargin at screen bmgn+2*ysize
set bmargin at screen bmgn+ysize
set lmargin at screen lmgn+xsize
set rmargin at screen lmgn+2*xsize
set format x ""
set format y ""
set label "(d) SmFe_{11}V   {/Times-Italic T}=400 K" font ",15" at lblx,lbly
#set label "[100]" font ",15" at 9,0.35
plot \
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 1 lc rgb "black" title "SmFe_{12}",\
     "SmFe11V_f_m0m4exp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "red" title "f-site",\
     "SmFe11V_i_m0m4exp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "dark-green" title "i-site",\
     "SmFe11V_j_m0m4exp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "blue" title "j-site"
unset label
set format x
set format y
set xtics font ",15"
set tmargin at screen bmgn+ysize
set bmargin at screen bmgn
set lmargin at screen lmgn
set rmargin at screen lmgn+xsize
set key reverse Left
keyx=28
keyy=1.45
set key keyx,keyy
set label "(e) SmFe_{11}Co    {/Times-Italic T}=0 K" font ",15" at lblx,lbly
#set label "[100]" font ",15" at 5,0.35
plot \
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 1 lc rgb "black" title "SmFe_{12}",\
     "SmFe11Co_f_m0m4exp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "red" title "8{/Times-Italic f}-site",\
     "SmFe11Co_i_m0m4exp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "dark-green" title "8{/Times-Italic i}-site",\
     "SmFe11Co_j_m0m4exp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "blue" title "8{/Times-Italic j}-site"
unset label
unset key
set tmargin at screen bmgn+ysize
set bmargin at screen bmgn
set lmargin at screen lmgn+xsize
set rmargin at screen lmgn+2*xsize
set format y ""
set label "(f) SmFe_{11}Co   {/Times-Italic T}=400 K" font ",15" at lblx,lbly
#set label "[100]" font ",15" at 3.2,0.35
set xlabel "{/Times-Italic B} [T]" font ",15" offset screen -0.16,0
plot \
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 1 lc rgb "black" \
     title "{/Times-Roman=12 SmFe_{12}}",\
     "SmFe11Co_f_m0m4exp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "red"\
     title "{/Times-Roman=15 8{/Times-Italic f}-site}",\
     "SmFe11Co_i_m0m4exp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "dark-green" \
     title "{/Times-Roman=15 8{/Times-Italic i}-site}",\
     "SmFe11Co_j_m0m4exp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "blue" \
     title "{/Times-Roman=15 8{/Times-Italic j}-site}"
unset multiplot
reset










set terminal postscript eps enhanced color font "Times-Roman"

set xrange [0:29]
set yrange [0:2.5]
set nokey
set ytics 0,1,2
set ytics font ",15"
set mxtics 2
set mytics 2
set size 0.75,0.75
set output "magcurve_SmFe11M_fullexp_lambda411_wK1Fe_mix.eps"
set multiplot
# set left and right margin
lmgn=0.06
rmgn=0.73
bmgn=0.09
tmgn=0.74
xsize=(rmgn-lmgn)/2.0
ysize=(tmgn-bmgn)/3.0
# set label position
lblx=1
lbly=2.15
set tmargin at screen bmgn+3*ysize
set bmargin at screen bmgn+2*ysize
set lmargin at screen lmgn
set rmargin at screen lmgn+xsize
set format x ""
set label "(a) SmFe_{11}Ti   {/Times-Italic T}=0 K" font ",15" at lblx,lbly
set label "[110]" font ",15" at 8.5,0.35
plot \
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T0.0_110_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 1 lc rgb "black" title "SmFe_{12}",\
     "SmFe11Ti_f_fullexp_lambda411_opencore_wK1Fe_magcurve_T0.0_110_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "red" title "f-site",\
     "SmFe11Ti_i_fullexp_lambda411_opencore_wK1Fe_magcurve_T0.0_110_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "dark-green" title "i-site",\
     "SmFe11Ti_j_fullexp_lambda411_opencore_wK1Fe_magcurve_T0.0_110_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "blue" title "j-site"
unset label
set tmargin at screen bmgn+3*ysize
set bmargin at screen bmgn+2*ysize
set lmargin at screen lmgn+xsize
set rmargin at screen lmgn+2*xsize
set format x ""
set format y ""
set label "(b) SmFe_{11}Ti   {/Times-Italic T}=400 K" font ",15" at lblx,lbly
set label "[110]" font ",15" at 6.8,0.35
plot \
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T400.0_110_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 1 lc rgb "black" title "SmFe_{12}",\
     "SmFe11Ti_f_fullexp_lambda411_opencore_wK1Fe_magcurve_T400.0_110_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "red" title "f-site",\
     "SmFe11Ti_i_fullexp_lambda411_opencore_wK1Fe_magcurve_T400.0_110_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "dark-green" title "i-site",\
     "SmFe11Ti_j_fullexp_lambda411_opencore_wK1Fe_magcurve_T400.0_110_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "blue" title "j-site"
unset label
set ylabel "{/Symbol-Oblique m}_0{/Times-Italic M} [T]" font ",15" offset screen 0.01,0
set format y
set tmargin at screen bmgn+2*ysize
set bmargin at screen bmgn+ysize
set lmargin at screen lmgn
set rmargin at screen lmgn+xsize
set label "(c) SmFe_{11}V    {/Times-Italic T}=0 K" font ",15" at lblx,lbly
set label "[110]" font ",15" at 10.5,0.35
plot \
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T0.0_110_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 1 lc rgb "black" title "SmFe_{12}",\
     "SmFe11V_f_fullexp_lambda411_opencore_wK1Fe_magcurve_T0.0_110_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "red" title "f-site",\
     "SmFe11V_i_fullexp_lambda411_opencore_wK1Fe_magcurve_T0.0_110_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "dark-green" title "i-site",\
     "SmFe11V_j_fullexp_lambda411_opencore_wK1Fe_magcurve_T0.0_110_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "blue" title "j-site"
unset label
unset ylabel
set tmargin at screen bmgn+2*ysize
set bmargin at screen bmgn+ysize
set lmargin at screen lmgn+xsize
set rmargin at screen lmgn+2*xsize
set format x ""
set format y ""
set label "(d) SmFe_{11}V   {/Times-Italic T}=400 K" font ",15" at lblx,lbly
set label "[110]" font ",15" at 9,0.35
plot \
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T400.0_110_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 1 lc rgb "black" title "SmFe_{12}",\
     "SmFe11V_f_fullexp_lambda411_opencore_wK1Fe_magcurve_T400.0_110_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "red" title "f-site",\
     "SmFe11V_i_fullexp_lambda411_opencore_wK1Fe_magcurve_T400.0_110_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "dark-green" title "i-site",\
     "SmFe11V_j_fullexp_lambda411_opencore_wK1Fe_magcurve_T400.0_110_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "blue" title "j-site"
unset label
set format x
set format y
set xtics font ",15"
set tmargin at screen bmgn+ysize
set bmargin at screen bmgn
set lmargin at screen lmgn
set rmargin at screen lmgn+xsize
set key reverse Left
keyx=28
keyy=1.45
set key keyx,keyy
set label "(e) SmFe_{11}Co    {/Times-Italic T}=0 K" font ",15" at lblx,lbly
set label "[100]" font ",15" at 5,0.35
plot \
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 1 lc rgb "black" title "SmFe_{12}",\
     "SmFe11Co_f_fullexp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "red" title "8{/Times-Italic f}-site",\
     "SmFe11Co_i_fullexp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "dark-green" title "8{/Times-Italic i}-site",\
     "SmFe11Co_j_fullexp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "blue" title "8{/Times-Italic j}-site"
unset label
unset key
set tmargin at screen bmgn+ysize
set bmargin at screen bmgn
set lmargin at screen lmgn+xsize
set rmargin at screen lmgn+2*xsize
set format y ""
set label "(f) SmFe_{11}Co   {/Times-Italic T}=400 K" font ",15" at lblx,lbly
set label "[100]" font ",15" at 3.2,0.35
set xlabel "{/Times-Italic B} [T]" font ",15" offset screen -0.16,0
plot \
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 1 lc rgb "black" \
     title "{/Times-Roman=12 SmFe_{12}}",\
     "SmFe11Co_f_fullexp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "red"\
     title "{/Times-Roman=15 8{/Times-Italic f}-site}",\
     "SmFe11Co_i_fullexp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "dark-green" \
     title "{/Times-Roman=15 8{/Times-Italic i}-site}",\
     "SmFe11Co_j_fullexp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "blue" \
     title "{/Times-Roman=15 8{/Times-Italic j}-site}"
unset multiplot
reset




set xrange [0:29]
set yrange [0:2.5]
set nokey
set ytics 0,1,2
set ytics font ",15"
set mxtics 2
set mytics 2
set size 0.75,0.75
set output "magcurve_SmFe11M_fullexp_lambda411_wK1Fe_100.eps"
set multiplot
# set left and right margin
lmgn=0.06
rmgn=0.73
bmgn=0.09
tmgn=0.74
xsize=(rmgn-lmgn)/2.0
ysize=(tmgn-bmgn)/3.0
# set label position
lblx=1
lbly=2.15
set tmargin at screen bmgn+3*ysize
set bmargin at screen bmgn+2*ysize
set lmargin at screen lmgn
set rmargin at screen lmgn+xsize
set format x ""
set label "(a) SmFe_{11}Ti   {/Times-Italic T}=0 K" font ",15" at lblx,lbly
#set label "[100]" font ",15" at 8.5,0.35
plot \
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 1 lc rgb "black" title "SmFe_{12}",\
     "SmFe11Ti_f_fullexp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "red" title "f-site",\
     "SmFe11Ti_i_fullexp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "dark-green" title "i-site",\
     "SmFe11Ti_j_fullexp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "blue" title "j-site"
unset label
set tmargin at screen bmgn+3*ysize
set bmargin at screen bmgn+2*ysize
set lmargin at screen lmgn+xsize
set rmargin at screen lmgn+2*xsize
set format x ""
set format y ""
set label "(b) SmFe_{11}Ti   {/Times-Italic T}=400 K" font ",15" at lblx,lbly
#set label "[100]" font ",15" at 6.8,0.35
plot \
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 1 lc rgb "black" title "SmFe_{12}",\
     "SmFe11Ti_f_fullexp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "red" title "f-site",\
     "SmFe11Ti_i_fullexp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "dark-green" title "i-site",\
     "SmFe11Ti_j_fullexp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "blue" title "j-site"
unset label
set ylabel "{/Symbol-Oblique m}_0{/Times-Italic M} [T]" font ",15" offset screen 0.01,0
set format y
set tmargin at screen bmgn+2*ysize
set bmargin at screen bmgn+ysize
set lmargin at screen lmgn
set rmargin at screen lmgn+xsize
set label "(c) SmFe_{11}V    {/Times-Italic T}=0 K" font ",15" at lblx,lbly
#set label "[100]" font ",15" at 10.5,0.35
plot \
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 1 lc rgb "black" title "SmFe_{12}",\
     "SmFe11V_f_fullexp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "red" title "f-site",\
     "SmFe11V_i_fullexp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "dark-green" title "i-site",\
     "SmFe11V_j_fullexp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "blue" title "j-site"
unset label
unset ylabel
set tmargin at screen bmgn+2*ysize
set bmargin at screen bmgn+ysize
set lmargin at screen lmgn+xsize
set rmargin at screen lmgn+2*xsize
set format x ""
set format y ""
set label "(d) SmFe_{11}V   {/Times-Italic T}=400 K" font ",15" at lblx,lbly
#set label "[100]" font ",15" at 9,0.35
plot \
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 1 lc rgb "black" title "SmFe_{12}",\
     "SmFe11V_f_fullexp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "red" title "f-site",\
     "SmFe11V_i_fullexp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "dark-green" title "i-site",\
     "SmFe11V_j_fullexp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "blue" title "j-site"
unset label
set format x
set format y
set xtics font ",15"
set tmargin at screen bmgn+ysize
set bmargin at screen bmgn
set lmargin at screen lmgn
set rmargin at screen lmgn+xsize
set key reverse Left
keyx=28
keyy=1.45
set key keyx,keyy
set label "(e) SmFe_{11}Co    {/Times-Italic T}=0 K" font ",15" at lblx,lbly
#set label "[100]" font ",15" at 5,0.35
plot \
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 1 lc rgb "black" title "SmFe_{12}",\
     "SmFe11Co_f_fullexp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "red" title "8{/Times-Italic f}-site",\
     "SmFe11Co_i_fullexp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "dark-green" title "8{/Times-Italic i}-site",\
     "SmFe11Co_j_fullexp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "blue" title "8{/Times-Italic j}-site"
unset label
unset key
set tmargin at screen bmgn+ysize
set bmargin at screen bmgn
set lmargin at screen lmgn+xsize
set rmargin at screen lmgn+2*xsize
set format y ""
set label "(f) SmFe_{11}Co   {/Times-Italic T}=400 K" font ",15" at lblx,lbly
#set label "[100]" font ",15" at 3.2,0.35
set xlabel "{/Times-Italic B} [T]" font ",15" offset screen -0.16,0
plot \
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 1 lc rgb "black" \
     title "{/Times-Roman=12 SmFe_{12}}",\
     "SmFe11Co_f_fullexp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "red"\
     title "{/Times-Roman=15 8{/Times-Italic f}-site}",\
     "SmFe11Co_i_fullexp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "dark-green" \
     title "{/Times-Roman=15 8{/Times-Italic i}-site}",\
     "SmFe11Co_j_fullexp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "blue" \
     title "{/Times-Roman=15 8{/Times-Italic j}-site}"
unset multiplot
reset










set terminal postscript eps enhanced color font "Times-Roman"

set xrange [0:29]
set yrange [0:2.5]
set nokey
set ytics 0,1,2
set ytics font ",15"
set mxtics 2
set mytics 2
set size 0.75,0.75
set output "magcurve_SmFe11M_m0m4exp_lambda411_wK1Fe_mix.eps"
set multiplot
# set left and right margin
lmgn=0.06
rmgn=0.73
bmgn=0.09
tmgn=0.74
xsize=(rmgn-lmgn)/2.0
ysize=(tmgn-bmgn)/3.0
# set label position
lblx=1
lbly=2.15
set tmargin at screen bmgn+3*ysize
set bmargin at screen bmgn+2*ysize
set lmargin at screen lmgn
set rmargin at screen lmgn+xsize
set format x ""
set label "(a) SmFe_{11}Ti   {/Times-Italic T}=0 K" font ",15" at lblx,lbly
#set label "[100]" font ",15" at 8.5,0.35
plot \
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 1 lc rgb "black" title "SmFe_{12}",\
     "SmFe11Ti_f_m0m4exp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 2 lw 3 lc rgb "red" title "f-site",\
     "SmFe11Ti_i_m0m4exp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 2 lw 3 lc rgb "dark-green" title "i-site",\
     "SmFe11Ti_j_m0m4exp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 2 lw 3 lc rgb "blue" title "j-site",\
     "SmFe11Ti_f_fullexp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "red" title "f-site",\
     "SmFe11Ti_i_fullexp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "dark-green" title "i-site",\
     "SmFe11Ti_j_fullexp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "blue" title "j-site"
unset label
set tmargin at screen bmgn+3*ysize
set bmargin at screen bmgn+2*ysize
set lmargin at screen lmgn+xsize
set rmargin at screen lmgn+2*xsize
set format x ""
set format y ""
set label "(b) SmFe_{11}Ti   {/Times-Italic T}=400 K" font ",15" at lblx,lbly
#set label "[100]" font ",15" at 6.8,0.35
plot \
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 1 lc rgb "black" title "SmFe_{12}",\
     "SmFe11Ti_f_m0m4exp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 2 lw 3 lc rgb "red" title "f-site",\
     "SmFe11Ti_i_m0m4exp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 2 lw 3 lc rgb "dark-green" title "i-site",\
     "SmFe11Ti_j_m0m4exp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 2 lw 3 lc rgb "blue" title "j-site",\
     "SmFe11Ti_f_fullexp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "red" title "f-site",\
     "SmFe11Ti_i_fullexp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "dark-green" title "i-site",\
     "SmFe11Ti_j_fullexp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "blue" title "j-site"
unset label
set ylabel "{/Symbol-Oblique m}_0{/Times-Italic M} [T]" font ",15" offset screen 0.01,0
set format y
set tmargin at screen bmgn+2*ysize
set bmargin at screen bmgn+ysize
set lmargin at screen lmgn
set rmargin at screen lmgn+xsize
set label "(c) SmFe_{11}V    {/Times-Italic T}=0 K" font ",15" at lblx,lbly
#set label "[100]" font ",15" at 10.5,0.35
plot \
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 1 lc rgb "black" title "SmFe_{12}",\
     "SmFe11V_f_m0m4exp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 2 lw 3 lc rgb "red" title "f-site",\
     "SmFe11V_i_m0m4exp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 2 lw 3 lc rgb "dark-green" title "i-site",\
     "SmFe11V_j_m0m4exp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 2 lw 3 lc rgb "blue" title "j-site",\
     "SmFe11V_f_fullexp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "red" title "f-site",\
     "SmFe11V_i_fullexp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "dark-green" title "i-site",\
     "SmFe11V_j_fullexp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "blue" title "j-site"
unset label
unset ylabel
set tmargin at screen bmgn+2*ysize
set bmargin at screen bmgn+ysize
set lmargin at screen lmgn+xsize
set rmargin at screen lmgn+2*xsize
set format x ""
set format y ""
set label "(d) SmFe_{11}V   {/Times-Italic T}=400 K" font ",15" at lblx,lbly
#set label "[100]" font ",15" at 9,0.35
plot \
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 1 lc rgb "black" title "SmFe_{12}",\
     "SmFe11V_f_m0m4exp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 2 lw 3 lc rgb "red" title "f-site",\
     "SmFe11V_i_m0m4exp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 2 lw 3 lc rgb "dark-green" title "i-site",\
     "SmFe11V_j_m0m4exp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 2 lw 3 lc rgb "blue" title "j-site",\
     "SmFe11V_f_fullexp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "red" title "f-site",\
     "SmFe11V_i_fullexp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "dark-green" title "i-site",\
     "SmFe11V_j_fullexp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "blue" title "j-site"
unset label
set format x
set format y
set xtics font ",15"
set tmargin at screen bmgn+ysize
set bmargin at screen bmgn
set lmargin at screen lmgn
set rmargin at screen lmgn+xsize
set key reverse Left
keyx=28
keyy=1.45
set key keyx,keyy
set label "(e) SmFe_{11}Co    {/Times-Italic T}=0 K" font ",15" at lblx,lbly
#set label "[100]" font ",15" at 5,0.35
plot \
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 1 lc rgb "black" title "SmFe_{12}",\
     "SmFe11Co_f_m0m4exp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 2 lw 3 lc rgb "red" notitle "8{/Times-Italic f}-site",\
     "SmFe11Co_i_m0m4exp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 2 lw 3 lc rgb "dark-green" notitle "8{/Times-Italic i}-site",\
     "SmFe11Co_j_m0m4exp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 2 lw 3 lc rgb "blue" notitle "8{/Times-Italic j}-site",\
     "SmFe11Co_f_fullexp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "red" title "8{/Times-Italic f}-site",\
     "SmFe11Co_i_fullexp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "dark-green" title "8{/Times-Italic i}-site",\
     "SmFe11Co_j_fullexp_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "blue" title "8{/Times-Italic j}-site"
unset label
unset key
set tmargin at screen bmgn+ysize
set bmargin at screen bmgn
set lmargin at screen lmgn+xsize
set rmargin at screen lmgn+2*xsize
set format y ""
set label "(f) SmFe_{11}Co   {/Times-Italic T}=400 K" font ",15" at lblx,lbly
#set label "[100]" font ",15" at 3.2,0.35
set xlabel "{/Times-Italic B} [T]" font ",15" offset screen -0.16,0
plot \
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 1 lc rgb "black" \
     title "{/Times-Roman=12 SmFe_{12}}",\
     "SmFe11Co_f_m0m4exp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 2 lw 3 lc rgb "red"\
     title "{/Times-Roman=15 8{/Times-Italic f}-site}",\
     "SmFe11Co_i_m0m4exp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 2 lw 3 lc rgb "dark-green" \
     title "{/Times-Roman=15 8{/Times-Italic i}-site}",\
     "SmFe11Co_j_m0m4exp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 2 lw 3 lc rgb "blue" \
     title "{/Times-Roman=15 8{/Times-Italic j}-site}",\
     "SmFe11Co_f_fullexp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "red"\
     title "{/Times-Roman=15 8{/Times-Italic f}-site}",\
     "SmFe11Co_i_fullexp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "dark-green" \
     title "{/Times-Roman=15 8{/Times-Italic i}-site}",\
     "SmFe11Co_j_fullexp_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" u 1:3 w l lt 1 lw 3 lc rgb "blue" \
     title "{/Times-Roman=15 8{/Times-Italic j}-site}"
unset multiplot
reset
















