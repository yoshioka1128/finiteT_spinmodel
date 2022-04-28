set terminal postscript eps enhanced color font "Times-Roman"
set size 0.5,0.5
set lmargin at screen 0.07
set key reverse Left
set output "magcurve_annum_100_T0_T400.eps"
set key at 15.5,0.8
set xrange [0:15]
set yrange [0:2]
set mxtics 2
set mytics 2
set ytics 0,0.5,2.5
set xlabel "{/Times-Italic B} [T] (|| {/Times-Italic a}-axis)" 
set ylabel "{/Symbol}{/Symbol-Oblique m}_0 {/Times-BoldItalic M}_{/Times-Roman s}({/Times-Italic T})\
 {/Symbol \327} {/Times-BoldItalic n}_{/Times-Italic a} [T]" 2,0
HN=13.8127238275133
K1=10.0249475838234
K2=-3.37717119734796     
Ms=1.82407479221322 
HN400=5.19963202289366
Ms400=1.30342073103541 
mu0=4.0*pi*0.1

#set arrow from 13.8,1.56 to 13.8,1.75 head nofilled size screen 0.01,15 lw 0.5
#set arrow from 5.2,1.58 to 5.2,1.4 head nofilled size screen 0.01,15 lw 0.5
#set label "{/Times-Italic B}_{/Times-Roman=8 N}(0)" font ",11" tc rgb "blue" at 13.5,1.7
#set label "{/Times-Italic B}_{/Times-Roman=8 FP}(0)" font ",11" tc rgb "blue" at 4.2,1.8
#set label "{/Times-Italic B}_{/Times-Roman=8 N}(400)" font ",11" tc rgb "red" at 4.2,1.6

#set label "{/Times-Italic B}_{/Times-Roman=8 N}(0)=13.8 T" font ",11" tc rgb "blue" at 11.5,1.47
#set label "{/Times-Italic B}_{/Times-Roman=8 N}(400)=5.2 T" font ",11" tc rgb "red" at 2.5,1.65

set label "    {/Times-Italic T}=0 K" tc rgb "blue" at 6.5,1.73
set label "{/Times-Italic T}=400 K"   tc rgb "red" at 12,1.22

set label "broken curves: statistics" font ",12"   at 7.95,0.4
set label "solid curves: model C"  font ",12"  at 8.5,0.25

ax=(-1+sqrt(-3*K1/K2-2))/3.0
a=(2.0*K1*ax+4.0*K2*ax**3)*mu0/Ms
#print a,ax*Ms

set arrow from 0,0 to HN,Ms nohead lt 5 lc rgb "black"
set arrow from 0,0 to HN400,Ms400 nohead lt 5 lc rgb "black"

set key at 14,0.6
set nokey
plot \
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" \
     u 1:3 w l lw 3 lt 2 lc rgb "blue" notitle,\
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T0.0_100_mixexcf2_rd1.0.txt" \
     w l lw 3 lt 1 lc rgb "blue" notitle,\
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" \
     u 1:3 w l lw 3 lt 2 lc rgb "red" notitle,\
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T400.0_100_mixexcf2_rd1.0.txt" \
     w l lw 3 lt 1 lc rgb "red" notitle,\
     "point2.txt" w p lw 2 pt 6 ps 1 lc rgb "black" title "{/Times-Roman=12 {/Times-Italic B}_N({/Times-Italic T})}"
reset




set size 0.5,0.5
set lmargin at screen 0.07
set key reverse Left
set output "magcurve_annum_K1K2model.eps"
set key at 15.5,0.8
set xrange [0:15]
set yrange [0:2]
set mxtics 2
set mytics 2
set ytics 0,0.5,2.5
set xlabel "Applied Field [T]" 
set ylabel "{/Symbol}{/Symbol-Oblique m}_0{/Times-Italic M_{{/Times-Roman s},x}} [T]" 2,0
#set ylabel "{/Symbol}{/Symbol m}{/Times-Italic M_{s,x}} [T]" 2,0
HN=13.8127238275133
K10=10.0249475838234
K20=-3.37717119734796     
Ms0=1.82407479221322 
HN400=5.19963202289366
Ms400=1.30342073103541 
mu0=4.0*pi*0.1

p=K10/(2.0*K20)/Ms0
q=-Ms0/(4.0*K20)

MT0(x)=(-q*x/2+sqrt(((q*x)/2.0)**2.0+(p/3.0)**3.0))**0.333333

x=0.1
print ((q*x)/2.0)**2,(p/3.0)**3

set arrow from 13.8,1.56 to 13.8,1.75 head nofilled size screen 0.01,15 lw 0.5
set arrow from 5.2,1.58 to 5.2,1.4 head nofilled size screen 0.01,15 lw 0.5
#set label "{/Times-Italic=12 y={/Times-Roman [}M{/Times-Roman (}T{/Times-Roman )}/H_{/Times-Roman=8 N}{/Times-Roman=12 ]}x" at 10.2,1.3
set label "{/Times-Italic=12 H_{/Times-Roman=8 N}{/Times-Roman=12 =13.8 T}" at 11.7,1.46
set label "{/Times-Italic=12 H_{/Times-Roman=8 N}{/Times-Roman=12 =5.2 T}" at 3,1.65

set label "    {/Times-Italic T}=0 K" tc rgb "blue" at 3.2,1.87
set label "{/Times-Italic T}=400 K"   tc rgb "red" at 2,1.4

set label "broken lines: numerical"   at 6,0.5
set label "solid lines: model C"   at 6.6,0.35

ax=(-1+sqrt(-3*K10/K20-2))/3
a=(2*K10*ax+4*K20*ax**3)*mu0/Ms0

set nokey
plot \
     MT0(x) w l lw 3 lc rgb "purple" ,\
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" \
     u 1:3 w l lw 3 lt 2 lc rgb "blue" title "{/Times-Roman=12 numerical}",\
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T0.0_100_mixexcf2_rd1.0.txt" \
     w l lw 3 lt 1 lc rgb "blue" title "{/Times-Roman=12 model C}",\
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T400.0_100_nJex2_rd1.0.txt" \
     u 1:3 w l lw 3 lt 2 lc rgb "red" title "{/Times-Roman=12 numerical}",\
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T400.0_100_mixexcf2_rd1.0.txt" \
     w l lw 3 lt 1 lc rgb "red" title "{/Times-Roman=12 model C}",\
     x/HN*Ms w l lw 2 lt 5 lc rgb "black" notitle,\
     x/HN400*Ms400 w l lw 2 lt 5 lc rgb "black" notitle,\
     "point2.txt" w p lw 2 pt 6 ps 2 lc rgb "black" notitle,\
     "point3.txt" w p lw 2 pt 2 ps 1 lc rgb "black" notitle
reset





set size 0.5,0.5
set lmargin at screen 0.07
set key reverse Left
set output "magcurve_annum_100_T0.eps"
set key at 15.5,0.8
set xrange [0:15]
set yrange [0:2]
set mxtics 2
set mytics 2
set ytics 0,0.5,2.5
set xlabel "Applied Field [T]" 
set ylabel "Magnetization [T]" 2,0
HN=13.8127238275133
Ms=1.82407479221322 
set arrow from 13.8,1.56 to 13.8,1.75 head nofilled size screen 0.01,15 lw 0.5

set label "{/Times-Italic=12 y={/Times-Roman [}M{/Times-Roman (}T{/Times-Roman )}/H_{/Times-Roman=8 N}{/Times-Roman=12 ]}x" at 10.2,1.3
set label "{/Times-Italic=12 H_{/Times-Roman=8 N}{/Times-Roman=12 =13.8 T}" at 11.7,1.46
plot \
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T0.0_100_nJex2_rd1.0.txt" \
     u 1:3 w l lw 3 lt 2 lc rgb "red" title "{/Times-Roman=12 numerical}",\
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T0.0_100_lowestJ_rd1.0.txt" \
     w l lw 1 lt 1 lc rgb "dark-green" title "{/Times-Roman=12 model A}",\
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T0.0_100_mixexcf1_rd1.0.txt" \
     w l lw 1 lt 1 lc rgb "blue" title "{/Times-Roman=12 model B}",\
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T0.0_100_mixexcf2_rd1.0.txt" \
     w l lw 3 lt 1 lc rgb "red" title "{/Times-Roman=12 model C}",\
     x/HN*Ms w l lw 1 lt 5 lc rgb "black" notitle,\
     "point2.txt" w p pt 2 ps 1 lc rgb "black" notitle
reset





set size 0.5,0.5
set lmargin at screen 0.07
set key reverse Left
set output "magcurve_annum_100_T300.eps"
set key at 15.5,0.8
set xrange [0:15]
set yrange [0:2]
set mxtics 2
set mytics 2
set ytics 0,0.5,2.5
set xlabel "Applied Field [T]" 
set ylabel "Magnetization [T]" 2,0
HN=13.8127238275133
Ms=1.82407479221322 
set arrow from 13.8,1.56 to 13.8,1.75 head nofilled size screen 0.01,15 lw 0.5
set label "{/Times-Italic=12 y={/Times-Roman [}M{/Times-Roman (}T{/Times-Roman )}/H_{/Times-Roman=8 N}{/Times-Roman=12 ]}x" at 10.2,1.3
set label "{/Times-Italic=12 H_{/Times-Roman=8 N}{/Times-Roman=12 =13.8 T}" at 11.7,1.46
plot \
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T300.0_100_nJex2_rd1.0.txt" u 1:3 w l lw 2 lt 2 lc rgb "red" title "{/Times-Roman=12 numerical}",\
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T300.0_100_lowestJ_rd1.0.txt" w l lt 1 lc rgb "dark-green" title "{/Times-Roman=12 model A}",\
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T300.0_100_mixexcf1_rd1.0.txt" w l lt 1 lc rgb "blue" title "{/Times-Roman=12 model B}",\
     "SmFe12_lambda411_opencore_wK1Fe_magcurve_T300.0_100_mixexcf2_rd1.0.txt" w l lw 2 lt 1 lc rgb "red" title "{/Times-Roman=12 model C}",\
     x/HN*Ms w l lw 1 lt 5 lc rgb "black" notitle,\
     "point2.txt" w p pt 2 ps 1 lc rgb "black" notitle
reset

