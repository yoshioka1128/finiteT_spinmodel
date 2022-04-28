set terminal postscript eps enhanced color font "Times-Roman,16"

set polar
set yrange [-7:7]
set xrange [-7:7]
set xtics 2
set mxtics 

set angle degree
set size square
set grid polar 45
set grid lt 1 lc rgb "black"
#set format x ""
#set format y ""

set style arrow 1 head filled size 0.7,30,70 lw 2	lc rgb "dark-violet"
set style arrow 2 head filled size 0.7,30,70 lw 2	lc rgb "blue"
set style arrow 3 head filled size 0.7,30,70 lw 2	lc rgb "red"
set style arrow 4 nohead filled size 0.7,30,70 lt 2 lw 2 	lc rgb "black"

set output "magstruct.eps"
set size 0.685,0.5
set noborder
set xtics scale 0,0
set ytics scale 0,0
set format x ""
set format y ""
set multiplot
set lmargin at screen 0.0
set rmargin at screen 0.54
set bmargin at screen 0.01
set tmargin at screen 0.41
set nokey

set label "{/Times-Bold-Italic m}@_{{/Times-Italic L},1}" at 3.9,1.3 tc rgb "red"
set label "{/Times-Bold-Italic m}@_{{/Times-Italic L},2}" at 1.5,4.4 tc rgb "red"
set label "{/Times-Bold-Italic m}@_{{/Times-Italic S},1}" at -5.7,-1.3 tc rgb "blue"
set label "{/Times-Bold-Italic m}@_{{/Times-Italic S},2}" at -3.2,-3.9 tc rgb "blue"
set label "10{/Times-Bold-Italic m}@_{1}" at 1.9,-1.2 tc rgb "dark-violet"
set label "10{/Times-Bold-Italic m}@_{2}" at -2.4,2.7 tc rgb "dark-violet"
set label "{/Times-Bold-Italic M}^{TM}" at 6.0,4.9 tc rgb "black"
set label "{/Times-Italic a}-axis" at 6.2,0.1
set label "{/Times-Italic b}-axis" at -1.2,6.7
set label "2.0" at 2.0/sqrt(2.0)-0.8,-2.0/sqrt(2.0)-0.9
set label "4.0" at 4.0/sqrt(2.0)-0.8,-4.0/sqrt(2.0)-0.9
set label "6.0" at 6.0/sqrt(2.0)-0.8,-6.0/sqrt(2.0)-0.9
set label "(a) {/Times-Italic B}=1 T (|| {/Times-Italic a}-axis)" font ",18" at -7,8.5

plot \
"magstruct_SmFe11V_i_S.txt"  every ::0::0 u ($1-1.0):($1-1.0):($3):($4) w vector arrowstyle 2,\
"magstruct_SmFe11V_i_L.txt"  every ::0::0 u ($1-1.0):($1-1.0):($3):($4) w vector arrowstyle 3,\
"magstruct_SmFe11V_i_J.txt"  every ::0::0 u ($1-1.0):($1-1.0):($3):(10*$4) w vector arrowstyle 1,\
"magstruct_SmFe11V_i_S.txt"  every ::0::0 u ($1-1.0):($1-1.0):($6):($7) w vector arrowstyle 2,\
"magstruct_SmFe11V_i_L.txt"  every ::0::0 u ($1-1.0):($1-1.0):($6):($7) w vector arrowstyle 3,\
"magstruct_SmFe11V_i_J.txt"  every ::0::0 u ($1-1.0):($1-1.0):($6):(10*$7) w vector arrowstyle 1,\
"magstruct_SmFe11V_i_Fe.txt" every ::0::0 u ($1-1.0):($1-1.0):($3):($4*7.5/41.64459) w vector arrowstyle 4

unset label
unset ylabel
set lmargin at screen 0.35
set rmargin at screen 0.89
set bmargin at screen 0.01
set tmargin at screen 0.41

set nokey

set label "{/Times-Bold-Italic m}@_{{/Times-Italic L},1}" at 4.05,-0.35 tc rgb "red"
set label "{/Times-Bold-Italic m}@_{{/Times-Italic L},2}" at 3.3,2.6 tc rgb "red"
set label "{/Times-Bold-Italic m}@_{{/Times-Italic S},1}" at -5.9,0.7 tc rgb "blue"
set label "{/Times-Bold-Italic m}@_{{/Times-Italic S},2}" at -5.4,-2 tc rgb "blue"
set label "10{/Times-Bold-Italic m}@_{1}" at 1.8,-1.3 tc rgb "dark-violet"
set label "10{/Times-Bold-Italic m}@_{2}" at 0,3 tc rgb "dark-violet"
set label "{/Times-Bold-Italic M}^{TM}" at 7.2,2.3 tc rgb "black"
set label "{/Times-Italic a}-axis" at 6.2,0.1
set label "{/Times-Italic b}-axis" at -1.2,6.7
set label "2.0" at 2.0/sqrt(2.0)-0.8,-2.0/sqrt(2.0)-0.9
set label "4.0" at 4.0/sqrt(2.0)-0.8,-4.0/sqrt(2.0)-0.9
set label "6.0" at 6.0/sqrt(2.0)-0.8,-6.0/sqrt(2.0)-0.9
set label "(b) {/Times-Italic B}=4 T (|| {/Times-Italic a}-axis)" font ",18" at -7,8.5
set xlabel "Magnetic Moments [{/Symbol-Oblique m}_B]" -11.6,1.8
plot \
"magstruct_SmFe11V_i_S.txt"  every ::3::3 u ($1-4.0):($1-4.0):($3):($4) w vector arrowstyle 2,\
"magstruct_SmFe11V_i_L.txt"  every ::3::3 u ($1-4.0):($1-4.0):($3):($4) w vector arrowstyle 3,\
"magstruct_SmFe11V_i_J.txt"  every ::3::3 u ($1-4.0):($1-4.0):($3):(10*$4) w vector arrowstyle 1,\
"magstruct_SmFe11V_i_S.txt"  every ::3::3 u ($1-4.0):($1-4.0):($6):($7) w vector arrowstyle 2,\
"magstruct_SmFe11V_i_L.txt"  every ::3::3 u ($1-4.0):($1-4.0):($6):($7) w vector arrowstyle 3,\
"magstruct_SmFe11V_i_J.txt"  every ::3::3 u ($1-4.0):($1-4.0):($6):(10*$7) w vector arrowstyle 1,\
"magstruct_SmFe11V_i_Fe.txt" every ::3::3 u ($1-4.0):($1-4.0):($3):($4) w vector arrowstyle 4
unset multiplot

reset