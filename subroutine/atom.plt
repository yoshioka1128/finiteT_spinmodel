# shape of unit cell
  ix=8
  iy=1
  iz=3
# enlargement ratio
  enlarge=3.0
# resolution
  resolution1=1600
  resolution2=750
# unit time
  dtime=0.4
# time setting
  nini=1 # initial
  nend=200 # final
  dt=1 # time step

set view equal xyz
set view 80,15,enlarge,1

ax=8.81
cx=12.21
set xrange [-2:ix*ax+2]
set yrange [-2:iy*ax+2]
set zrange [-2:iz*cx+2]

set ticslevel 0
set nokey

# initialize
ix0=0
iy0=0
iz0=0

load "cellz_9_1_1.plt"

unset border
set xtics -100,1,-99
set ytics -100,1,-99
set ztics -100,1,-99
unset colorbox
unset cblabel
unset key

set style arrow 1 head filled size 5,2,5 lw 2 lc rgb "red"
set style arrow 2 head filled size 5,2,5 lw 2 lc rgb "blue"
set style arrow 3 head filled size 5,2,5 lw 2 lc rgb "grey40"



# magnetization
set terminal gif animate optimize size resolution1, resolution2
set output "./mag_atom.gif"
set tics font 'Times,18'
set origin 0,0
load "mag_atom.plt"
n = 0

# anisotropy eng
#set terminal gif animate optimize size 2000, 480
#set cbrange [0:20]
#set palette define (0 "blue",20 "red") #HA rm Nd and surface Fe
#set output "./100gb_EA.gif"
#set tics font 'Times,18'
#set origin -2.9,-2.9
#load "EA_9_1_1.plt"
#n=0
#
## exchange eng
#set terminal gif animate optimize size 2000, 480
#set cbrange [0:20]
#set palette define (0 "blue",20 "red") #HA rm Nd and surface Fe
#set output "./100gb_Eex.gif"
#set tics font 'Times,18'
#set origin -2.9,-2.9
#load "Eex_9_1_1.plt"
#n=0
#
## Zeeman eng
#set terminal gif animate optimize size 2000, 480
#set cbrange [0:20]
#set palette define (0 "blue",20 "red") #HA rm Nd and surface Fe
#set output "./100gb_EH.gif"
#set tics font 'Times,18'
#set origin -2.9,-2.9
#load "EH_9_1_1.plt"
#n=0
#
## eng sum
#set terminal gif animate optimize size 2000, 480
#set cbrange [0:20]
#set palette define (0 "blue",20 "red") #HA rm Nd and surface Fe
#set output "./100gb_Esum.gif"
#set tics font 'Times,18'
#set origin -2.9,-2.9
#load "Esum_9_1_1.plt"
#n=0

reset
