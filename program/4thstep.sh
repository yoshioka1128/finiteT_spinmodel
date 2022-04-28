#!/bin/bash
source get_seedname.sh
#wannier90.x "$filename"
cp "$filename""_hr.dat" fort.7
cp ~/WIEN2k/forcfp/fort.90 ./
cp ~/WIEN2k/forcfp/case.inbkq ./
~/WIEN2k/forcfp/bkq<case.inbkq>"$filename"".outbkq"
~/research/program/BkqtoAlm
cp fort.66 "$filename""_Alm_K1.txt"
cp fort.67 "$filename""_Ki.txt"
