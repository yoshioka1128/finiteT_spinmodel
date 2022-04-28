#!/bin/bash
filename=$(basename $(pwd))
echo "seedname is: $filename"
export filename

x lapw5 -d
sed -i "s/clmval/vcoul/" lapw5.def
sed -i "s/clmvaldn/vcoul/" lapw5.def
ifort ~/research/program/arrange_rho_single.f90

cp ~/research/files/case.in5c ./
mv case.in5c "$filename".in5c
sed -i "s/0 0 1 1/1 0 0.5 1/" "$filename".in5c # (001) middle
sed -i "s/1 0 1 1/1 1 0.5 1/" "$filename".in5c
sed -i "s/0 1 1 1/0 0 0.5 1/" "$filename".in5c
lapw5 lapw5.def
./a.out
mv "$filename".rho_arrange "$filename".pot_arrange_001middle

cp ~/research/files/case.in5c ./
mv case.in5c "$filename".in5c
sed -i "s/0 0 1 1/1 0 0 1/" "$filename".in5c # (110) middle
sed -i "s/1 0 1 1/0 1 0 1/" "$filename".in5c
sed -i "s/0 1 1 1/1 0 1 1/" "$filename".in5c
lapw5 lapw5.def
./a.out
mv "$filename".rho_arrange "$filename".pot_arrange_110

cp ~/research/files/case.in5c ./
mv case.in5c "$filename".in5c
sed -i "s/0 0 1 1/0 0 0 1/" "$filename".in5c # (1-10) middle
sed -i "s/1 0 1 1/0 1 0 1/" "$filename".in5c
sed -i "s/0 1 1 1/0 0 1 1/" "$filename".in5c
lapw5 lapw5.def
./a.out
mv "$filename".rho_arrange "$filename".pot_arrange_1-10

cp ~/research/files/case.in5c ./
mv case.in5c "$filename".in5c
sed -i "s/0 0 1 1/0.5 0 0 1/" "$filename".in5c # (100) middle
sed -i "s/1 0 1 1/0.5 1 0 1/" "$filename".in5c
sed -i "s/0 1 1 1/0.5 0 1 1/" "$filename".in5c
lapw5 lapw5.def
./a.out
mv "$filename".rho_arrange "$filename".pot_arrange_100middle

cp ~/research/files/case.in5c ./
mv case.in5c "$filename".in5c
sed -i "s/0 0 1 1/1 0.5 0 1/" "$filename".in5c # (010) middle
sed -i "s/1 0 1 1/0 0.5 0 1/" "$filename".in5c
sed -i "s/0 1 1 1/1 0.5 1 1/" "$filename".in5c
lapw5 lapw5.def
./a.out
mv "$filename".rho_arrange "$filename".pot_arrange_010middle
