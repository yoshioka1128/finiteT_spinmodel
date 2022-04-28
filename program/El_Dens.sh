#!/bin/bash
filename=$(basename $(pwd))
echo "seedname is: $filename"
export filename

# 001
x lapw5 -d
sed -i "s/clmval'/clmsum'/" lapw5.def
cp ~/research/files/case.in5c ./
mv case.in5c "$filename".in5c
sed -i "s/0 0 1 1/0 0 0.5 1/" "$filename".in5c
sed -i "s/1 0 1 1/1 0 0.5 1/" "$filename".in5c
sed -i "s/0 1 1 1/0 1 0.5 1/" "$filename".in5c
lapw5 lapw5.def
mv "$filename".rho "$filename".rhosum_001

x lapw5 -d
sed -i "s/clmval'/clmup'/" lapw5.def
cp ~/research/files/case.in5c ./
mv case.in5c "$filename".in5c
sed -i "s/0 0 1 1/0 0 0.5 1/" "$filename".in5c
sed -i "s/1 0 1 1/1 0 0.5 1/" "$filename".in5c
sed -i "s/0 1 1 1/0 1 0.5 1/" "$filename".in5c
lapw5 lapw5.def
mv "$filename".rho "$filename".rhoup_001

x lapw5 -d
sed -i "s/clmval'/clmdn'/" lapw5.def
cp ~/research/files/case.in5c ./
mv case.in5c "$filename".in5c
sed -i "s/0 0 1 1/0 0 0.5 1/" "$filename".in5c
sed -i "s/1 0 1 1/1 0 0.5 1/" "$filename".in5c
sed -i "s/0 1 1 1/0 1 0.5 1/" "$filename".in5c
lapw5 lapw5.def
mv "$filename".rho "$filename".rhodn_001

# 110
x lapw5 -d
sed -i "s/clmval'/clmsum'/" lapw5.def
cp ~/research/files/case.in5c ./
mv case.in5c "$filename".in5c
sed -i "s/0 0 1 1/0 1 0 1/" "$filename".in5c
sed -i "s/1 0 1 1/1 0 0 1/" "$filename".in5c
sed -i "s/0 1 1 1/0 1 1 1/" "$filename".in5c
lapw5 lapw5.def
mv "$filename".rho "$filename".rhosum_110

x lapw5 -d
sed -i "s/clmval'/clmup'/" lapw5.def
cp ~/research/files/case.in5c ./
mv case.in5c "$filename".in5c
sed -i "s/0 0 1 1/0 1 0 1/" "$filename".in5c
sed -i "s/1 0 1 1/1 0 0 1/" "$filename".in5c
sed -i "s/0 1 1 1/0 1 1 1/" "$filename".in5c
lapw5 lapw5.def
mv "$filename".rho "$filename".rhoup_110

x lapw5 -d
sed -i "s/clmval'/clmdn'/" lapw5.def
cp ~/research/files/case.in5c ./
mv case.in5c "$filename".in5c
sed -i "s/0 0 1 1/0 1 0 1/" "$filename".in5c
sed -i "s/1 0 1 1/1 0 0 1/" "$filename".in5c
sed -i "s/0 1 1 1/0 1 1 1/" "$filename".in5c
lapw5 lapw5.def
mv "$filename".rho "$filename".rhodn_110

# 100
x lapw5 -d
sed -i "s/clmval'/clmsum'/" lapw5.def
cp ~/research/files/case.in5c ./
mv case.in5c "$filename".in5c
sed -i "s/0 0 1 1/0.5 0 0 1/" "$filename".in5c
sed -i "s/1 0 1 1/0.5 1 0 1/" "$filename".in5c
sed -i "s/0 1 1 1/0.5 0 1 1/" "$filename".in5c
lapw5 lapw5.def
mv "$filename".rho "$filename".rhosum_100

x lapw5 -d
sed -i "s/clmval'/clmup'/" lapw5.def
cp ~/research/files/case.in5c ./
mv case.in5c "$filename".in5c
sed -i "s/0 0 1 1/0.5 0 0 1/" "$filename".in5c
sed -i "s/1 0 1 1/0.5 1 0 1/" "$filename".in5c
sed -i "s/0 1 1 1/0.5 0 1 1/" "$filename".in5c
lapw5 lapw5.def
mv "$filename".rho "$filename".rhoup_100

x lapw5 -d
sed -i "s/clmval'/clmdn'/" lapw5.def
cp ~/research/files/case.in5c ./
mv case.in5c "$filename".in5c
sed -i "s/0 0 1 1/0.5 0 0 1/" "$filename".in5c
sed -i "s/1 0 1 1/0.5 1 0 1/" "$filename".in5c
sed -i "s/0 1 1 1/0.5 0 1 1/" "$filename".in5c
lapw5 lapw5.def
mv "$filename".rho "$filename".rhodn_100


ifort ~/research/program/arrange_rho.f90
./a.out
