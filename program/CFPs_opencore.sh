#!/bin/bash
filename=$(basename $(pwd))
echo "seedname is: $filename"
export filename

ifort ~/research/program/CFPs_opencore.f90
./a.out