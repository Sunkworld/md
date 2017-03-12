#!/bin/bash

gammalist='0.01 0.1 0.2 0.4 0.6 0.8 1.0 1.2 1.4 1.6 1.8 2.0 4.0 10 40 100 1000 10000'
for i in $gammalist
do
  mkdir $i
  cd $i
  mkdir 0.02 0.05 0.1
  for t in 0.02 0.05 0.1
  do
    cd $t
    cp ../../mol.phys.13.f90 ../../random.f90 .
    sed -i 's/\:\ h\ =\ 0.1d0/\:\ h\ =\ '$t'd0/g' mol.phys.13.f90
    sed -i 's/\:\ gamma\ =\ 0.8d0/\:\ gamma\ =\ '$i'd0/g' mol.phys.13.f90
    gfortran random.f90 mol.phys.13.f90 -o test.x
    ./test.x
    cat result.maindat >> ../../result.maindat
    cd ..
  done
  cd ..
done
