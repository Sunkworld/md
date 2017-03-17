#!/bin/bash
gammalist='0.01 0.02 0.05 0.1 0.2 0.4 0.6 0.8 1.0 1.2 1.4 1.6 1.8 2.0 4.0 10 20 40 50 100 200 300 1000 10000'
for i in $gammalist
do
  cd $i
  for j in 0.02 0.05 0.1
  do
    cd $j
    cp ../../autocorr.f90 ../../fftw3.f .
    gfortran autocorr.f90 -o autocorr.x -lfftw3
    ./autocorr.x
    echo "gamma="$i", dt="$j" Completed"
    cd ..
  done
  cd ..
done
