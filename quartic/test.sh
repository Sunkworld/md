#!/bin/bash
methodlist='BAOAB ABOBA OBABO OABAO OBAB OABA BABO ABAO'
tlist='0.05'
gammalist=' 0.4 0.6 0.8 1.0 1.2 1.4 1.6 1.8 2.0' 

for method in $methodlist
do
  cd $method
  for j in unphysical
  do
  cd $j
  mkdir corfile
  for i in $gammalist
  do
    mkdir $i
    cd $i
    for t in $tlist
    do
      mkdir $t
      cd $t
      cp ../../$method-mol.f90 ../../../../random.f90 .
      sed -i 's/\:\ h\ =.*/\:\ h\ =\ '$t'd0/g' $method-mol.f90
      sed -i 's/\:\ gamma\ =.*/\:\ gamma\ =\ '$i'd0/g' $method-mol.f90
      gfortran random.f90 $method-mol.f90 -o test.x
      ./test.x
      cat result.maindat >> ../../result.maindat
      cp ../../../../autocorr.f90 ../../../../fftw3.f .
      gfortran autocorr.f90 -o autocorr.x -lfftw3
      ./autocorr.x
      echo "gamma="$i", dt="$t" Completed"
      cp mdcor ../../corfile/$t-$i-mdcor
      cd ..
    done
    cd ..
    rm -rf $i
  done
  cd ..
  done
  cd ..
done



