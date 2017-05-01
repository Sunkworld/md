#!/bin/bash
mkdir stats
methodlist='BKB KBK'
cd stats
for i in $methodlist
do
  mkdir $i
  cd $i
  for j in physical unphysical
  do
    mkdir $j
  done
  cd ..
done
cd ..
for i in $methodlist
do
  cd $i
  for j in physical unphysical
  do
    cd $j
#    cp ../../stat.py .
#    ./stat.py
cd corfile
    cp ../../../Source1.f90 .
    ifort Source1.f90 -o 1.x;./1.x
cd ..
    cp corfile/*.csv result.maindat ../../stats/$i/$j
    cd ..
    echo "$i-$j completed"
  done
  cd ..
done
