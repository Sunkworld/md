#!/bin/bash
methodlist='BAOAB ABOBA OBABO OABAO OBAB OABA BABO ABAO'
tlist='0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.2 1.5 1.8'

runmd(){
  cd $1
    for t in $tlist
    do
      mkdir $t
      cd $t
      cp ../$method-mol.f90 ../../random.f90 .
      sed -i 's/\:\ h\ =.*/\:\ h\ =\ '$t'd0/g' $method-mol.f90
      ifort random.f90 $method-mol.f90 -o test.x
      ./test.x
      cat result.maindat >> ../result.maindat
      echo "method="$method", dt="$t" Completed"
      cd ..
    done
  cd ..
}
for method in $methodlist
do
  runmd $method & 
done
wait
mkdir stats
for method in $methodlist
do
  cp $method/result.maindat stats/$method.dat
done



