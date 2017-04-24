#!/bin/bash
methodlist='ASA'
tlist=(0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.2 1.5 1.8)
ASAlist=(1.98699 1.95134 1.90078 1.84299 1.78342 1.72524 1.67001 1.61842 1.57069 1.52695 1.45267 1.405 1.65024)
SASlist=(1.99189 1.96999 1.93986 1.90712 1.87552 1.84667 1.82046 1.79577 1.77119 1.74566 1.69191 1.63931 1.76332)
BKBlist=(1.99508 1.98125 1.96093 1.93728 1.91342 1.89186 1.87422 1.86081 1.85036 1.83981 1.80108 1.66809 1.64533)
KBKlist=(1.99999 1.99983 1.99919 1.99771 1.99523 1.99188 1.98828 1.98545 1.98472 1.98765 2.00984 2.08134 2.15428)
runmd(){
  cd $1
    for ((i=0;i<13;i++))
    do
      mkdir ${tlist[i]}
      cd ${tlist[i]}
      cp ../$method-mol.f90 ../../random.f90 .
      sed -i 's/\:\ h\ =.*/\:\ h\ =\ '${tlist[i]}'d0/g' $method-mol.f90
	  eval gamma=\${${method}list[i]}
	  sed -i 's/\:\ gamma\ =.*/\:\ gamma\ ='$gamma'd0/g' $method-mol.f90
      ifort random.f90 $method-mol.f90 -mkl -o test.x
      ./test.x
      cat result.maindat >> ../result.maindat
      echo "method="$method", dt="${tlist[i]}" Completed"
      cd ..
    done
  cd ..
}
for method in ${methodlist[@]}
do
  runmd $method & 
done
wait
mkdir stats
for method in ${methodlist[@]}
do
  cp $method/result.maindat stats/$method.dat
done



