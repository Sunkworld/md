#!/bin/bash
methodlist='KBK BKB'
tlist=(0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.2 1.5 1.8)
#ASAlist=(0.9995838014 0.998340789 0.9962874679 0.9934505517 0.9898658462 0.9855768252 0.9806330145 0.9750882994 0.9689992641 0.9624236501 0.9480414979 0.9241962407 0.8987410396)
#SASlist=(1.00208 1.00829 1.0185 1.03248 1.04979 1.06973 1.09131 1.11333 1.13463 1.1544 1.19073 1.29267 1.64373)
BKBlist=(1.00042 1.00167 1.00375 1.00665 1.01032 1.01468 1.01955 1.02463 1.02947 1.03341 1.03512 1.00847 0.935453)
KBKlist=(1.00042 1.00165 1.00367 1.00642 1.00983 1.01383 1.01835 1.02333 1.02881 1.03489 1.05003 1.0901 1.17466)
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



