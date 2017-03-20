#!/bin/bash
gammalist='0.01 0.02 0.05 0.1 0.2 0.4 0.6 0.8 1.0 1.2 1.4 1.6 1.8 2.0 4.0 10 20 40 50 100 200 300 1000 10000'
gammalist='0.01'
for gamma in  $gammalist
do
  for t in 0.02 0.05 0.1
  do
    echo "Reading the picture for gamma=" $gamma "dt=" $t 
    xmgrace -nxy corfile/$t-$gamma-mdcor
	read -p 'Input the cutoff:' cutoff
	echo $t,$gamma`cat corfile/$t-$gamma-mdcor | grep ' '$cutoff' '` | tr ' ' ',' >> $t.csv
  done
done

