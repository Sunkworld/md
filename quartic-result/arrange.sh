#!/bin/bash
methodlist='ABAO ABOBA BABO BAOAB OABA OABAO OBAB OBABO BKB KBK ASA SAS'
methodlist1='ABAO ABOBA BABO BAOAB OABA OABAO OBAB OBABO'
methodlist2='BKB KBK ASA SAS'
HOME_DIR=`pwd`
mkdir stats
for m in $methodlist
do
  cd $m
  for j in physical unphysical
  do
    if [ -d $j ]
    then
      cd $j
      method=`pwd | awk -F '/' '{print $(NF-1)}'`
#      if [ "$j" = "physical" ]
#      then
        cat 0.05_pot.csv | awk '{print log($1)/2.302585092994,",",$3,$4}' > tmp1_pot.csv
        cat 0.05_ham.csv | awk '{print log($1)/2.302585092994,",",$3,$4}' > tmp1_ham.csv
#      else
#        cat 0.05_ham.csv | awk '{print $3,$4}' > tmp1_ham.csv
#        cat 0.05_pot.csv | awk '{print $3,$4}' > tmp1_pot.csv
#      fi
      cat 0.1_ham.csv | awk '{print $3,$4}' > tmp2_ham.csv
      cat 0.3_ham.csv | awk '{print $3,$4,","}' > tmp3_ham.csv
      paste -d ','  tmp1_ham.csv tmp2_ham.csv tmp3_ham.csv > ${method}-${j}-ham.csv
      cat 0.1_pot.csv | awk '{print $3,$4}' > tmp2_pot.csv
      cat 0.3_pot.csv | awk '{print $3,$4,","}' > tmp3_pot.csv
      paste -d ','  tmp1_pot.csv tmp2_pot.csv tmp3_pot.csv > ${method}-${j}-pot.csv
      sed -i '' '1i\
      log(g),0.05,0.05-error,0.1,0.1-error,0.3,0.3-error,
      ' ${method}-${j}-ham.csv
      sed -i '' '1i\
      log(g),0.05,0.05-error,0.1,0.1-error,0.3,0.3-error,
      ' ${method}-${j}-pot.csv
      rm tmp*
      mv ${method}-* ${HOME_DIR}/stats
      cd ..
    fi
  done
  cd ..
done
cd stats
for m in $methodlist1
do
  paste -d ',' $m-physical-ham.csv $m-unphysical-ham.csv > $m-ham.csv
  paste -d ',' $m-physical-pot.csv $m-unphysical-pot.csv > $m-pot.csv
done
for m in $methodlist2
do
  mv $m-physical-ham.csv $m-ham.csv
  mv $m-physical-pot.csv $m-pot.csv
done
ls | grep physical | xargs rm
