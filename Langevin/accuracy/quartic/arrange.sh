#!/bin/bash
#methodlist='ABAO ABOBA BABO BAOAB OABA OABAO OBAB OBABO BKB KBK ASA SAS'
methodlist=(BAOAB OBAB BABO OBABO ABOBA OABA ABAO OABAO ASA SAS BKB KBK)
methodlist1='ABAO ABOBA BABO BAOAB OABA OABAO OBAB OBABO'
methodlist2='BKB KBK ASA SAS'
HOME_DIR=`pwd`
mkdir stats

for i in fixed_pot 
do
  cd $i
  for j in unphysical
  do
    cd $j/stats
    for ((n=0;n<12;n++))
    do
      method=${methodlist[${n}]}
#      if [ $n -eq 0 ]
#      then
        cat $method.csv | awk -F ',' '{print $1,",",$3,",",$4,","}' >tmp${n}_p.csv
        cat $method.csv | awk -F ',' '{print $1,",",$5,",",$6,","}' >tmp${n}_k.csv
#      else
#        cat $method.csv | awk -F ',' '{print $3,",",$4,","}' >tmp${n}_p.csv
#        cat $method.csv | awk -F ',' '{print $5,",",$6,","}' >tmp${n}_k.csv
#      fi
    done
    cp ${HOME_DIR}/exact_*.csv .
    paste tmp0_k.csv tmp1_k.csv tmp2_k.csv tmp3_k.csv tmp4_k.csv tmp5_k.csv tmp6_k.csv tmp7_k.csv tmp8_k.csv tmp9_k.csv tmp10_k.csv tmp11_k.csv exact_k.csv>$i-$j-tmp-k.csv
    paste tmp0_p.csv tmp1_p.csv tmp2_p.csv tmp3_p.csv tmp4_p.csv tmp5_p.csv tmp6_p.csv tmp7_p.csv tmp8_p.csv tmp9_p.csv tmp10_p.csv tmp11_p.csv exact_p.csv>$i-$j-tmp-p.csv
    if [ "$j" = "physical" ]
    then
      sed -i '' '1i\
        dt,middle,middle-error,end,end-error,beginning,beginning-error,side,side-error,PV-middle,PV-middle-error,PV-end,PV-end-error,PV-beginning,PV-beginning-error,PV-side,PV-side-error,middle-pT,middle-pT-error,side-pT,side-pT-error,middle-xT,middle-xT-error,side-xT,side-xT-error,dt,exact
      '  $i-$j-tmp-k.csv
      sed -i '' '1i\
        dt,middle,middle-error,end,end-error,beginning,beginning-error,side,side-error,PV-middle,PV-middle-error,PV-end,PV-end-error,PV-beginning,PV-beginning-error,PV-side,PV-side-error,middle-pT,middle-pT-error,side-pT,side-pT-error,middle-xT,middle-xT-error,side-xT,side-xT-error,dt,exact
      '  $i-$j-tmp-p.csv

      cat $i-$j-tmp-k.csv | awk -F ',' '{
      if (NF>2) {print $0}
      else {print ",,,,,,,,,,,,,,,,,,,,,,,,,",$0}}'>$i-$j-k.csv
      cat $i-$j-tmp-p.csv | awk -F ',' '{
      if (NF>2) {print $0}
      else {print ",,,,,,,,,,,,,,,,,,,,,,,,,",$0}}'>$i-$j-p.csv
      rm $i-$j-tmp-k.csv $i-$j-tmp-p.csv
    else
      sed -i '' '1i\
        dt,middle,middle-error,dt,end,end-error,dt,beginning,beginning-error,dt,side,side-error,dt,PV-middle,PV-middle-error,dt,PV-end,PV-end-error,dt,PV-beginning,PV-beginning-error,dt,PV-side,PV-side-error,dt,exact
      '  $i-$j-tmp-k.csv
      sed -i '' '1i\
        dt,middle,middle-error,dt,end,end-error,dt,beginning,beginning-error,dt,side,side-error,dt,PV-middle,PV-middle-error,dt,PV-end,PV-end-error,dt,PV-beginning,PV-beginning-error,dt,PV-side,PV-side-error,dt,exact
      '  $i-$j-tmp-p.csv

      cat $i-$j-tmp-k.csv | awk -F ',' '{
      if (NF>2) {print $0}
      else {print ",,,,,,,,,,,,,,,,,,,,,,,,",$0}}'>$i-$j-k.csv
      cat $i-$j-tmp-p.csv | awk -F ',' '{
      if (NF>2) {print $0}
      else {print ",,,,,,,,,,,,,,,,,,,,,,,,",$0}}'>$i-$j-p.csv
      rm $i-$j-tmp-k.csv $i-$j-tmp-p.csv

    fi
    mv $i-$j-* ${HOME_DIR}/stats
    rm tmp*
    cd ../../
  done
  cd ..
done
