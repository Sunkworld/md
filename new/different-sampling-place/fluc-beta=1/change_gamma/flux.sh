#!/bin/bash
name="AOAB"
gammalist=(0.5 1 2 4 6 8 10 12 15)
measure=(fp fk fh)
mkdir betastats
for gamma in ${gammalist[@]}
do
	mkdir gamma=$gamma
	cd gamma=$gamma
	mkdir $name
	cp ../$name-mol.f90 $name
	cp ../random.f90 ../test.sh .
	sed -i 's/gamma\ =.*/gamma\ = '$gamma'/g' $name/$name-mol.f90
	nohup ./test.sh &
	cd ..
done
wait
for gamma in ${gammalist[@]}
do
	cp gamma=$gamma/stats/$name.dat betastats/gamma=$gamma.dat
done
cd betastats
len=${#gammalist[@]}
str=("" "" "")
for ((i=0;i<3;i++))
do
  str[${i}]="paste -d ''"
done
title=""
for ((n=0;n<len;n++))
do
	gamma=${gammalist[${n}]}
	if [ $n -eq 0 ]
	then 
		cat gamma=$gamma.dat |awk  '{print $1","$3","$4","}' >tmp${n}_fp.csv
		cat gamma=$gamma.dat |awk  '{print $1","$5","$6","}' >tmp${n}_fk.csv
		cat gamma=$gamma.dat |awk  '{print $1","$7","$8","}' >tmp${n}_fh.csv
	title+=",$gamma,,"
	else
		cat gamma=$gamma.dat |awk  '{print $3","$4","}' >tmp${n}_fp.csv
		cat gamma=$gamma.dat |awk  '{print $5","$6","}' >tmp${n}_fk.csv
		cat gamma=$gamma.dat |awk  '{print $7","$8","}' >tmp${n}_fh.csv
	title+="$gamma,,"
	fi
	for ((i=0;i<3;i++))
	do
       	  str[${i}]+=" tmp${n}_${measure[${i}]}.csv"
	done
done
for ((i=0;i<3;i++))
do
  str[${i}]+=">${measure[${i}]}.csv"
  echo ${str[${i}]} | sh
  sed -i "1i $title" ${measure[${i}]}.csv
done
rm tmp*
