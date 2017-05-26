#!/bin/bash
methodlist=(BAOAB OBAB BABO OBABO ABOBA OABA ABAO OABAO ASA SAS BKB KBK)
for ((n=0;n<12;n++))
do
	method=${methodlist[${n}]}
	if [ ! -f $method.dat ]; then exit;fi;
	if [ $n -eq 0 ]
	then
		cat $method.dat | awk '{print $1,",",$3,",",$4,","}' >tmp${n}_p.csv
		cat $method.dat | awk  '{print $1,",",$5,",",$6,","}' >tmp${n}_k.csv
		cat $method.dat | awk  '{print $1,",",$7,",",$8,","}' >tmp${n}_t.csv
    else
        cat $method.dat | awk  '{print $3,",",$4,","}' >tmp${n}_p.csv
        cat $method.dat | awk  '{print $5,",",$6,","}' >tmp${n}_k.csv
        cat $method.dat | awk  '{print $7,",",$8,","}' >tmp${n}_t.csv
    fi
done
paste -d '' tmp0_k.csv tmp1_k.csv tmp2_k.csv tmp3_k.csv tmp4_k.csv tmp5_k.csv tmp6_k.csv tmp7_k.csv tmp8_k.csv tmp9_k.csv tmp10_k.csv tmp11_k.csv >k.csv
paste -d '' tmp0_p.csv tmp1_p.csv tmp2_p.csv tmp3_p.csv tmp4_p.csv tmp5_p.csv tmp6_p.csv tmp7_p.csv tmp8_p.csv tmp9_p.csv tmp10_p.csv tmp11_p.csv >p.csv
paste -d '' tmp0_t.csv tmp1_t.csv tmp2_t.csv tmp3_t.csv tmp4_t.csv tmp5_t.csv tmp6_t.csv tmp7_t.csv tmp8_t.csv tmp9_t.csv tmp10_t.csv tmp11_t.csv >t.csv
rm tmp*
