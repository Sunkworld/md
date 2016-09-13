#!/bin/bash
f="q_mp_q_1e6.txt"
g="q_mp_q_2e7.txt"
h="q_mp_q_1e8.txt"
v="H-H"
nr=`cat $f | wc -l`
nr2=`cat $g | wc -l`
nr3=`cat $h | wc -l`
# write the file for gnuplot
echo "set term png enhanced

set output '${f}.png'
set xlabel 'x'
set ylabel 'p(x)'
set xrange[-3:3]
set yrange[0:0.45]
gauss(x)=1/sqrt(2*pi)*exp(-x*x/2)
bin_width=0.1 # set bin width    
bin_num(x)=floor(x/bin_width)  
rounded(x)=bin_width*(bin_num(x) + 0.5)
plot '${f}' u (rounded(\$1)):(1/bin_width/$nr) t 'mp-1e6' smooth frequency with lines, '${g}' u (rounded(\$1)):(1/bin_width/$nr2) t 'mp-2e7' smooth frequency with lines, '${h}' u (rounded(\$1)):(1/bin_width/$nr3) t 'mp-1e8' smooth frequency with lines#, gauss(x) t 'Exact'
" |  gnuplot
