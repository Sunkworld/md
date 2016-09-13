#!/bin/bash
f="q_mp.txt"
v="H-H"
nr=`cat $f | wc -l`
# write the file for gnuplot
echo "set term png enhanced

set output '${f}.png'
set xlabel 'x'
set ylabel 'p(x)'
set xrange[-3:3]
set yrange[0:0.5]
gauss(x)=1/sqrt(2*pi)*exp(-x*x/2)
bin_width=0.1 # set bin width    
bin_num(x)=floor(x/bin_width)  
rounded(x)=bin_width*(bin_num(x) + 0.5)
plot '${f}' u (rounded(\$1)):(1/bin_width/$nr) t 'mp(gamma=1)' smooth frequency with lines#, gauss(x) t 'Exact'
" |  gnuplot
