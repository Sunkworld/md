#!/bin/bash
f="dt=0.60-mp-x.txt"
# write the file for gnuplot
echo "set term png enhanced

set output '${f}.png'
set xlabel 'x'
set ylabel 'p(x)'
set xrange[-3:3]
set yrange[0:0.5]
plot '${f}' u 1:2 t 'mp' smooth frequency with lines, '' u 1:3 t 'Exact' smooth frequency with lines
" |  gnuplot
