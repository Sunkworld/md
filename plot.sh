#!/bin/bash
f="baoab.txt"
g="mp.txt"
# write the file for gnuplot
echo "set term png enhanced

set output '${f}.png'
set xlabel 'dt'
set ylabel 'E_p'

f(x)=0.404715
plot '${f}' u 1:3:4 t 'baoab' with yerrorlines, '${g}' u 1:3:4 t 'mp' with yerrorlines, f(x)
 " |  gnuplot
