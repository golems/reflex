#!/usr/bin/gnuplot -persist

set term postscript enhanced eps color
set output "cp.eps"
set xlabel  "time (s)"
set title "Cart-Pendulum"

plot 'cp_x0.dat' title 'x' with lines, \
     'cp_x1.dat' title 'dx' with lines, \
     'cp_x2.dat' title '{/Symbol f}' with lines, \
     'cp_x3.dat' title 'd{/Symbol f}' with lines, \
     'cp_u.dat' title 'u' with lines, \
     'cp_xh0.dat' title 'est. x' with lines, \
     'cp_xh1.dat' title 'est. dx' with lines, \
     'cp_xh2.dat' title 'est. {/Symbol f}' with lines, \
     'cp_xh3.dat' title 'est. d{/Symbol f}' with lines
