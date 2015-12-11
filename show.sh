#!/bin/sh

./run
gnuplot -p -e "plot 'points.txt', 'line.txt' w l, sin(2*pi*x)"
gnuplot -p -e "set logscale y; set format y '10^{%L}'; plot 'error.txt' u 1:2 w l, 'error.txt' u 1:3 w l"
#gnuplot -p -e "set logscale y; set format y '10^{%L}'; plot 'errdep.txt' w l"
