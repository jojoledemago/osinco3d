#!/bin/gnuplot
reset
set view map
unset surface
set size ratio -1 
set tics out
set multiplot layout 2,2
set contour
set cntrparam levels 5
set isosamples 150
set key out
set title "ux"
splot 'outputs/solution_xy.dat' u 1:2:3 w l t ''
set title "vx"
splot 'outputs/solution_xy.dat' u 1:2:4 w l t ''
set title "uz"
splot 'outputs/solution_xy.dat' u 1:2:5 w l t ''
set title "pp"
splot 'outputs/solution_xy.dat' u 1:2:6 w l t ''
#splot 'solution_inflow.dat' u 2:1:3

