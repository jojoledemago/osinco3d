#!/bin/gnuplot
unset multiplot
reset
set view map
unset surface
set size ratio -1 
set tics out
set contour
set xrange[0:pi]
set yrange[0:pi]
#set cntrparam levels incremental -10.6, 0.05, -10
#set cntrparam levels 6
set isosamples 150
set key out
#splot 'solution_xy.dat' u 1:2:3 w l t 'Ux'
#splot 'solution_xy.dat' u 1:2:4 w l t 'Uy'
#splot 'solution_xy.dat' u 1:2:5 w l t 'Uz'
#splot 'solution_xy.dat' u 1:2:(0.5*sqrt($7*$7+$8*$8+$9*$9)) w l t 'Vorticity'
#splot 'solution_xy.dat' u 1:2:9 w l t 'ROT_Z'
splot 'outputs/solution_xy.dat' u 1:2:10 w l t 'Q-criterion'
#splot 'solution_xy.dat' u 1:2:11 w l t 'DIV(U)'
#splot 'solution_xy.dat' u 1:2:6 w l t 'P'
#splot 'solution_xz.dat' u 1:2:9 w l t ''
#splot 'solution_yz.dat' u 2:1:8 w l t 'rot_x'
#splot 'solution_inflow.dat' u 2:1:3

