#!/bin/gnuplot
reset
unset multiplot
# visualization
#
#set terminal jpeg size 1080,720 # Sortie jpeg
#set output 'f_plot.png' # Nom du fichier de sortie
set nokey # Sans titre
set pm3d map # Trace d'une carte
#set palette defined
set palette rgbformulae 7,5,15 # Degrade bleu-rouge
#set palette defined ( 0 0 0 0, 1 1 1 1)
set size ratio -1 # Echelles orthonormees
#set xrange[0:pi]
#set yrange[0:pi]
set grid
#set view 0,0 # Vue de dessus
#set cbrange[0.5:1] # Etendue de la coloration
set tics out # Tics vers l'exterieur
#unset xtics # Pas de graduation selon z
#unset ytics
unset ztics
splot 'outputs/solution_xy.dat' u 1:2:(sqrt($7*$7+$8*$8+$9*$9))
#splot 'outputs/solution_xy.dat' u 1:2:4
#splot 'outputs/solution_xz.dat' u 1:2:(sqrt($7*$7+$8*$8+$9*$9))
#splot 'outputs/solution_xz.dat' u 1:2:5
#splot 'outputs/solution_yz.dat' u 2:1:(sqrt($7*$7+$8*$8+$9*$9))
#splot 'outputs/solution_yz.dat' u 2:1:12
