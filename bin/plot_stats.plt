#!/usr/bin/gnuplot -persist
#
#    
#    	G N U P L O T
#    	Version 5.4 patchlevel 2    last modified 2021-06-01 
#    
#    	Copyright (C) 1986-1993, 1998, 2004, 2007-2021
#    	Thomas Williams, Colin Kelley and many others
#    
#    	gnuplot home:     http://www.gnuplot.info
#    	faq, bugs, etc:   type "help FAQ"
#    	immediate help:   type "help"  (plot window: hit 'h')
reset
# set terminal qt 0 font "Sans,9"
# set output
set grid xtics nomxtics ytics nomytics noztics nomztics nortics nomrtics \
 nox2tics nomx2tics noy2tics nomy2tics nocbtics nomcbtics
set grid layerdefault   lt 0 linecolor 0 linewidth 0.500,  lt 0 linecolor 0 linewidth 0.500
unset raxis
set key notitle
set key fixed right top vertical Right noreverse enhanced autotitle box
set key noinvert samplen 4 spacing 1 width 1 height 0.5
set key maxcolumns 0 maxrows 0
set key noopaque
set y2tics
set xlabel "t" 
set xlabel  font "" textcolor lt -1 norotate
set xrange [ * : * ] noreverse writeback
set ylabel "{/Symbol e}_t" 
set ylabel  font "" textcolor lt -1 rotate
set y2label "<U_i>^2" 
set y2label  font "" textcolor lt -1 rotate
set yrange [ 0 : 0.016 ] noreverse writeback
set y2range [ 0 : 0.16 ] noreverse writeback
set cblabel "" 
set cblabel  font "" textcolor lt -1 rotate
set cbrange [ * : * ] noreverse writeback
set locale "fr_FR.UTF-8"
plot "outputs/stats.dat" u 1:3 w l dt 5 lw 3 lt 2 t '{/Symbol e}_t osinco3d', "outputs/stats.dat" u 1:6 w l dt 5 lw 3 lt 4 t '<U_x>^2 osinco3d' axes x1y2, "outputs/stats.dat" u 1:8 w l dt 5 lw 3 lt 6 t '<U_z>^2 osinco3d' axes  x1y2
#    EOF
