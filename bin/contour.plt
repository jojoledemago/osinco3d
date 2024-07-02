#!/bin/gnuplot
# Visualization of contour plot

# Set output terminal to JPEG with specified size
set terminal jpeg size 1620,1080
# Set the output file name
set output 'tampon.jpeg'
# Set the view to 2D map projection
set view map
# Disable surface drawing
unset surface
# Set the plot to be square
set size ratio -1
# Set tics to be drawn outside the plot border
set tics out
# Enable contour drawing
set contour
# Set contour levels to be drawn at increments
set cntrparam levels incremental 0, 0.1, 1
# Set the number of samples for iso lines
set isosamples 150
# Display the key (legend) outside the plot with a box
set key outside box
# Unset tics for x, y, and z axes to remove axis labels
unset xtics
unset ytics
unset ztics
# Plot the data from 'tampon.dat' as lines with the title 'X-velocity'
splot 'tampon.dat' with lines title 'ROT_Z'
# Exit gnuplot
quit

