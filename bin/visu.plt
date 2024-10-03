#!/bin/gnuplot
# visualization
#
set terminal jpeg size 1620,1080 # Sortie jpeg
set output 'tampon.jpeg' # Nom du fichier de sortie
set nokey # Sans titre
set pm3d map # Trace d'une carte
set palette model RGB defined (-1 'black', 0 'red', 1 'yellow', 2 'white') #heat scale with black
set size ratio -1 # Echelles orthonormées
set tics out # Tics vers l'extérieur
#set cbrange[0:1] # Etendue de la coloration
unset xtics # Pas de graduation selon z
unset ytics
unset ztics
unset colorbox # Supprimer l'échelle de coloration
splot 'tampon.dat'

#set palette defined
#set palette defined ( 0 0 0 0, 1 1 1 1)
#set palette rgbformulae 7,5,15 # Degrade noir-jaune
#set palette model RGB defined (-1 'black', 0 'red', 1 'yellow', 2 'white') #heat scale with black
#set palette model RGB defined (0 'white', 1 'yellow', 2 'red')
#set palette model RGB defined (0 'white', 1 'yellow', 2 'red') #heat scale
#set palette model RGB defined (0 'red', 1 'yellow', 2 'whitr') #heat scale
#set palette model RGB defined (-1 'black', 0 'red', 1 'yellow', 2 'white') #heat scale with black
#set palette rgbformulae 33,13,10 # Degrade bleu-rouge
#set cbrange[-1.5:0] # Etendue de la coloration
#set xrange[0:16] # Intervalle en x
#set yrange[0:0.2] # Intervalle en y
